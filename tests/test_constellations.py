"""The per-constellation tables of DESIGN 2.2.1 against the code's own tables.

Each system's table in docs/DESIGN.md has one row per signal the code knows
for it: the signals of its SPARTN bias masks (``model.types.SPARTN_SIGNALS``),
of the IGS SSR identifiers the adapter reads (gnsswire's ``IGS_SIGNALS`` and
``adapters.igs_ssr._RTCM_SIGNALS``, through ``adapters.igs_ssr.signal``) and
of the virtual base station's MSM signal IDs (``vrs.msm.SIGNAL_ID``). A
signal added to the code without its row, a row without its signal, or a
cell that no longer says what the code does fails here, system by system and
signal by signal.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from gnsswire import rtcm3

from ssrcore.adapters import igs_ssr
from ssrcore.model.types import SPARTN_SIGNALS, Signal, System
from ssrcore.vrs import msm

DESIGN = Path(__file__).resolve().parents[1] / "docs" / "DESIGN.md"

# [design: the titles DESIGN 2.2.1 gives the tables, and the name of each system in
# gnsswire's IGS SSR tables (IGS SSR v1.00 Table 5)]
TITLES = {System.GPS: "GPS", System.GLONASS: "GLONASS", System.GALILEO: "Galileo",
          System.BEIDOU: "BeiDou", System.QZSS: "QZSS"}
IGS_NAMES = {System.GPS: "GPS", System.GLONASS: "GLONASS", System.GALILEO: "Galileo",
             System.BEIDOU: "BDS", System.QZSS: "QZSS"}
COLUMNS = ["Signal", "Band, MHz", "SPARTN", "IGS SSR", "Line-of-sight record", "Virtual base station"]
# [spec: SPARTN Table 6.2: the code and phase bias masks of each system]
CODE_MASK = {System.GPS: 27, System.GLONASS: 28, System.GALILEO: 105, System.BEIDOU: 106, System.QZSS: 107}
PHASE_MASK = {System.GPS: 25, System.GLONASS: 26, System.GALILEO: 102, System.BEIDOU: 103, System.QZSS: 104}


def section() -> str:
    """[design: DESIGN 2.2.1, from its heading to the next heading of its level or above]"""
    text = DESIGN.read_text(encoding="utf-8")
    start = text.index("#### 2.2.1 ")
    end = re.compile(r"^#{2,4} ", re.M).search(text, start + 1)
    return text[start:end.start() if end else len(text)]


def tables() -> dict[System, list[dict[str, str]]]:
    """[design: each table follows a line that is its system's title in bold; a row is
    a dict of its cells by column, the signal's backticks removed]"""
    out: dict[System, list[dict[str, str]]] = {}
    current = None
    for line in section().splitlines():
        title = re.fullmatch(r"\*\*(\w+)\*\*.*", line)
        if title:
            current = next((s for s, t in TITLES.items() if t == title[1]), None)
            continue
        if current is None or not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells == COLUMNS or set(line) <= set("|- "):
            continue
        assert len(cells) == len(COLUMNS), line
        row = dict(zip(COLUMNS, cells))
        row["Signal"] = row["Signal"].strip("`")
        out.setdefault(current, []).append(row)
    return out


def igs_identifiers(system: System) -> dict[str, int]:
    """[spec: IGS SSR v1.00 Table 7 IDF024 and RTCM 10403.3 DF380, DF382, as the adapter
    reads them (IG7)] the identifier of each RINEX code the adapter knows for a system."""
    ids = dict(rtcm3.IGS_SIGNALS[IGS_NAMES[system]])
    ids.update({i: code for (s, i), code in igs_ssr._RTCM_SIGNALS.items() if s == system})
    return {code: i for i, code in ids.items()}


def named(system: System, code: str) -> bool:
    """[spec: DESIGN 2.4 row λs] whether the code has a frequency for the signal's band."""
    try:
        Signal(system, code)
    except ValueError:
        return False
    return True


def bias_source(system: System, code: str) -> bool:
    """[design: DESIGN 2.2.1, a bias comes from a SPARTN mask position or an IGS SSR
    identifier the adapter takes]"""
    i = igs_identifiers(system).get(code)
    return code in SPARTN_SIGNALS[system] or (i is not None and igs_ssr.signal(system, i) is not None)


TABLES = tables()
SYSTEMS = list(System)


def rows(system: System) -> dict[str, dict[str, str]]:
    """[design: a system's rows by signal]"""
    found = TABLES.get(system, [])
    by_code = {r["Signal"]: r for r in found}
    assert len(by_code) == len(found), f"{system.name}: a signal has two rows"
    return by_code


def test_each_system_has_its_table_and_the_others_their_line():
    """[design: DESIGN 2.2.1, one table per system the engine names, and the line for SBAS
    and every other system]"""
    assert set(TABLES) == set(System)
    assert "**SBAS and every other system** are not supported" in section()


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_one_row_per_signal_the_code_knows(system):
    """[design: DESIGN 2.2.1, the rows are the union of the three tables of the code]"""
    known = set(SPARTN_SIGNALS[system]) | set(igs_identifiers(system)) | set(msm.SIGNAL_ID[system])
    documented = set(rows(system))
    assert documented == known, (f"{system.name}: rows without a signal {sorted(documented - known)}, "
                                 f"signals without a row {sorted(known - documented)}")


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_band_column_is_the_record_frequency(system):
    """[spec: Signal.frequency, DESIGN 2.4 row λs; GLONASS ICD 3.3.1.1 for f0 + K·Δf]"""
    for code, row in rows(system).items():
        cell = row["Band, MHz"]
        if not named(system, code):
            assert cell.endswith(": none"), (code, cell)
            continue
        if system == System.GLONASS:
            m = re.fullmatch(r"\w+, ([\d.]+) \+ ([\d.]+)·K", cell)
            assert m, (code, cell)
            f0, df = (round(float(x) * 1e6) for x in m.groups())
            assert Signal(system, code).frequency(0) == f0, (code, cell)
            assert Signal(system, code).frequency(1) - f0 == df, (code, cell)
        else:
            m = re.fullmatch(r"[\w+]+, ([\d.]+)", cell)
            assert m, (code, cell)
            assert round(float(m[1]) * 1e6) == Signal(system, code).frequency(), (code, cell)


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_spartn_column_is_the_bias_mask(system):
    """[spec: SPARTN Table 6.2 SF025 to SF028, SF102 to SF107, through SPARTN_SIGNALS]"""
    for code, row in rows(system).items():
        cell = row["SPARTN"]
        if code not in SPARTN_SIGNALS[system]:
            assert cell == "—", (code, cell)
            continue
        position = SPARTN_SIGNALS[system].index(code)
        assert cell == (f"bit {position}: C{code} (SF{CODE_MASK[system]:03d}), "
                        f"L{code} (SF{PHASE_MASK[system]:03d})"), (code, cell)
        assert Signal.from_spartn(system, position) == Signal(system, code)


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_igs_column_is_the_adapter_mapping(system):
    """[spec: IGS SSR v1.00 Table 7 IDF024; RTCM 10403.3 DF380, DF382; rule IG7 through
    adapters.igs_ssr.signal]"""
    ids = igs_identifiers(system)
    for code, row in rows(system).items():
        cell = row["IGS SSR"]
        if code not in ids:
            assert cell == "—", (code, cell)
            continue
        i = ids[code]
        m = re.fullmatch(r"ID (\d+)( \(RTCM DF38[02]\))?(: left out, no frequency)?", cell)
        assert m and int(m[1]) == i, (code, cell)
        assert bool(m[2]) == ((system, i) in igs_ssr._RTCM_SIGNALS), (code, cell)
        taken = igs_ssr.signal(system, i)
        assert bool(m[3]) == (taken is None), (code, cell)
        assert taken is None or taken == Signal(system, code), (code, cell)


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_igs_identifiers_left_out_are_not_in_a_table(system):
    """[spec: rule IG7: an identifier no table names is left out]"""
    ids = set(igs_identifiers(system).values())
    for i in range(32):
        if i not in ids:
            assert igs_ssr.signal(system, i) is None, (system.name, i)


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_record_column_follows_frequency_and_bias_source(system):
    """[design: DESIGN 2.2.1 and 2.4: a signal without a frequency cannot be named; one
    without a bias source has no b_s, B_s]"""
    for code, row in rows(system).items():
        cell = row["Line-of-sight record"]
        if not named(system, code):
            expected = "none: no frequency"
        elif bias_source(system, code):
            expected = "b_s, B_s, I_s, λs·w"
        else:
            expected = "I_s, λs·w; no bias"
        assert cell == expected, (code, cell)


@pytest.mark.parametrize("system", SYSTEMS, ids=lambda s: s.name)
def test_virtual_base_station_column_is_the_msm_signal_list(system):
    """[spec: RTCM 10403.3 Tables 3.5-91, 3.5-96, 3.5-99, 3.5-105, 3.5-108 through
    vrs.msm.SIGNAL_ID; design: DESIGN 2.6.3, a signal is sent with its code bias and
    frequency]"""
    for code, row in rows(system).items():
        cell = row["Virtual base station"]
        signal_id = msm.SIGNAL_ID[system].get(code)
        if signal_id is None:
            expected = "none: no MSM ID" if bias_source(system, code) else "—"
        elif not named(system, code):
            expected = f"ID {signal_id}: never sent, no frequency"
        elif not bias_source(system, code):
            expected = f"ID {signal_id}: never sent, no bias"
        else:
            expected = f"ID {signal_id}"
            assert msm.signal_id(Signal(system, code)) == signal_id
        assert cell == expected, (code, cell)


def test_a_signal_added_to_the_code_without_its_row_fails(monkeypatch):
    """[design: the check catches what it is for: a new MSM signal with no row]"""
    monkeypatch.setitem(msm.SIGNAL_ID, System.GPS, {**msm.SIGNAL_ID[System.GPS], "1Z": 99})
    with pytest.raises(AssertionError, match="signals without a row"):
        test_one_row_per_signal_the_code_knows(System.GPS)
