"""The kernels on microcontrollers, and one query in C (DESIGN 4, Q17).

- Every kernel compiles unchanged, as C99 with the strict flags of the host
  build (``tools/portability.py`` ``STRICT``), for each target whose
  compiler is on this machine: the host, Arm Cortex-M4 with the flags of
  DESIGN 4, ESP32-S3, ESP32, ESP32-C3 and ESP32-C6. A target without its
  compiler skips, naming it. For each ELF target the objects hold no mutable
  data, call no other kernel, and need nothing but the math library and the
  compiler's run-time routines; the size and routine table is printed
  (``uv run pytest tests/test_portability.py -rP`` shows it).
- The compiler's own include trace shows that a kernel opens no header but
  ``<stdint.h>``, ``<math.h>`` and kernel headers.
- ``tools/query`` runs one query's kernel calls end to end in plain C with
  the inputs of the public test
  ``test_engine.test_every_term_of_a_gps_record_is_the_toolkits_value``;
  every input is the engine's decoded input, and every value equals the
  engine's record to the last bit.
- The counting build (``tools/query/count.hpp``) gives the same values and
  counts the double operations of each part of the query.
- The query and the kernels build for the ESP32-S3 into a static archive and
  a relocatable object with the math library linked in.
"""

from __future__ import annotations

import math
import struct
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import portability  # noqa: E402

import test_engine as te  # noqa: E402
from ssrcore.adapters import spartn_atmosphere, spartn_ocb  # noqa: E402
from ssrcore.model.types import System  # noqa: E402
from ssrcore.toolkit import astro, geodesy  # noqa: E402

TARGETS = {t.name: t for t in portability.targets()}


def same(a: float, b: float) -> bool:
    """Bit equality of two doubles (0.0 and -0.0 differ)."""
    return struct.pack("<d", a) == struct.pack("<d", b)


def target_or_skip(name: str) -> portability.Target:
    t = TARGETS[name]
    if t.compiler is None:
        pytest.skip(t.absent)
    return t


# ---- every kernel, every target -------------------------------------------------------


@pytest.mark.parametrize("name", list(TARGETS))
def test_every_kernel_compiles_unchanged(name, tmp_path):
    """DESIGN 3.1 and 4: the kernels as written, C99, strict flags, for the target."""
    t = target_or_skip(name)
    failed, lines = [], [f"{t.name}: {t.what}", f"  {'kernel':<12}{'code':>8}{'const':>8}"]
    needed: set[str] = set()
    total = [0, 0]
    for k in portability.kernels():
        obj = tmp_path / f"{k}.o"
        r = portability.compile_kernel(t, k, obj)
        if r.returncode != 0:
            failed.append(f"{k}:\n{r.stderr}")
            continue
        if not t.elf:
            continue
        f = portability.footprint(t, obj)
        assert f.mutable == 0, f"{k} holds mutable data on {name} (DESIGN 3.1: no static mutable state)"
        assert not [s for s in f.undefined if s.startswith("ssr_")], f"{k} calls another kernel"
        needed |= set(f.undefined)
        total[0] += f.code
        total[1] += f.constant
        lines.append(f"  {k:<12}{f.code:>8}{f.constant:>8}")
    assert not failed, "\n".join(failed)
    if t.elf:
        math_, soft, runtime, other = portability.classify(needed)
        assert not other, f"symbols outside the math library and the compiler's run-time: {other}"
        lines += [f"  {'total':<12}{total[0]:>8}{total[1]:>8}",
                  f"  math library: {' '.join(math_)}",
                  f"  software double: {' '.join(soft) or 'none'}",
                  f"  other run-time: {' '.join(runtime) or 'none'}"]
        print("\n".join(lines))


def test_kernels_open_only_the_allowed_headers():
    """DESIGN 3.1, from the preprocessor's own trace (-H) rather than the source text:
    a kernel source and a kernel header include only <stdint.h>, <math.h> and kernel
    headers."""
    t = target_or_skip("host")
    include = str(portability.KERNEL_INCLUDE)
    wrong = []
    for k in portability.kernels():
        pairs = portability.included_headers(t, k)
        assert pairs, f"no include trace for {k}"
        for includer, included in pairs:
            ours = includer.endswith(f"{k}.c") or includer.startswith(include)
            allowed = included.startswith(include) or Path(included).name in ("stdint.h", "math.h")
            if ours and not allowed:
                wrong.append(f"{k}: {includer} includes {included}")
    assert not wrong, "\n".join(wrong)


# ---- one query in C ---------------------------------------------------------------------


@pytest.fixture(scope="module")
def host_query(tmp_path_factory):
    target_or_skip("host")
    directory = tmp_path_factory.mktemp("query")
    exe = portability.build_query(directory)
    return exe, portability.query_values(exe)


@pytest.fixture(scope="module")
def engine_record():
    e = te.engine_with()
    rec = e.query(te.G05, te.T, te.RECEIVER)
    assert rec.missing == {}
    return e, rec


def test_the_c_programs_inputs_are_the_engines(host_query, engine_record):
    """tools/query/scene.c holds exactly what the adapters gave the engine."""
    _, c = host_query
    e, _ = engine_record
    eph = e.ephemerides.issues(te.G05)[0].ephemeris
    block = spartn_ocb.ocb_input(te.decoded(te.ocb(System.GPS, te.TS, [te.GPS_SAT], eos=True))).satellites[0]
    hp = spartn_atmosphere.from_hpac(te.decoded(te.hpac(System.GPS, te.TS, (5,)))).areas[0]
    iono = hp.ionosphere.satellites[0]
    expected = {f"in.receiver[{i}]": te.RECEIVER[i] for i in range(3)}
    expected.update({f"in.{n}": getattr(eph, n) for n in (
        "sqrt_a", "e", "m0", "delta_n", "omega", "omega0", "omega_dot", "i0", "idot",
        "cuc", "cus", "crc", "crs", "cic", "cis", "af0", "af1", "af2")})
    expected.update({"in.toe": eph.toe.seconds, "in.toc": eph.toc.seconds, "in.t_rx.week": te.T.week,
                     "in.t_rx.second": te.T.second, "in.t_rx.fraction": te.T.fraction,
                     "in.yaw_deg": block.orbit.yaw, "in.d_c": block.clock.correction,
                     "in.hydrostatic": hp.troposphere.hydrostatic})
    expected.update({f"in.d_o[{i}]": v for i, v in
                     enumerate((block.orbit.radial, block.orbit.along, block.orbit.cross))})
    expected.update({f"in.tropo_c[{i}]": v for i, v in enumerate(hp.troposphere.coefficients)})
    expected.update({f"in.tropo_residuals[{i}]": v for i, v in enumerate(hp.troposphere.residuals)})
    expected.update({f"in.iono_c[{i}]": v for i, v in enumerate(iono.coefficients)})
    expected.update({f"in.iono_residuals[{i}]": v for i, v in enumerate(iono.residuals)})
    wrong = [f"{k}: C {c[k]!r}, engine {v!r}" for k, v in expected.items() if not same(c[k], float(v))]
    assert not wrong, "\n".join(wrong)


def test_the_c_program_gives_the_engines_record_to_the_last_bit(host_query, engine_record):
    """tools/query: the kernel calls of one query in C, with the engine's arithmetic between
    them, give every value of the record bit for bit (same compiler flags, same libm)."""
    _, c = host_query
    _, rec = engine_record
    s, g, m, tr, io = rec.satellite, rec.geometry, rec.modelled, rec.troposphere, rec.ionosphere
    sky = astro.sky(te.T.week, te.T.seconds, -float(te.LEAP))
    expected = {
        "query.t_tx.week": rec.t_tx.week, "query.t_tx.second": rec.t_tx.second,
        "query.t_tx.fraction": rec.t_tx.fraction, "query.tau": te.T - rec.t_tx,
        "query.broadcast_clock": s.broadcast.clock, "query.clock": s.clock, "query.clock_drift": s.clock_drift,
        "query.range": g.range, "query.range_rate": g.range_rate, "query.elevation": g.elevation,
        "query.azimuth": g.azimuth, "query.sagnac": m.earth_rotation, "query.shapiro": m.shapiro,
        "query.relativity_clock": m.relativity_clock, "query.solid_tide": m.solid_tide,
        "query.windup": m.windup, "query.ure": rec.clock_correction.ure,
        "query.hydrostatic_mapping": tr.hydrostatic_mapping, "query.wet_mapping": tr.wet_mapping,
        "query.tropo_slant": tr.slant, "epoch.hydrostatic_zenith": tr.hydrostatic_zenith,
        "epoch.wet_zenith": tr.wet_zenith, "epoch.tropo_polynomial": tr.polynomial, "epoch.tropo_grid": tr.grid,
        "epoch.tropo_sigma": tr.sigma, "query.iono_polynomial": io.polynomial, "query.iono_grid": io.grid,
        "query.slant_tec": io.slant_tec, "epoch.iono_sigma": io.sigma,
        "epoch.tt_centuries": sky.tt_centuries, "epoch.gmst": sky.gmst,
    }
    for name, vector in (("query.broadcast_position", s.broadcast.position), ("query.velocity", s.velocity),
                         ("query.position", s.position), ("query.orbit_ecef", rec.orbit.ecef),
                         ("query.unit", g.unit), ("epoch.tide_displacement", m.tide_displacement),
                         ("epoch.sun", sky.sun), ("epoch.moon", sky.moon),
                         ("receiver.geodetic", geodesy.ecef_to_geodetic(te.RECEIVER))):
        expected.update({f"{name}[{i}]": v for i, v in enumerate(vector)})
    for k, signal in enumerate(rec.signals):
        expected[f"query.code_delay[{k}]"] = io.code_delay[signal]
        expected[f"query.phase_delay[{k}]"] = io.phase_delay(signal)
    wrong = [f"{k}: C {c[k]!r}, engine {v!r}" for k, v in expected.items() if not same(c[k], float(v))]
    assert not wrong, "\n".join(wrong)
    assert c["query.iterations"] == 3.0
    assert math.degrees(c["query.elevation"]) == pytest.approx(70.7, abs=0.2)


def test_the_c_program_times_its_parts(host_query):
    """main.c's ``time N`` (the figures of DESIGN 2.8 come from a larger N)."""
    exe, _ = host_query
    t = portability.time_query(exe, 200)
    assert 0.0 < t["query_seconds"] < 1e-3 and 0.0 < t["epoch_seconds"] < 1e-2
    print(f"host: {t['query_seconds'] * 1e6:.2f} us per query, {t['epoch_seconds'] * 1e6:.2f} us per epoch")


# ---- the operation counts ---------------------------------------------------------------


def test_the_counting_build_counts_and_agrees_bit_for_bit(host_query, tmp_path):
    """count.hpp: the same values as the C build, and the operations of each part."""
    if portability.cxx_compiler() is None:
        pytest.skip("no C++ compiler for the counting build")
    _, c = host_query
    n, values = portability.counts(portability.build_count(tmp_path))
    wrong = [f"{k}: counted {v!r}, C {c[k]!r}" for k, v in values.items() if not same(v, c[k])]
    assert values and not wrong, "\n".join(wrong)
    units = {(part, unit) for part, unit, _ in n}
    # the per-epoch work (DESIGN 2.8) is not in the query part, and the query part is not empty
    assert ("query", "astro") not in units and ("query", "tides") not in units
    assert ("epoch", "astro") in units and ("epoch", "tides") in units
    per_query: dict[str, int] = {}
    for (part, _, kind), v in n.items():
        if part == "query":
            per_query[kind] = per_query.get(kind, 0) + v
    assert per_query.get("mul", 0) > 100 and per_query.get("sqrt", 0) > 0
    print(portability.count_table(n))


# ---- the firmware artefact --------------------------------------------------------------


def test_the_query_builds_for_the_esp32s3(tmp_path):
    """tools/query and the kernels as a static archive and a relocatable object for the
    ESP32-S3, the math library and the compiler's routines linked in; nothing left
    undefined but what the C library of a firmware provides."""
    t = target_or_skip("esp32s3")
    a = portability.build_archive(t, tmp_path)
    assert a["archive"].is_file() and a["linked"].is_file()
    whole = a["whole"]
    assert whole.mutable == 0
    assert set(whole.undefined) <= {"__errno", "memcpy"}, whole.undefined
    rt = a["routines"]
    print(f"esp32s3 query and kernels with libm and libgcc: code {whole.code}, constant {whole.constant} bytes; "
          f"routines alone code {rt.code}, constant {rt.constant}; undefined {' '.join(whole.undefined)}")
