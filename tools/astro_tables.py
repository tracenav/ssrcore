"""Generate the truncated ELP 2000-82B and VSOP87D tables of kernels/src/astro.c.

    uv run python tools/astro_tables.py refs/Q11_sun_moon_tides           # rewrite astro.c
    uv run python tools/astro_tables.py refs/Q11_sun_moon_tides --check   # compare only

The sources are the files of the authors as IMCCE distributes them (see
THIRD_PARTY_NOTICES.md and refs/README.md): ELP1 to ELP36 of ELP 2000-82B
and VSOP87D.ear of VSOP87. They are read here and never committed. The
selection is that of the
comment at the head of astro.c: every ELP term with |A| >= 0.2 (arcsec or
km), every VSOP87D Earth term with |A| >= 1e-6 (rad or au), in file order,
values copied verbatim (Python's shortest round-trip form of each double).

[design: DESIGN 2.7 `astro`; the generator is committed so that the tables
can be checked against their sources (tests/test_astro_tables.py)]
"""

from __future__ import annotations

import sys
from pathlib import Path

ELP_MIN = 0.2
VSOP_MIN = 1e-6
BEGIN = "/* BEGIN GENERATED TABLES */\n"
END = "/* END GENERATED TABLES */\n"
TARGET = Path(__file__).resolve().parent.parent / "kernels" / "src" / "astro.c"
MAIN_NAMES = {1: "ELP_MAIN_LONGITUDE", 2: "ELP_MAIN_LATITUDE", 3: "ELP_MAIN_DISTANCE"}


def read_elp(directory: Path) -> dict[int, list]:
    """The 36 ELP 2000-82B files, each a list of rows after its header line.

    [spec: ELP 2000-82B notice, file formats: ELP1-3 four 3-digit multipliers
    of D, l', l, F then the amplitude (and six derivatives, not used); ELP4-9
    and 22-36 five multipliers then phase (deg) and amplitude; ELP10-21
    eleven multipliers then phase and amplitude]
    """
    series = {}
    for k in range(1, 37):
        rows = []
        for line in (directory / f"ELP{k}").read_text().splitlines()[1:]:
            if not line.strip():
                continue
            width = 4 if k <= 3 else 5 if (k <= 9 or k >= 22) else 11
            ints = [int(line[i * 3:(i + 1) * 3]) for i in range(width)]
            values = line[width * 3:].split()
            if k <= 3:
                rows.append((ints, float(values[0])))
            else:
                rows.append((ints, float(values[0]), float(values[1])))
        series[k] = rows
    return series


def read_vsop(path: Path) -> dict[tuple[int, int], list[tuple[float, float, float]]]:
    """VSOP87D.ear: the terms A, B, C of each (coordinate, power of time) block.

    [spec: VSOP87 vsop87.doc, file format: a header line naming the
    coordinate (column 42) and the power of time (column 60), then one term
    per line with A, B, C in columns 80-97, 98-111, 112-131]
    """
    series: dict[tuple[int, int], list] = {}
    key = None
    for line in path.read_text().splitlines():
        if line.startswith(" VSOP87"):
            key = (int(line[41]), int(line[59]))
            series[key] = []
        elif line.strip():
            series[key].append((float(line[79:97]), float(line[97:111]), float(line[111:131])))
    return series


def _perturbation_layout(k: int) -> tuple[int, int]:
    """Coordinate (0 longitude, 1 latitude, 2 distance) and power of time of ELP file k.

    [spec: ELP 2000-82B notice: files 4-9 Earth figure, 10-21 planetary,
    22-36 tidal, relativistic and solar eccentricity; within each group the
    three coordinates in turn, and the files of t¹ (and t² for 34-36) after
    those of t⁰]
    """
    coordinate = (k - 1) % 3
    if 4 <= k <= 9:
        power = 0 if k <= 6 else 1
    elif 10 <= k <= 21:
        power = 0 if k in (10, 11, 12, 16, 17, 18) else 1
    elif 22 <= k <= 33:
        power = 0 if k in (22, 23, 24, 28, 29, 30, 31, 32, 33) else 1
    else:
        power = 2
    return coordinate, power


def _arguments(k: int, ints: list[int]) -> list[int]:
    """The multipliers of ELP file k as the kernel's 13 arguments
    (zeta, D, l', l, F, Me, V, T, Ma, J, S, U, N).

    [spec: ELP 2000-82B notice: files 10-15 give eight planetary multipliers
    then D, l, F; files 16-21 seven planetary multipliers then D, l', l, F;
    the others zeta, D, l', l, F]
    """
    n = [0] * 13
    if 10 <= k <= 15:
        n[5:13] = ints[0:8]
        n[1], n[3], n[4] = ints[8], ints[9], ints[10]
    elif 16 <= k <= 21:
        n[5:12] = ints[0:7]
        n[1], n[2], n[3], n[4] = ints[7], ints[8], ints[9], ints[10]
    else:
        n[0:5] = ints[0:5]
    return n


def generate(sources: Path) -> str:
    """The text between the markers of astro.c, from the files in ``sources``.

    [design: the selection stated in astro.c and in this module's docstring]
    """
    elp = read_elp(sources)
    lines = []
    for k in (1, 2, 3):
        rows = [(ints, a) for ints, a in elp[k] if abs(a) >= ELP_MIN]
        lines.append(f"/* ELP{k}: {len(rows)} of {len(elp[k])} terms. */")
        lines.append(f"static const elp_main_term {MAIN_NAMES[k]}[{len(rows)}] = {{")
        for ints, a in rows:
            lines.append("    {{{{{:3d}, {:3d}, {:3d}, {:3d}}}, {}}},".format(*ints, repr(a)))
        lines.append("};")
        lines.append("")

    perturbations = []
    for k in range(4, 37):
        coordinate, power = _perturbation_layout(k)
        for ints, phase, a in elp[k]:
            if abs(a) >= ELP_MIN:
                perturbations.append((k, coordinate, power, _arguments(k, ints), phase, a))
    lines.append(f"/* ELP4 to ELP36: {len(perturbations)} terms. */")
    lines.append(f"static const elp_perturbation_term ELP_PERTURBATIONS[{len(perturbations)}] = {{")
    for k, coordinate, power, n, phase, a in perturbations:
        lines.append("    {{{}, {}, {{{}}}, {}, {}}}, /* ELP{} */".format(
            coordinate, power, ", ".join(f"{x:2d}" for x in n), repr(phase), repr(a), k))
    lines.append("};")
    lines.append("")

    vsop = read_vsop(sources / "VSOP87D.ear")
    rows = [(coordinate - 1, alpha, a, b, c)
            for (coordinate, alpha), terms in vsop.items() for a, b, c in terms if abs(a) >= VSOP_MIN]
    total = sum(len(terms) for terms in vsop.values())
    lines.append(f"/* VSOP87D Earth: {len(rows)} of {total} terms. */")
    lines.append(f"static const vsop_term VSOP87D_EARTH[{len(rows)}] = {{")
    for coordinate, alpha, a, b, c in rows:
        lines.append(f"    {{{coordinate}, {alpha}, {a!r}, {b!r}, {c!r}}},")
    lines.append("};")
    return "\n".join(lines) + "\n"


def committed(target: Path = TARGET) -> str:
    """The text between the markers of the committed astro.c. [design: as ``generate``]"""
    text = target.read_text()
    return text.split(BEGIN, 1)[1].split(END, 1)[0]


def main(argv: list[str]) -> int:
    """[design: command line of the module docstring]"""
    if len(argv) < 2:
        print(__doc__)
        return 2
    tables = generate(Path(argv[1]))
    if "--check" in argv[2:]:
        same = tables == committed()
        print("astro.c tables match their sources" if same else "astro.c tables differ from their sources")
        return 0 if same else 1
    text = TARGET.read_text()
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    TARGET.write_text(head + BEGIN + tables + END + tail)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
