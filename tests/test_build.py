"""The build chain and the distribution's metadata.

The kernel is SPARTN's integer-to-value conversion, Vdecoded = Vencoded *
Res + Rngmin (SPARTN v2.0.3 section 6.3). Python evaluates the same
expression in IEEE 754 double with two roundings, so the kernel must agree
bit for bit; it does not if the compiler fuses the multiply and add into
one FMA (one rounding), which is what -ffp-contract=off prevents.

The licence is MIT in LICENSE and in the package metadata; LICENSE names
its holder, Trace Navigation LLC. THIRD_PARTY_NOTICES.md, the attribution of
the published values the package holds, is the second licence file, so it is
in the source distribution and in the wheel's metadata. RELEASES.md has the form of the
sibling's release notes and names the package's version first.
"""

import re
import struct
import tomllib
from pathlib import Path

from ssrcore._kernels.scale import decode

ROOT = Path(__file__).resolve().parents[1]
PROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]

# SF020, satellite corrections: 14 bits, range ±16.382 m, resolution 0.002 m (Table 6.2).
SF020_BITS, SF020_RES, SF020_MIN = 14, 0.002, -16.382


def bits(x: float) -> bytes:
    return struct.pack("<d", x)


def test_sf020_range_ends_and_middle():
    assert decode(0, SF020_RES, SF020_MIN) == SF020_MIN
    assert abs(decode(8191, SF020_RES, SF020_MIN)) < 1e-12
    assert abs(decode(16382, SF020_RES, SF020_MIN) - 16.382) < 1e-12


def test_kernel_matches_ieee_double_bit_for_bit():
    for encoded in range(1 << SF020_BITS):
        expected = encoded * SF020_RES + SF020_MIN
        assert bits(decode(encoded, SF020_RES, SF020_MIN)) == bits(expected), encoded


def test_licence_is_mit_held_by_trace_navigation():
    licence = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert licence.startswith("MIT License\n\nCopyright (c) 2026 Trace Navigation LLC\n")
    assert PROJECT["license"] == "MIT"
    assert PROJECT["license-files"][0] == "LICENSE"


def test_third_party_notices_are_a_licence_file():
    assert PROJECT["license-files"] == ["LICENSE", "THIRD_PARTY_NOTICES.md"]
    notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    assert notices.startswith("# Third-party notices\n")


RELEASE_HEADING = re.compile(r"Version \d+\.\d+\.\d+ \(\d{4}-(\d\d|XX)-(\d\d|XX)\)")


def _version(text: str) -> tuple[int, ...]:
    """[design: X.Y.Z as integers, so that versions compare as numbers]"""
    return tuple(int(part) for part in text.split("."))


def test_release_notes_have_their_form():
    # Per release, newest first: the heading, a line of "=" as long, then
    # "* " items with "  * " sub-items; nothing else.
    lines = (ROOT / "RELEASES.md").read_text(encoding="utf-8").splitlines()
    headings = [i for i, line in enumerate(lines) if line.startswith("Version ")]
    assert headings and headings[0] == 0
    for start, end in zip(headings, [*headings[1:], len(lines)]):
        heading, underline, *items = lines[start:end]
        assert RELEASE_HEADING.fullmatch(heading), heading
        assert underline == "=" * len(heading), heading
        assert items and items[0].startswith("* "), heading
        for item in items:
            assert re.fullmatch(r"(  )?\* \S.*", item), item
    versions = [_version(lines[i].split()[1]) for i in headings]
    assert versions == sorted(versions, reverse=True) and len(set(versions)) == len(versions)


def test_release_notes_name_the_current_version_first():
    first = (ROOT / "RELEASES.md").read_text(encoding="utf-8").splitlines()[0]
    assert first.split()[1] == PROJECT["version"], first
