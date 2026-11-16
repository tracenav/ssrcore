"""The ``scale`` kernel against SPARTN v2.0.3 section 6.3 and Table 6.2.

Every ranged field of Table 6.2 is transcribed here a second time, as the
table prints it (bits, range, resolution, special value), independently of the
kernel's table. The test derives each range end from that text and checks
that the kernel gives the lower end at code 0, the upper end at the code the
range implies, refuses the codes past it, and reports the special value.
"""

import math
import struct

import pytest

from ssrcore._kernels import scale as kernel
from ssrcore.toolkit import ssr
from ssrcore.toolkit.ssr import Status

# (SF, bits, range as printed, resolution, special value or None)
# [spec: SPARTN v2.0.3 Table 6.2]
TABLE_6_2 = [
    (5, 9, "0 to 511", 1, None),
    (18, 8, "0 to 255", 1, None),
    (19, 7, "0 to 127", 1, None),
    (20, 14, "±16.382", 0.002, None),
    (21, 6, "0 to 354", 6, 0x3F),
    (29, 11, "±20.46", 0.02, None),
    (30, 5, "1 to 32", 1, None),
    (31, 8, "0 to 255", 1, None),
    (32, 11, "±90.00", 0.1, None),
    (33, 12, "±180.00", 0.1, None),
    (34, 3, "1 to 8", 1, None),
    (35, 3, "1 to 8", 1, None),
    (36, 5, "0.1 to 3.2", 0.1, None),
    (37, 5, "0.1 to 3.2", 0.1, None),
    (39, 7, "0 to 127", 1, None),
    (43, 8, "±0.508", 0.004, None),
    (45, 7, "±0.252", 0.004, None),
    (46, 7, "±0.063", 0.001, None),
    (47, 9, "±0.0510", 0.0002, None),
    (48, 9, "±1.020", 0.004, None),
    (49, 9, "±0.255", 0.001, None),
    (50, 11, "±0.2046", 0.0002, None),
    (52, 6, "±0.124", 0.004, 0x3F),
    (53, 8, "±0.508", 0.004, 0xFF),
    (57, 12, "±81.88", 0.04, None),
    (58, 12, "±16.376", 0.008, None),
    (59, 13, "±8.190", 0.002, None),
    (60, 14, "±327.64", 0.04, None),
    (61, 14, "±65.528", 0.008, None),
    (62, 15, "±32.766", 0.002, None),
    (64, 4, "±0.28", 0.04, 0xF),
    (65, 7, "±2.52", 0.04, 0x7F),
    (66, 10, "±20.44", 0.04, 0x3FF),
    (67, 14, "±327.64", 0.04, 0x3FFF),
    (68, 4, "0 to 15", 1, None),
    (71, 2, "1 to 4", 1, None),
    (72, 2, "0 to 3", 1, None),
    (73, 8, "±85.00", 1.0, None),
    (74, 9, "±180.00", 1.0, None),
    (75, 4, "1 to 16", 1, None),
    (76, 4, "1 to 16", 1, None),
    (80, 12, "±511.75", 0.25, None),
    (82, 7, "±15.75", 0.25, 0x7F),
    (83, 11, "±255.75", 0.25, 0x7FF),
    (99, 10, "0 to 1023", 1, None),
    (100, 8, "0 to 255", 1, None),
    (101, 8, "0 to 255", 1, None),
]

# Geographic fields, whose range is not the symmetric span of their width.
GEOGRAPHIC = {32, 33, 73, 74}


def ends(text: str) -> tuple[float, float]:
    if text.startswith("±"):
        x = float(text[1:])
        return -x, x
    lo, hi = text.split(" to ")
    return float(lo), float(hi)


def top_code(text: str, res: float) -> int:
    lo, hi = ends(text)
    n = (hi - lo) / res
    assert abs(n - round(n)) < 1e-6, "range and resolution of the table disagree"
    return round(n)


def bits(x: float) -> bytes:
    return struct.pack("<d", x)


@pytest.mark.parametrize("sf, width, text, res, special", TABLE_6_2, ids=[f"SF{r[0]:03d}" for r in TABLE_6_2])
def test_range_ends(sf, width, text, res, special):
    lo, hi = ends(text)
    top = top_code(text, res)
    assert top < 1 << width

    low = ssr.field(sf, 0)
    assert low.status is Status.OK and low.value == lo
    high = ssr.field(sf, top)
    assert high.status is Status.OK
    assert high.value == pytest.approx(hi, rel=1e-12, abs=1e-12)

    for code in range(top + 1, 1 << width):
        expected = Status.SPECIAL if code == special else Status.OUTSIDE
        assert ssr.field(sf, code) == ssr.FieldValue(expected, None), code
    assert ssr.field(sf, 1 << width).status is Status.OUTSIDE
    if special is not None:
        assert ssr.field(sf, special).status is Status.SPECIAL


@pytest.mark.parametrize("sf, width, text, res, special", TABLE_6_2, ids=[f"SF{r[0]:03d}" for r in TABLE_6_2])
def test_every_code_is_section_6_3(sf, width, text, res, special):
    """Inside the range, the kernel's value is Vencoded * Res + Rngmin in IEEE double, bit for bit."""
    lo, _ = ends(text)
    for code in range(top_code(text, res) + 1):
        if code == special:
            continue
        assert bits(ssr.field(sf, code).value) == bits(code * res + lo), code


def test_symmetric_fields_span_their_width():
    """A cross-check of the transcription: every ±x field but the geographic
    ones uses all codes of its width, but the last, for the range."""
    for sf, width, text, res, special in TABLE_6_2:
        if text.startswith("±") and sf not in GEOGRAPHIC:
            assert top_code(text, res) == (1 << width) - 2, sf
            assert special in (None, (1 << width) - 1), sf


def test_satellite_correction_sf020_example():
    """6.3 worked by hand: code 8191 is the middle of ±16.382 m, so 0 m; one step is 2 mm."""
    assert ssr.field(20, 8191).value == pytest.approx(0.0, abs=1e-12)
    assert ssr.field(20, 8192).value == pytest.approx(0.002, abs=1e-12)
    assert ssr.decode(8192, 0.002, -16.382) == ssr.field(20, 8192).value


def test_yaw_degrees():
    assert ssr.field(21, 0).value == 0.0
    assert ssr.field(21, 59).value == 354.0
    assert ssr.field(21, 60).status is Status.OUTSIDE


CONTINUITY = [0.0, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0, 320.0]


@pytest.mark.parametrize("sf", [15, 22])
def test_continuity_seconds(sf):
    """[spec: Table 6.2 SF015, SF022] seconds by indicator."""
    assert [ssr.field(sf, k).value for k in range(8)] == CONTINUITY
    assert all(ssr.field(sf, k).status is Status.OK for k in range(8))
    assert ssr.field(sf, 8).status is Status.OUTSIDE


@pytest.mark.parametrize("sf, steps", [
    (24, [0.01, 0.02, 0.05, 0.1, 0.3, 1.0]),
    (42, [0.010, 0.020, 0.040, 0.080, 0.160, 0.320]),
    (55, [0.03, 0.05, 0.07, 0.14, 0.28, 0.56, 1.12, 2.24, 4.48, 8.96, 17.92, 35.84, 71.68, 143.36]),
])
def test_quality_lists(sf, steps):
    """[spec: Table 6.2 SF024, SF042, SF055] 0 unknown, steps, last code "greater than" the last step."""
    assert ssr.field(sf, 0) == ssr.FieldValue(Status.SPECIAL, None)
    for k, v in enumerate(steps, start=1):
        assert ssr.field(sf, k) == ssr.FieldValue(Status.OK, v)
    top = len(steps) + 1
    assert ssr.field(sf, top) == ssr.FieldValue(Status.ABOVE, steps[-1])
    assert ssr.field(sf, top + 1).status is Status.OUTSIDE


def test_quality_doubles_from_step_four_of_sf055():
    """An independent check of the SF055 transcription: from 0.07 TECU each step doubles (to rounding)."""
    values = [ssr.field(55, k).value for k in range(3, 15)]
    for a, b in zip(values, values[1:]):
        assert b == pytest.approx(2 * a, abs=0.0051)


def test_bpac_lists():
    """[spec: Table 6.2 SF070 (km), SF077, SF078 (degrees)]"""
    assert [ssr.field(70, k).value for k in range(4)] == [350.0, 400.0, 450.0, 500.0]
    for sf in (77, 78):
        assert [ssr.field(sf, k).value for k in range(4)] == [2.5, 5.0, 10.0, 15.0]


@pytest.mark.parametrize("sf", [8, 9, 10, 11, 13, 14, 16, 40, 41, 54, 84, 1, 200])
def test_codes_are_not_converted(sf):
    assert ssr.field(sf, 0) == ssr.FieldValue(Status.NO_FIELD, None)


def test_hydrostatic_delay_offset():
    """[spec: Table 6.2 SF043, "offset from 2.3m"] 2.3 m ± 0.508 m.

    Independent check of the offset's sign and size: the nominal dry zenith
    delay at sea level, alpha = 2.3 m of TM-23 5.4.2.2 eq. (5.66), is the
    middle code; at 2000 m height, alpha exp(-beta H) = 1.82 m with
    beta = 0.116e-3 /m, still inside the span.
    """
    assert ssr.quantity(43, 0).value == pytest.approx(1.792, abs=1e-12)
    assert ssr.quantity(43, 127).value == pytest.approx(2.3, abs=1e-12)
    assert ssr.quantity(43, 254).value == pytest.approx(2.808, abs=1e-12)
    assert ssr.quantity(43, 0).value < 2.3 * math.exp(-0.116e-3 * 2000.0) < ssr.quantity(43, 254).value
    assert ssr.quantity(43, 255).status is Status.OUTSIDE


def test_troposphere_t00_offset():
    """[spec: Table 6.2 SF045, SF048, "offset by 0.252m"; revision 1.6.2 note]

    The small T00 then spans 0 to 0.504 m: a wet zenith delay is not negative,
    and the nominal 0.1 m of TM-23 5.4.2.2 eq. (5.66) lies inside. Both sizes
    share the offset, so equal decoded values give equal T00.
    """
    assert 0.0 < 0.1 < ssr.quantity(45, 126).value
    assert ssr.quantity(45, 0).value == pytest.approx(0.0, abs=1e-12)
    assert ssr.quantity(45, 126).value == pytest.approx(0.504, abs=1e-12)
    assert ssr.quantity(48, 0).value == pytest.approx(-0.768, abs=1e-12)
    assert ssr.quantity(48, 510).value == pytest.approx(1.272, abs=1e-12)
    assert ssr.quantity(48, 255).value == pytest.approx(ssr.quantity(45, 63).value, abs=1e-12)


def test_quantity_equals_field_without_offset():
    for sf, width, text, res, special in TABLE_6_2:
        if sf in (43, 45, 48):
            continue
        for code in (0, top_code(text, res)):
            assert ssr.quantity(sf, code) == ssr.field(sf, code)
    assert ssr.quantity(52, 0x3F).status is Status.SPECIAL
    assert ssr.quantity(24, 7) == ssr.FieldValue(Status.ABOVE, 1.0)


def test_status_codes_of_the_binding():
    # DESIGN 3.1: 0 is success, every other condition a distinct positive code.
    assert (kernel.OK, kernel.SPECIAL, kernel.OUTSIDE, kernel.ABOVE, kernel.NO_FIELD) == (0, 1, 2, 3, 4)
