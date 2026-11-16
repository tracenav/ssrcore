"""State-space corrections: SPARTN field scaling, orbit (radial, along, cross) and clock.

DESIGN.md section 2.7. Typed wrappers over the kernels ``scale`` (SPARTN
v2.0.3 section 6.3 and Table 6.2) and ``correction`` (SPARTN v2.0.3 sections 8.7 and
8.8; RTCM 10403.3 3.5.13.5 and 3.5.13.6 for corrections with rates, rule C11); no
arithmetic and no state here.

Field values are in the unit Table 6.2 gives them (m, TECU, degrees, seconds,
km, or a count). [design: each field can then be checked against the table;
a formula that needs radians or metres converts where it uses the value]
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

from ssrcore._kernels import scale as _scale
from ssrcore._kernels import correction as _correction
from ssrcore.toolkit import KernelError, Vector


#: [spec: SPARTN v2.0.3 Table 8.0] speed of light, m/s.
SPEED_OF_LIGHT: float = _correction.SPEED_OF_LIGHT


class Status(IntEnum):
    """Outcome of a field conversion. [design: DESIGN 3.1 status codes, one per reason a value is absent]"""

    OK = _scale.OK
    """The value is the field's value."""
    SPECIAL = _scale.SPECIAL
    """A special value with no number (invalid, unknown); value is None."""
    OUTSIDE = _scale.OUTSIDE
    """The encoded value is wider than the field or beyond its range; value is None."""
    ABOVE = _scale.ABOVE
    """The field says "greater than" value (last step of SF024, SF042, SF055)."""
    NO_FIELD = _scale.NO_FIELD
    """The field is a code (flag, mask, enumeration) or unknown, not converted; value is None."""


@dataclass(frozen=True)
class FieldValue:
    """A converted SPARTN field. [design: the status says why a value is absent]"""

    status: Status
    value: float | None


def decode(encoded: int, resolution: float, range_min: float) -> float:
    """Vdecoded = Vencoded * Res + Rngmin, without special-value checks.

    [spec: SPARTN v2.0.3 section 6.3]
    """
    return _scale.decode(encoded, resolution, range_min)


def field(sf: int, encoded: int) -> FieldValue:
    """Value of field SF<sf> (``sf=20`` for SF020) from its encoded value.

    Ranged fields are checked for special values, then for the range, then
    converted by section 6.3; fields whose special-value column gives a
    quantity per code (SF015, SF022, SF024, SF042, SF055, SF070, SF077,
    SF078) are read from that list.

    [spec: SPARTN v2.0.3 section 6.3, Table 6.2]
    """
    status, value = _scale.field(sf, encoded)
    return FieldValue(Status(status), value)


def quantity(sf: int, encoded: int) -> FieldValue:
    """The quantity field SF<sf> carries: ``field`` plus the offset of the field's notes.

    SF043 is the area average hydrostatic delay offset from 2.3 m; SF045 and
    SF048 are the troposphere coefficient T00 offset by 0.252 m. Other
    fields are as ``field``.

    [spec: SPARTN v2.0.3 Table 6.2, notes of SF043, SF045, SF048]
    """
    status, value = _scale.quantity(sf, encoded)
    return FieldValue(Status(status), value)


def _checked(status: int) -> None:
    """[design: DESIGN 3.1, a kernel failure is raised at the toolkit]"""
    if status != _correction.OK:
        raise KernelError("correction", status, "orbit frame undefined: r_dot or r x r_dot is zero")


def orbit_frame(r: Vector, r_dot: Vector) -> tuple[Vector, Vector, Vector]:
    """(e_radial, e_along, e_cross) of the broadcast state r (m), r_dot (m/s), ECEF.

    e_along = r_dot/|r_dot|, e_cross = r x r_dot/|r x r_dot|,
    e_radial = e_along x e_cross. [spec: SPARTN v2.0.3 section 8.8]
    Raises :class:`~ssrcore.toolkit.KernelError` when r_dot or r x r_dot is
    zero. [design: the frame is undefined there, which no broadcast orbit
    gives; the caller reports the satellite's orbit as missing]
    """
    status, er, ea, ec = _correction.orbit_frame(r, r_dot)
    _checked(status)
    return er, ea, ec


def orbit_correction(r: Vector, r_dot: Vector, d_o: Vector) -> Vector:
    """dX = [e_radial e_along e_cross] . dO in ECEF, m, for dO = (radial, along, cross), m.

    [spec: SPARTN v2.0.3 section 8.8] Raises as ``orbit_frame``.
    """
    status, d_x = _correction.orbit(r, r_dot, d_o)
    _checked(status)
    return d_x


def corrected_position(r: Vector, r_dot: Vector, d_o: Vector) -> Vector:
    """Corrected satellite nominal antenna phase centre X = X_b - dX, ECEF, m.

    [spec: SPARTN v2.0.3 section 8.8] Raises as ``orbit_frame``.
    """
    status, x = _correction.position(r, r_dot, d_o)
    _checked(status)
    return x


def corrected_clock(t_b: float, d_c: float) -> float:
    """Corrected satellite clock t_s = t_b - dC/c, s, for broadcast clock t_b (s) and dC (m).

    The relativistic clock term of the same section is
    ``ssrcore.toolkit.relativity``'s and is kept apart in the record
    (DESIGN 2.4). [spec: SPARTN v2.0.3 section 8.7]
    """
    return _correction.clock(t_b, d_c)


def orbit_at(d_o: Vector, rate: Vector, dt: float) -> Vector:
    """Orbit correction at a time, δO(t) = δO + δȮ·(t − t0), m, for δO (m) and δȮ (m/s)
    as (radial, along, cross) at the reference time t0 and ``dt`` = t − t0 (s).

    [spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6); IGS SSR v1.00 4.1] Rule C11; SPARTN
    corrections carry no rates and are not evaluated (C3).
    """
    if len(d_o) != 3 or len(rate) != 3:
        raise ValueError("d_o and rate are (radial, along, cross)")
    return _correction.orbit_at(d_o, rate, dt)


def clock_at(c0: float, c1: float, c2: float, dt: float) -> float:
    """Clock correction at a time, δC(t) = C0 + C1·(t − t0) + C2·(t − t0)², m, for C0 (m),
    C1 (m/s), C2 (m/s²) at the reference time t0 and ``dt`` = t − t0 (s).

    [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8); IGS SSR v1.00 4.2] Rule C11.
    """
    return _correction.clock_at(c0, c1, c2, dt)


def clock_rate_at(c1: float, c2: float, dt: float) -> float:
    """Rate of the clock correction, C1 + 2·C2·(t − t0), m/s.

    [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8), its derivative in t] Rule C11; the
    corrected clock drift is ``corrected_clock(drift, clock_rate_at(...))``.
    """
    return _correction.clock_rate_at(c1, c2, dt)
