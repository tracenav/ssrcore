"""Sun and Moon positions, for tides, yaw and phase wind-up.

DESIGN.md section 2.7 and Q11. The Moon is ELP 2000-82B (Chapront-Touzé and
Chapront), the Sun VSOP87 version D (Bretagnon and Francou), both truncated
as ``kernels/src/astro.c`` states; the terrestrial frame follows IERS
Conventions (2010) chapter 5. Times are a GPS week and second of week, and
UT1 − GPS in seconds, which the caller knows from the leap seconds and, if
it has it, UT1 − UTC.

:func:`sky` is the one function that composes them for an epoch: the
consumers of the Sun and Moon (``toolkit.tides`` and ``toolkit.windup``, and
the nominal yaw) take its :class:`Sky` as an argument, so the astronomy
is evaluated once per epoch and the same values reach every term.
"""

from dataclasses import dataclass

from ssrcore._kernels import astro as _astro
from ssrcore.toolkit import Vector


def tt_centuries(week: int, seconds: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4; TM-23 3.1.1; TN36 ch. 5 sections
    5.5.7, 5.7.1] TT in Julian centuries from J2000.0."""
    return _astro.tt_centuries(week, seconds)


def gmst(week: int, seconds: float, ut1_minus_gps: float) -> float:
    """[spec: TN36 ch. 5 eqs. (5.15), (5.32)] Greenwich mean sidereal time, rad."""
    return _astro.gmst(week, seconds, ut1_minus_gps)


def moon_ecliptic(t: float) -> Vector:
    """[spec: ELP 2000-82B sections 2 to 4, 8] Geocentric Moon: longitude and
    latitude on the mean ecliptic and equinox of date (rad), distance (m)."""
    return _astro.moon_ecliptic(t)


def sun_ecliptic(t: float) -> Vector:
    """[spec: VSOP87 version D] Geocentric Sun: longitude and latitude on the
    mean ecliptic and equinox of date (rad), distance (m)."""
    return _astro.sun_ecliptic(t)


def nutation(t: float) -> tuple[float, float]:
    """[spec: TN36 ch. 5 eq. (5.35), extract of Tables 5.3a and 5.3b]
    Nutation in longitude and obliquity, rad."""
    return _astro.nutation(t)


def true_equatorial(lbr: Vector, t: float) -> Vector:
    """[spec: TN36 ch. 5 section 5.4.5, eq. (5.40)] Ecliptic-of-date
    coordinates to the true equator and equinox of date, m."""
    dpsi, deps = _astro.nutation(t)
    return _astro.ecliptic_to_true_equator(lbr, dpsi, _astro.mean_obliquity(t), deps)


@dataclass(frozen=True, slots=True)
class Sky:
    """[design: one value for what the tide and wind-up need of an epoch] The
    geocentric Sun and Moon in ECEF (m), TT in Julian centuries from J2000.0
    and Greenwich mean sidereal time (rad)."""

    sun: Vector
    moon: Vector
    tt_centuries: float
    gmst: float


def sky(week: int, seconds: float, ut1_minus_gps: float) -> Sky:
    """[spec: ELP 2000-82B; VSOP87 version D; TN36 ch. 5 sections 5.4.3,
    5.4.5; IAU 1982 GMST with the TN36 5.5.7 terms] The Sun, the Moon, TT
    and GMST at a GPS time; ``ut1_minus_gps`` is UT1 − GPS in seconds."""
    sun, moon = _astro.sun_moon_ecef(week, seconds, ut1_minus_gps)
    return Sky(sun, moon, _astro.tt_centuries(week, seconds), _astro.gmst(week, seconds, ut1_minus_gps))
