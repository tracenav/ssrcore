"""Broadcast orbits and clocks: Keplerian (GPS, Galileo, BeiDou, QZSS) and GLONASS.

DESIGN.md section 2.7. Typed wrappers of the kernels ``kepler`` and
``glonass``; no state, no arithmetic of their own.

Units are SI and angles radians (DESIGN 4): broadcast angles in semicircles
are converted with :func:`semicircles`, which uses the value of pi the ICDs
fix for orbit computation. Times are (week, seconds of week); the time of
evaluation and the ephemeris reference time must be in the same time scale
and week count (GPS time for GPS and QZSS, GST for Galileo, BDT for BeiDou,
and for GLONASS whatever scale the caller brought tb to, DESIGN T4).
Validity in time is not judged here; it is the model's (DESIGN 2.5).
"""

from __future__ import annotations

from dataclasses import dataclass

from ssrcore._kernels import glonass as _glonass
from ssrcore._kernels import kepler as _kepler
from ssrcore.toolkit import KernelError, Vector
from ssrcore.toolkit.timescale import System

#: System codes of :class:`KeplerOrbit` (they select the ICD's mu and omega_e):
#: the members of :class:`~ssrcore.toolkit.timescale.System`. QZSS uses the
#: GPS constants of IS-GPS-200N with IRN-IS-200N-003 and -004 until IS-QZSS-PNT-003 is staged
#: (DESIGN Q14).
GPS = System.GPS
GALILEO = System.GALILEO
BEIDOU = System.BEIDOU
QZSS = System.QZSS

#: Largest GLONASS integration step, s (RTCM 10403.3 3.5.13).
GLONASS_MAX_STEP = _glonass.MAX_STEP
#: Largest |t - tb| the GLONASS integration accepts, s.
GLONASS_MAX_SPAN = _glonass.MAX_SPAN

_KEPLER_STATUS = {
    _kepler.BAD_ELEMENTS: "bad elements (unknown system, e outside [0, 1) or sqrt(A) not positive)",
    _kepler.NOT_CONVERGED: "Kepler's equation did not converge",
}
_GLONASS_STATUS = {_glonass.TOO_FAR: "|t - tb| beyond the integration bound"}


@dataclass(frozen=True, slots=True)
class KeplerOrbit:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-III; Galileo OS SIS ICD
    2.2 5.1.1 Table 67; BDS-SIS-ICD-B1I-3.0 5.2.4.12] Broadcast orbit elements in SI
    units, angles in radians; ``geo`` for a BeiDou GEO satellite."""

    system: int
    toe_week: int
    toe: float
    sqrt_a: float
    e: float
    m0: float
    delta_n: float
    omega: float
    omega0: float
    omega_dot: float
    i0: float
    idot: float
    cuc: float
    cus: float
    crc: float
    crs: float
    cic: float
    cis: float
    geo: bool = False


@dataclass(frozen=True, slots=True)
class KeplerClock:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.3;
    BDS-SIS-ICD-B1I-3.0 5.2.4.9] Clock polynomial: toc and af0 (s),
    af1 (s/s), af2 (s/s^2)."""

    toc_week: int
    toc: float
    af0: float
    af1: float
    af2: float


@dataclass(frozen=True, slots=True)
class GlonassOrbit:
    """[spec: GLONASS ICD 5.1 4.4] Broadcast state at tb in PZ-90: position
    (m), velocity (m/s), luni-solar acceleration (m/s^2), tau_n (s) and
    gamma_n; tb as (week, seconds of week) in the caller's time scale."""

    tb_week: int
    tb: float
    position: Vector
    velocity: Vector
    lunisolar: Vector
    tau_n: float = 0.0
    gamma_n: float = 0.0


@dataclass(frozen=True, slots=True)
class OrbitState:
    """[design: one value for a satellite state] Position (m) and velocity
    (m/s) in the system's Earth-fixed frame; for a Keplerian orbit also the
    eccentric anomaly (rad), which the relativistic clock term of
    IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1 uses."""

    position: Vector
    velocity: Vector
    eccentric_anomaly: float | None = None


def semicircles(value: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2; Galileo OS SIS ICD 2.2 Table
    68; BDS-SIS-ICD-B1I-3.0 Table 5-11] Radians from semicircles, with
    pi = 3.1415926535898."""
    return _kepler.semicircles(value)


def beidou_geo(prn: int) -> bool:
    """[spec: BDS-SIS-ICD-B1I-3.0 4.3 Table 4-1] Whether a BeiDou SV ID is a
    GEO satellite (1 to 5, 59 to 63)."""
    return _kepler.beidou_geo(prn)


def eccentric_anomaly(mean_anomaly: float, e: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV] Eccentric anomaly
    (rad) by the table's iteration."""
    status, ek = _kepler.eccentric_anomaly(mean_anomaly, e)
    if status != _kepler.OK:
        raise KernelError("kepler", status, _KEPLER_STATUS[status])
    return ek


def kepler_state(orbit: KeplerOrbit, week: int, sow: float) -> OrbitState:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV; Galileo OS SIS ICD
    2.2 5.1.1 Table 68; BDS-SIS-ICD-B1I-3.0 5.2.4.12 Table 5-11] Antenna phase
    centre position and velocity at (week, sow) of the orbit's time scale."""
    o = orbit
    status, p, v, ek = _kepler.state(
        o.system, o.geo, o.toe_week, o.toe, o.sqrt_a, o.e, o.m0, o.delta_n, o.omega,
        o.omega0, o.omega_dot, o.i0, o.idot, o.cuc, o.cus, o.crc, o.crs, o.cic, o.cis,
        week, sow)
    if status != _kepler.OK:
        raise KernelError("kepler", status, _KEPLER_STATUS[status])
    return OrbitState(p, v, ek)


def kepler_clock(clock: KeplerClock, week: int, sow: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4
    Eq. 14; BDS-SIS-ICD-B1I-3.0 5.2.4.9] Clock offset (s) af0 + af1 dt + af2 dt^2,
    without the relativistic term (DESIGN 2.4, ``toolkit.relativity``) and
    without group delays (DESIGN E4)."""
    c = clock
    return _kepler.clock(c.af0, c.af1, c.af2, c.toc_week, c.toc, week, sow)


def kepler_clock_drift(clock: KeplerClock, week: int, sow: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] Clock drift (s/s) af1 + 2 af2
    dt, without the relativistic term."""
    c = clock
    return _kepler.clock_drift(c.af1, c.af2, c.toc_week, c.toc, week, sow)


def glonass_state(orbit: GlonassOrbit, week: int, sow: float) -> OrbitState:
    """[spec: GLONASS ICD 5.1 A.3.1.2 with the corrections of RTCM 10403.3
    3.5.13; RK4 at steps of at most 30 s] Position and velocity in PZ-90 at
    (week, sow) of the scale of tb."""
    o = orbit
    status, p, v = _glonass.state(o.position, o.velocity, o.lunisolar, o.tb_week, o.tb,
                                  week, sow)
    if status != _glonass.OK:
        raise KernelError("glonass", status, _GLONASS_STATUS[status])
    return OrbitState(p, v)


def glonass_clock(orbit: GlonassOrbit, week: int, sow: float) -> float:
    """[spec: GLONASS ICD 5.1 3.3.3, 4.4] Satellite clock offset (s) from
    GLONASS time, -tau_n + gamma_n (t - tb); tau_c is the time scale's."""
    o = orbit
    return _glonass.clock(o.tau_n, o.gamma_n, o.tb_week, o.tb, week, sow)


def glonass_clock_drift(orbit: GlonassOrbit) -> float:
    """[spec: GLONASS ICD 5.1 3.3.3, 4.4] Clock drift (s/s): gamma_n, the
    time derivative of -tau_n + gamma_n (t - tb)."""
    return orbit.gamma_n
