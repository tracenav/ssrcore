"""Relativistic clock correction, Shapiro delay and Earth rotation during signal travel.

DESIGN.md section 2.7; SPARTN v2.0.3 sections 8.7 and 8.16; Ashby 2003
sections 2, 5 and 10. Positions are ECEF metres, velocities ECEF metres per
second. Each function is a typed wrapper of ``ssrcore._kernels.relativity``
and holds no state.

Signs follow the observation equations of DESIGN 2.4: ``clock`` is the
``dt_rel`` that is added to the satellite clock offset (so it enters a range
as ``-c dt_rel``); ``shapiro`` and ``sagnac`` are added to the geometric
range.
"""

from collections.abc import Sequence

from ssrcore._kernels import relativity as _k
from ssrcore.toolkit import KernelError

# c, GM and omega_E are those of ssrcore.toolkit.ssr.SPEED_OF_LIGHT and
# ssrcore.toolkit.geodesy.WGS84_GM, WGS84_OMEGA_E (DESIGN 4).


def clock(position: Sequence[float], velocity: Sequence[float]) -> float:
    """Relativistic clock correction -2 r.v / c^2 of a satellite, seconds.

    [spec: SPARTN v2.0.3 section 8.7; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1]
    ``position`` and ``velocity`` at the transmit time, ECEF or ECI alike. SPARTN 8.7
    applies it to every system but GLONASS, whose broadcast clock already
    contains it; that choice is the caller's.
    """
    px, py, pz = position
    vx, vy, vz = velocity
    return _k.clock(px, py, pz, vx, vy, vz)


def clock_kepler(mu: float, e: float, sqrt_a: float, eccentric_anomaly: float) -> float:
    """The same correction from broadcast elements, F e sqrt(A) sin E_k, seconds.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4 Eqs.
    15, 16; BDS-SIS-ICD-B1I-3.0 5.2.4.9] F = -2 sqrt(mu) / c^2 with the ``mu``
    (m^3/s^2) of the system's ICD; ``sqrt_a`` in m^1/2, ``eccentric_anomaly``
    in radians.
    """
    return _k.clock_kepler(mu, e, sqrt_a, eccentric_anomaly)


def shapiro(satellite: Sequence[float], receiver: Sequence[float]) -> float:
    """Shapiro path delay from satellite to receiver, metres.

    [spec: Ashby 2003 section 10, Eq. (96), second term]
    (2 GM / c^2) ln((r1 + r2 + l) / (r1 + r2 - l)). [design: the first term of
    Eq. (96), Phi_0 l / c^2, is not included: Ashby names the second term the
    Shapiro delay, and the first is the rate of coordinate time on the geoid
    applied to the whole path, not a propagation effect of the signal]

    Raises :class:`~ssrcore.toolkit.KernelError` when the straight path passes through the Earth's
    centre, where Eq. (96) is undefined.
    """
    sx, sy, sz = satellite
    rx, ry, rz = receiver
    status, delay = _k.shapiro(sx, sy, sz, rx, ry, rz)
    if status != _k.OK:
        raise KernelError("relativity", status, "path through the Earth's centre: Shapiro delay undefined")
    return delay


def sagnac(satellite: Sequence[float], receiver: Sequence[float]) -> float:
    """Sagnac term (Earth rotation during signal travel), metres.

    [spec: Ashby 2003 section 2] c tau - |r_s - r_r|, for the satellite in
    the ECEF frame of the transmit time and the receiver in that of the
    reception time, where c tau = |R_z(omega_E tau) r_s - r_r| is the
    distance in the frame of reception. Ashby's Eq. (7),
    omega_E (x_s y_r - y_s x_r) / c, is its first-order term; the exact form
    is computed [design: the remainder reaches 0.28 mm for a geostationary
    orbit, above 0.1 mm (DESIGN 2.4, row E)].
    Raises :class:`~ssrcore.toolkit.KernelError` when the two positions
    coincide.
    """
    sx, sy, sz = satellite
    rx, ry, rz = receiver
    status, delta = _k.sagnac(sx, sy, sz, rx, ry, rz)
    if status != _k.OK:
        raise KernelError("relativity", status, "satellite and receiver coincide: no path")
    return delta
