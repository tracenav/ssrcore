"""Carrier phase wind-up and the satellite antenna frame with yaw.

DESIGN.md section 2.7. The arithmetic is the C kernel ``kernels/src/windup.c``;
this module gives it typed names. Vectors are ECEF (x, y, z) tuples: positions
in metres, velocities in metres per second, directions as unit vectors.
Angles are in radians; wind-up is in cycles.

A typical use for one line of sight:

    x, y, _ = satellite_frame(r_sat, v_sat, yaw)       # SPARTN 8.14
    w = windup(r_sat, r_rcv, east, north, x, y)        # TM-23 5.5, in (-1/2, 1/2]
    w = unwrap(w_previous, w)                          # continuous along an arc

``east`` and ``north`` are the receiver's local unit vectors (the geodesy
toolkit gives them); a receiver whose antenna is turned passes its own
dipole directions instead.
"""

from collections.abc import Sequence

from ssrcore._kernels import windup as _k
from ssrcore.toolkit import KernelError, Vector


def _check(status: int, what: str) -> None:
    """[design: a degenerate geometry is a caller's error, not data] Raise on a kernel failure."""
    if status != _k.OK:
        raise KernelError("windup", status, f"{what}: degenerate geometry (a vector to normalise has zero length)")


def satellite_frame(position: Sequence[float], velocity: Sequence[float],
                    yaw: float = 0.0) -> tuple[Vector, Vector, Vector]:
    """[spec: SPARTN v2.0.3 section 8.14 and its figure; Table 6.2 SF021] Satellite antenna frame.

    ``position`` and ``velocity`` are the satellite's ECEF state; the velocity
    is made inertial with the Earth's rotation rate
    :data:`ssrcore.toolkit.geodesy.WGS84_OMEGA_E` (e_y uses the ECI velocity,
    8.14; one rate for every system, DESIGN 4). Returns the unit vectors (x, y, z): e_z towards the
    Earth's centre, e_y = e_z x v / |e_z x v|, e_x = e_y x e_z, with x and y
    turned by ``yaw`` (radians, right-handed about e_z; zero when SPARTN
    sends no yaw, 8.14).
    [design: the sense of the yaw rotation, which SPARTN does not state]
    """
    status, x, y, z = _k.frame_spartn(position, velocity, yaw)
    _check(status, "satellite_frame")
    return x, y, z


def satellite_frame_sun(position: Sequence[float], sun: Sequence[float]) -> tuple[Vector, Vector, Vector]:
    """[spec: TM-23/1 section 5.6.3, eqs. (5.76) to (5.78)] Nominal frame steered to the Sun.

    ``sun`` is the Sun's ECEF position (``toolkit.astro.sky(...).sun``). Returns the
    unit vectors (i, j, k): k towards the Earth's centre, j = unit(k x e) with
    e the direction to the Sun, i = j x k. This is the frame TM-23 5.5 uses
    for wind-up; SPARTN corrections use :func:`satellite_frame` instead.
    """
    status, i, j, k = _k.frame_sun(position, sun)
    _check(status, "satellite_frame_sun")
    return i, j, k


def windup(satellite: Sequence[float], receiver: Sequence[float],
           receiver_a: Sequence[float], receiver_b: Sequence[float],
           satellite_a: Sequence[float], satellite_b: Sequence[float]) -> float:
    """[spec: TM-23/1 section 5.5, eqs. (5.72) and (5.74), after Wu et al. 1993] Wind-up in cycles.

    The fractional wind-up in (-1/2, 1/2] for the line of sight from the
    satellite to the receiver (ECEF positions), with the receiver dipoles
    ``receiver_a``, ``receiver_b`` (east and north, TM-23 5.5) and the
    satellite dipoles ``satellite_a``, ``satellite_b`` (x and y of the
    satellite frame). Multiply by the wavelength for metres.
    """
    status, cycles = _k.cycles(satellite, receiver, receiver_a, receiver_b, satellite_a, satellite_b)
    _check(status, "windup")
    return cycles


def unwrap(previous: float, current: float) -> float:
    """[spec: TM-23/1 section 5.5, eqs. (5.71) and (5.73)] Wind-up continuous with ``previous``.

    Adds to ``current`` the whole number of cycles nearest to
    ``previous - current``. The first value of an arc is used as it is
    (N initialised as zero, TM-23 5.5); the caller keeps ``previous``.
    """
    return _k.unwrap(previous, current)
