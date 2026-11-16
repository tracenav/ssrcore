"""WGS-84 geodesy: ECEF and geodetic coordinates, local frames, azimuth and elevation.

DESIGN.md section 2.7. Thin, typed wrappers of the C kernel
``kernels/src/geodesy.c`` (through ``ssrcore._kernels.geodesy``), which
implements ESA GNSS Data Processing Vol. I (TM-23) section 3.1.4 and its
Appendix B. Angles are in radians, lengths in metres; geodetic coordinates
are (latitude, longitude, ellipsoidal height). Every function takes the
ellipsoid as (a, f) and defaults to WGS-84.
"""

from __future__ import annotations

from ssrcore._kernels import geodesy as _k
from ssrcore.toolkit import KernelError, Vector

#: [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; TM-23 Table 3.2] WGS-84 semi-major
#: axis, m.
WGS84_A: float = _k.WGS84_A
#: [spec: TM-23 Table 3.2] WGS-84 flattening, 1/298.257223563.
WGS84_F: float = _k.WGS84_F
#: [spec: TM-23 Table 3.2; Ashby 2003 Eq. (13)] WGS-84 GM, m^3/s^2.
WGS84_GM: float = _k.WGS84_GM
#: [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; Ashby 2003 section 2] Earth
#: rotation rate, rad/s.
WGS84_OMEGA_E: float = _k.WGS84_OMEGA_E
#: [design: DESIGN 4] The double nearest pi.
PI: float = _k.PI


def eccentricity_squared(f: float = WGS84_F) -> float:
    """[spec: TM-23 B.3] Squared first eccentricity e² = 2f − f²."""
    return _k.eccentricity_squared(f)


def geodetic_to_ecef(geodetic: Vector, a: float = WGS84_A, f: float = WGS84_F) -> Vector:
    """[spec: TM-23 B.1.1] ECEF (x, y, z) of (latitude, longitude, height)."""
    lat, lon, h = geodetic
    return _k.geodetic_to_ecef(lat, lon, h, a, f)


def ecef_to_geodetic(ecef: Vector, a: float = WGS84_A, f: float = WGS84_F) -> Vector:
    """[spec: TM-23 B.1.2] (latitude, longitude, height) of an ECEF position.

    Longitude is in (−π, π], and 0 on the polar axis. Raises
    :class:`~ssrcore.toolkit.KernelError` where the iteration gives no finite solution (the
    centre of the ellipsoid, non-finite input).
    """
    x, y, z = ecef
    status, geodetic = _k.ecef_to_geodetic(x, y, z, a, f)
    if status != _k.OK:
        raise KernelError("geodesy", status, f"no geodetic coordinates for ECEF {ecef!r}")
    return geodetic


def ecef_to_enu(vector: Vector, latitude: float, longitude: float) -> Vector:
    """[spec: TM-23 B.2.2] (east, north, up) of an ECEF vector at (latitude, longitude)."""
    dx, dy, dz = vector
    return _k.ecef_to_enu(latitude, longitude, dx, dy, dz)


def enu_to_ecef(vector: Vector, latitude: float, longitude: float) -> Vector:
    """[spec: TM-23 B.2.1] ECEF components of an (east, north, up) vector at (latitude, longitude)."""
    e, n, u = vector
    return _k.enu_to_ecef(latitude, longitude, e, n, u)


def enu_basis(latitude: float, longitude: float) -> tuple[Vector, Vector, Vector]:
    """[spec: TM-23 B.2.1, eq. B.9] The unit vectors (ê, n̂, û) of the local frame, in ECEF.

    They are the columns of the matrix B.8, so each is the image of a unit
    ENU vector under :func:`enu_to_ecef`.
    """
    return (
        _k.enu_to_ecef(latitude, longitude, 1.0, 0.0, 0.0),
        _k.enu_to_ecef(latitude, longitude, 0.0, 1.0, 0.0),
        _k.enu_to_ecef(latitude, longitude, 0.0, 0.0, 1.0),
    )


def line_of_sight(
    satellite: Vector,
    receiver: Vector,
    receiver_geodetic: Vector | None = None,
    a: float = WGS84_A,
    f: float = WGS84_F,
) -> tuple[Vector, float, float]:
    """[spec: TM-23 B.3, eqs. B.13 and B.14] (unit vector, elevation, azimuth) from receiver to satellite.

    Both positions are ECEF. The unit vector points from the receiver to the
    satellite; elevation is in [−π/2, π/2] and azimuth in [0, 2π), clockwise
    from north, in the local frame at the receiver's geodetic latitude and
    longitude. ``receiver_geodetic`` may be given when the caller already
    has it; otherwise it is computed from ``receiver`` on (a, f). Raises
    :class:`~ssrcore.toolkit.KernelError` when the two positions coincide.
    """
    if receiver_geodetic is None:
        receiver_geodetic = ecef_to_geodetic(receiver, a, f)
    lat, lon, _ = receiver_geodetic
    status, unit, elevation, azimuth = _k.line_of_sight(satellite, receiver, lat, lon)
    if status != _k.OK:
        raise KernelError("geodesy", status, f"no line of sight from {receiver!r} to {satellite!r}")
    return unit, elevation, azimuth


def elevation_azimuth(
    satellite: Vector,
    receiver: Vector,
    receiver_geodetic: Vector | None = None,
    a: float = WGS84_A,
    f: float = WGS84_F,
) -> tuple[float, float]:
    """[spec: TM-23 B.3] (elevation, azimuth) of the satellite seen from the receiver; see :func:`line_of_sight`."""
    _, elevation, azimuth = line_of_sight(satellite, receiver, receiver_geodetic, a, f)
    return elevation, azimuth


def range_rate(satellite: Vector, velocity: Vector, receiver: Vector) -> tuple[float, float]:
    """[spec: TM-23 B.3, eq. B.13] (ρ, ρ̇): the distance from ``receiver`` to
    ``satellite`` (ECEF, m) and its rate (m/s) for the satellite velocity
    ``velocity`` and a receiver at rest in the frame.

    [design: the rate the query can give for free; a moving receiver adds
    −u·v_r] Raises :class:`~ssrcore.toolkit.KernelError` when the two
    positions coincide.
    """
    status, rho, rate = _k.range_rate(satellite, velocity, receiver)
    if status != _k.OK:
        raise KernelError("geodesy", status, f"no range from {receiver!r} to {satellite!r}")
    return rho, rate


def rotate_z(vector: Vector, angle: float) -> Vector:
    """[spec: Ashby 2003 section 2, R_z] The coordinates of ``vector`` in a frame
    turned by ``angle`` (rad) about z: the satellite position at transmission
    in the Earth-fixed frame of reception, for ``angle`` = ω_E τ."""
    return _k.rotate_z(vector, angle)
