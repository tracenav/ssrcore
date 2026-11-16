"""The geodesy kernel and toolkit: TM-23 section 3.1.4 and Appendix B.

TM-23 Vol. I gives the formulae of 3.1.4 in Appendix B but its numerical
examples are in Vol. II (programs car2geo.f, geo2car.f), which is not
staged. The results are checked instead by independent routes:

- published WGS-84 values derived from a and f (semi-minor axis, first
  eccentricity squared; NIMA TR8350.2, Table 3.1), reached through the
  conversions themselves (the pole of the ellipsoid lies at b);
- a second, closed-form conversion from Cartesian to ellipsoidal coordinates
  (Heikkinen 1982, as given by Zhu 1994), written here independently of the
  iteration of TM-23 B.1.2;
- the ellipsoid normal from the gradient of the ellipsoid's equation, which
  must be the up vector of B.9 at the foot point;
- the rotation matrices of TM-23 eq. 3.6 composed as in B.7 and B.10,
  against the closed matrices B.8 and B.11 the kernel uses;
- round trips to 0.1 mm.
"""

import math

import pytest

from ssrcore._kernels import geodesy as k
from ssrcore.toolkit import geodesy as g
from ssrcore.toolkit import KernelError

A = 6378137.0
F = 1.0 / 298.257223563
B = A * (1.0 - F)
E2 = 2.0 * F - F * F

# The ARP of the RTCM 10403.3 section 4.2 example of message 1005, m.
RTCM_1005_ARP = (1114104.5999, -4850729.7108, 3975521.4643)

TENTH_MM = 1e-4


def heikkinen(x, y, z, a=A, f=F):
    """Closed-form ECEF to (lat, lon, h): Heikkinen (1982), in the form of
    J. Zhu, IEEE Trans. AES 30 (1994), independent of TM-23 B.1.2."""
    b = a * (1.0 - f)
    e2 = 2.0 * f - f * f
    ep2 = (a * a - b * b) / (b * b)
    p = math.hypot(x, y)
    ff = 54.0 * b * b * z * z
    gg = p * p + (1.0 - e2) * z * z - e2 * (a * a - b * b)
    c = e2 * e2 * ff * p * p / gg**3
    s = (1.0 + c + math.sqrt(c * c + 2.0 * c)) ** (1.0 / 3.0)
    pp = ff / (3.0 * (s + 1.0 / s + 1.0) ** 2 * gg * gg)
    q = math.sqrt(1.0 + 2.0 * e2 * e2 * pp)
    r0 = -pp * e2 * p / (1.0 + q) + math.sqrt(
        0.5 * a * a * (1.0 + 1.0 / q) - pp * (1.0 - e2) * z * z / (q * (1.0 + q)) - 0.5 * pp * p * p
    )
    u = math.sqrt((p - e2 * r0) ** 2 + z * z)
    v = math.sqrt((p - e2 * r0) ** 2 + (1.0 - e2) * z * z)
    z0 = b * b * z / (a * v)
    h = u * (1.0 - b * b / (a * v))
    lat = math.atan2(z + ep2 * z0, p)
    return lat, math.atan2(y, x), h


def r1(t):
    """TM-23 eq. 3.6: rotation about the first axis."""
    c, s = math.cos(t), math.sin(t)
    return ((1.0, 0.0, 0.0), (0.0, c, s), (0.0, -s, c))


def r3(t):
    """TM-23 eq. 3.6: rotation about the third axis."""
    c, s = math.cos(t), math.sin(t)
    return ((c, s, 0.0), (-s, c, 0.0), (0.0, 0.0, 1.0))


def matmul(m, n):
    return tuple(tuple(sum(m[i][j] * n[j][c] for j in range(3)) for c in range(3)) for i in range(3))


def apply(m, v):
    return tuple(sum(m[i][j] * v[j] for j in range(3)) for i in range(3))


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def close(u, v, tol):
    return all(abs(a - b) <= tol for a, b in zip(u, v))


LATITUDES = [-90.0, -89.9999, -60.0, -45.0, -10.0, 0.0, 1e-9, 23.5, 45.0, 52.0, 89.9999, 90.0]
LONGITUDES = [-180.0, -135.0, -90.0, -0.5, 0.0, 30.0, 90.0, 179.9999]
HEIGHTS = [-6000.0, 0.0, 0.0001, 8848.0, 400e3, 20200e3, 35786e3]


# ---- constants -----------------------------------------------------------------


def test_wgs84_constants_of_the_documents():
    # IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; TM-23 Table 3.2
    assert g.WGS84_A == 6378137.0
    assert g.WGS84_F == 1.0 / 298.257223563  # TM-23 Table 3.2


def test_eccentricity_squared_published_value():
    # NIMA TR8350.2 Table 3.1: e^2 = 0.00669437999014.
    assert abs(g.eccentricity_squared() - 0.00669437999014) < 5e-15
    # TM-23 B.3: e^2 = (a^2 - b^2) / a^2.
    assert abs(g.eccentricity_squared() - (A * A - B * B) / (A * A)) < 1e-17


# ---- ellipsoidal to Cartesian (B.1.1) -----------------------------------------


def test_axes_of_the_ellipsoid():
    # The equator at longitude 0 and 90 deg lies at a; the poles at b,
    # published as 6356752.3142 m (NIMA TR8350.2 Table 3.1).
    assert close(g.geodetic_to_ecef((0.0, 0.0, 0.0)), (A, 0.0, 0.0), 1e-9)
    x, y, z = g.geodetic_to_ecef((0.0, math.pi / 2, 0.0))
    assert abs(x) < 1e-9 and abs(y - A) < 1e-9 and z == 0.0
    x, y, z = g.geodetic_to_ecef((math.pi / 2, 0.0, 0.0))
    assert abs(x) < 1e-9 and abs(z - 6356752.3142) < 1e-4
    _, _, z = g.geodetic_to_ecef((-math.pi / 2, 0.0, 100.0))
    assert abs(z + B + 100.0) < 1e-8


@pytest.mark.parametrize("lat_deg", [-80.0, -33.0, 0.0, 12.0, 45.0, 67.0])
def test_surface_points_lie_on_the_ellipsoid(lat_deg):
    x, y, z = g.geodetic_to_ecef((math.radians(lat_deg), math.radians(17.0), 0.0))
    assert abs((x * x + y * y) / (A * A) + z * z / (B * B) - 1.0) < 1e-15


@pytest.mark.parametrize("lat_deg", [-75.0, -20.0, 0.0, 35.0, 60.0])
def test_up_vector_is_the_ellipsoid_normal(lat_deg):
    # The gradient of x^2/a^2 + y^2/a^2 + z^2/b^2 at the foot point is the
    # outward normal; B.9's u must be its direction, and the point at height
    # h must lie h along it from the foot.
    lat, lon = math.radians(lat_deg), math.radians(-71.0)
    x, y, z = g.geodetic_to_ecef((lat, lon, 0.0))
    grad = (x / (A * A), y / (A * A), z / (B * B))
    norm = math.sqrt(dot(grad, grad))
    normal = tuple(c / norm for c in grad)
    _, _, up = g.enu_basis(lat, lon)
    assert close(up, normal, 1e-15)
    above = g.geodetic_to_ecef((lat, lon, 1234.5))
    assert close(above, tuple(p + 1234.5 * n for p, n in zip((x, y, z), normal)), 1e-8)


# ---- Cartesian to ellipsoidal (B.1.2) ------------------------------------------


def test_rtcm_1005_example_against_closed_form():
    lat, lon, h = g.ecef_to_geodetic(RTCM_1005_ARP)
    ref = heikkinen(*RTCM_1005_ARP)
    assert abs(lat - ref[0]) < 1e-12
    assert abs(lon - ref[1]) < 1e-15
    assert abs(h - ref[2]) < TENTH_MM
    # The example station is in North America: about 38.8 N, 77.1 W, near sea level.
    assert 38.0 < math.degrees(lat) < 39.5 and -78.0 < math.degrees(lon) < -76.5 and -100 < h < 200
    assert close(g.geodetic_to_ecef((lat, lon, h)), RTCM_1005_ARP, TENTH_MM)


@pytest.mark.parametrize("lat_deg", [-89.0, -60.0, -45.0, -10.0, 0.5, 23.5, 45.0, 52.0, 88.0])
@pytest.mark.parametrize("h", HEIGHTS)
def test_against_closed_form(lat_deg, h):
    x, y, z = g.geodetic_to_ecef((math.radians(lat_deg), math.radians(-120.0), h))
    lat, lon, hh = g.ecef_to_geodetic((x, y, z))
    ref = heikkinen(x, y, z)
    assert abs(lat - ref[0]) < 1e-11
    assert abs(lon - ref[1]) < 1e-14
    assert abs(hh - ref[2]) < TENTH_MM


@pytest.mark.parametrize("lat_deg", LATITUDES)
@pytest.mark.parametrize("lon_deg", LONGITUDES)
@pytest.mark.parametrize("h", HEIGHTS)
def test_round_trips(lat_deg, lon_deg, h):
    geo = (math.radians(lat_deg), math.radians(lon_deg), h)
    ecef = g.geodetic_to_ecef(geo)
    lat, lon, hh = g.ecef_to_geodetic(ecef)
    # Geodetic -> ECEF -> geodetic: the same point to 0.1 mm.
    assert abs(hh - h) < TENTH_MM
    assert abs(lat - geo[0]) * (A + h) < TENTH_MM
    if abs(lat_deg) != 90.0:
        dlon = math.remainder(lon - geo[1], 2.0 * math.pi)
        assert abs(dlon) * (A + h) * math.cos(geo[0]) < TENTH_MM
    # ECEF -> geodetic -> ECEF.
    assert close(g.geodetic_to_ecef((lat, lon, hh)), ecef, TENTH_MM)


def test_poles_and_axis():
    lat, lon, h = g.ecef_to_geodetic((0.0, 0.0, B + 10.0))
    assert lat == math.pi / 2 and lon == 0.0 and abs(h - 10.0) < 1e-9
    lat, lon, h = g.ecef_to_geodetic((0.0, 0.0, -B))
    assert lat == -math.pi / 2 and abs(h) < 1e-9


def test_longitude_quadrants():
    for lon_deg in (-179.0, -91.0, -1.0, 1.0, 91.0, 179.0, 180.0):
        x, y, z = g.geodetic_to_ecef((0.3, math.radians(lon_deg), 50.0))
        assert abs(math.degrees(g.ecef_to_geodetic((x, y, z))[1]) - lon_deg) < 1e-9


def test_other_ellipsoid():
    # CGCS2000 (TM-23 Table 3.4): a = 6378137, f = 1/298.257222101.
    f = 1.0 / 298.257222101
    ecef = g.geodetic_to_ecef((0.7, 2.0, 100.0), A, f)
    lat, lon, h = g.ecef_to_geodetic(ecef, A, f)
    assert abs(lat - 0.7) < 1e-13 and abs(lon - 2.0) < 1e-14 and abs(h - 100.0) < 1e-6
    ref = heikkinen(*ecef, a=A, f=f)
    assert abs(lat - ref[0]) < 1e-12


def test_no_convergence_is_reported():
    for bad in [(0.0, 0.0, 0.0), (math.nan, 0.0, 0.0), (math.inf, 1.0, 1.0)]:
        status, _ = k.ecef_to_geodetic(*bad, A, F)
        assert status == k.NO_SOLUTION
        with pytest.raises(KernelError) as raised:
            g.ecef_to_geodetic(bad)
        assert raised.value.status == k.NO_SOLUTION


# ---- local frame (B.2) --------------------------------------------------------


def test_basis_at_origin_of_latitude_and_longitude():
    # B.9 at lat = lon = 0: e = y, n = z, u = x.
    e, n, u = g.enu_basis(0.0, 0.0)
    assert e == (0.0, 1.0, 0.0) and n == (0.0, 0.0, 1.0) and u == (1.0, 0.0, 0.0)


@pytest.mark.parametrize("lat_deg,lon_deg", [(0.0, 0.0), (52.0, 4.0), (-33.9, 151.2), (89.0, -170.0), (-90.0, 0.0)])
def test_matrices_are_the_composed_rotations(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    to_ecef = matmul(r3(-(math.pi / 2 + lon)), r1(-(math.pi / 2 - lat)))  # B.7
    to_enu = matmul(r1(math.pi / 2 - lat), r3(math.pi / 2 + lon))  # B.10
    v = (1234.5, -987.25, 4321.0)
    assert close(g.enu_to_ecef(v, lat, lon), apply(to_ecef, v), 1e-9)
    assert close(g.ecef_to_enu(v, lat, lon), apply(to_enu, v), 1e-9)


@pytest.mark.parametrize("lat_deg,lon_deg", [(12.0, 34.0), (-60.0, -150.0), (90.0, 0.0)])
def test_basis_is_orthonormal_right_handed_and_inverse(lat_deg, lon_deg):
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    e, n, u = g.enu_basis(lat, lon)
    for a_ in (e, n, u):
        assert abs(dot(a_, a_) - 1.0) < 1e-15
    assert abs(dot(e, n)) < 1e-15 and abs(dot(n, u)) < 1e-15 and abs(dot(u, e)) < 1e-15
    assert close(cross(e, n), u, 1e-15)
    v = (3.0, -4.0, 12.0)
    back = g.enu_to_ecef(g.ecef_to_enu(v, lat, lon), lat, lon)
    assert close(back, v, 1e-14)
    assert abs(dot(g.ecef_to_enu(v, lat, lon), g.ecef_to_enu(v, lat, lon)) - 169.0) < 1e-12


# ---- elevation and azimuth (B.3) ----------------------------------------------


def place(receiver_geo, elevation, azimuth, distance):
    """A satellite at the given elevation and azimuth: B.14 read backwards."""
    lat, lon, _ = receiver_geo
    enu = (
        distance * math.cos(elevation) * math.sin(azimuth),
        distance * math.cos(elevation) * math.cos(azimuth),
        distance * math.sin(elevation),
    )
    d = g.enu_to_ecef(enu, lat, lon)
    rcv = g.geodetic_to_ecef(receiver_geo)
    return rcv, tuple(r + c for r, c in zip(rcv, d))


def test_cardinal_directions():
    geo = (math.radians(48.0), math.radians(11.0), 500.0)
    for az_deg, expected in ((0.0, 0.0), (90.0, 90.0), (180.0, 180.0), (270.0, 270.0)):
        rcv, sat = place(geo, 0.0, math.radians(az_deg), 1000e3)
        el, az = g.elevation_azimuth(sat, rcv)
        assert abs(el) < 1e-12 and abs(math.degrees(az) - expected) < 1e-9
    rcv, sat = place(geo, math.pi / 2, 0.0, 20000e3)
    el, _ = g.elevation_azimuth(sat, rcv)
    assert abs(el - math.pi / 2) < 1e-12


@pytest.mark.parametrize("el_deg", [-10.0, 0.0, 5.0, 30.0, 60.0, 89.9])
@pytest.mark.parametrize("az_deg", [0.0, 1e-7, 45.0, 135.0, 225.0, 315.0, 359.9999])
def test_recovers_placed_angles(el_deg, az_deg):
    geo = (math.radians(-34.0), math.radians(-58.0), 25.0)
    rcv, sat = place(geo, math.radians(el_deg), math.radians(az_deg), 22000e3)
    unit, el, az = g.line_of_sight(sat, rcv)
    assert abs(el - math.radians(el_deg)) < 1e-12
    assert 0.0 <= az < 2.0 * math.pi
    assert abs(math.remainder(az - math.radians(az_deg), 2.0 * math.pi)) < 1e-11
    d = tuple(s - r for s, r in zip(sat, rcv))
    norm = math.sqrt(dot(d, d))
    assert close(unit, tuple(c / norm for c in d), 1e-15)


def test_elevation_is_angle_above_the_tangent_plane():
    # Independent route: elevation is the complement of the angle between
    # the line of sight and the ellipsoid normal (gradient of the ellipsoid).
    rcv_geo = (math.radians(40.0), math.radians(-105.0), 1600.0)
    rcv = g.geodetic_to_ecef(rcv_geo)
    sat = (-12_000e3, -14_000e3, 18_000e3)
    foot = g.geodetic_to_ecef((rcv_geo[0], rcv_geo[1], 0.0))
    grad = (foot[0] / (A * A), foot[1] / (A * A), foot[2] / (B * B))
    d = tuple(s - r for s, r in zip(sat, rcv))
    cos_zenith = dot(grad, d) / math.sqrt(dot(grad, grad) * dot(d, d))
    el, _ = g.elevation_azimuth(sat, rcv)
    assert abs(el - (math.pi / 2 - math.acos(cos_zenith))) < 1e-12


def test_given_geodetic_equals_computed():
    rcv = RTCM_1005_ARP
    sat = (15_000e3, -20_000e3, 10_000e3)
    computed = g.line_of_sight(sat, rcv)
    given = g.line_of_sight(sat, rcv, g.ecef_to_geodetic(rcv))
    assert computed == given


def test_zero_range_is_reported():
    status, _, _, _ = k.line_of_sight(RTCM_1005_ARP, RTCM_1005_ARP, 0.1, 0.2)
    assert status == k.ZERO_RANGE
    with pytest.raises(KernelError) as raised:
        g.elevation_azimuth(RTCM_1005_ARP, RTCM_1005_ARP)
    assert raised.value.status == k.ZERO_RANGE


# ---- range, its rate, and the frame of reception (DESIGN 2.4) ----------------------------


def test_range_is_the_norm_and_its_rate_the_projected_velocity():
    """TM-23 B.3, eq. B.13: ρ = |r_s − r_r|; the rate is u·v for a receiver at rest."""
    rho, rate = g.range_rate((3.0e6, 4.0e6, 12.0e6), (1000.0, -2000.0, 500.0), (0.0, 0.0, 0.0))
    assert rho == 13.0e6
    assert rate == pytest.approx((3 * 1000.0 - 4 * 2000.0 + 12 * 500.0) / 13.0, rel=1e-15)
    # the rate is the derivative of the range: a central difference over 1 ms agrees
    r, v, rcv = (15.6e6, -9.1e6, 18.3e6), (1500.0, 2700.0, -300.0), (4.0e6, 0.5e6, 4.9e6)
    step = 1e-3
    ahead = g.range_rate(tuple(a + b * step for a, b in zip(r, v)), v, rcv)[0]
    behind = g.range_rate(tuple(a - b * step for a, b in zip(r, v)), v, rcv)[0]
    assert g.range_rate(r, v, rcv)[1] == pytest.approx((ahead - behind) / (2 * step), abs=1e-6)
    with pytest.raises(KernelError):
        g.range_rate(rcv, v, rcv)


def test_rotate_z_turns_the_frame_not_the_point():
    """Ashby 2003 section 2: in a frame turned by +a about z a fixed point's longitude
    decreases by a; the norm and z are kept."""
    import math
    x, y, z = g.rotate_z((1.0, 0.0, 2.0), math.pi / 2)
    assert (x, y, z) == pytest.approx((0.0, -1.0, 2.0), abs=1e-16)
    v = (15.6e6, -9.1e6, 18.3e6)
    w = g.rotate_z(v, 7.2921151467e-5 * 0.07)
    assert math.hypot(*w) == pytest.approx(math.hypot(*v), rel=1e-15) and w[2] == v[2]
