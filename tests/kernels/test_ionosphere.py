"""The ionosphere kernel and ``toolkit.ionosphere``: SPARTN v2.0.3 8.11.

TEC to metres per frequency against the published factors of TM-23 Table 4.2 and
the cancellation in the ionosphere-free combination (TM-23 5.4.1.1), an
independent route. The polynomial and the grid residual are the area
kernel's (tests/kernels/test_area.py); the end-to-end test here joins them.
"""

import math

import pytest

from ssrcore._kernels import ionosphere as k
from ssrcore.toolkit import area, ionosphere
from ssrcore.toolkit.area import Area

C = (12.4, -1.6, 0.8, 0.05)  # C00 TECU, C01 TECU/deg, C10 TECU/deg, C11 TECU/deg^2


# ---- TEC to metres (8.11) --------------------------------------------------------

F0 = 10.23e6
GPS_L1, GPS_L2, GPS_L5 = 154 * F0, 120 * F0, 115 * F0
GAL_E5B = 118 * F0
GLO_G1, GLO_G2 = 1602e6, 1246e6


@pytest.mark.parametrize("f, alpha", [
    (GPS_L1, 0.1624), (GPS_L2, 0.2674), (GPS_L5, 0.2912),
    (GAL_E5B, 0.2766), (GLO_G1, 0.1570), (GLO_G2, 0.2596),
])
def test_one_tecu_matches_tm23_table_4_2(f, alpha):
    # TM-23 Table 4.2: alpha_f = 40.3 10^16 / f^2, metres of delay per TECU,
    # published to four decimals.
    code, phase = ionosphere.delay(1.0, f)
    assert round(code, 4) == alpha
    assert phase == -code


def test_code_delayed_phase_advanced_and_grid_added():
    code, phase = ionosphere.delay(20.0, GPS_L1, tec_grid=-0.5)
    assert code == pytest.approx(40.3e16 / GPS_L1 ** 2 * 19.5)
    assert code > 0 and phase == -code
    assert ionosphere.delay(20.0, GPS_L1) == ionosphere.delay(20.0, GPS_L1, 0.0)


def test_ionosphere_free_combination_cancels_the_delay():
    # TM-23 5.4.1.1: (f1^2 P1 - f2^2 P2) / (f1^2 - f2^2) removes the
    # first-order ionosphere, the 1/f^2 dependence of 8.11.
    for f1, f2 in ((GPS_L1, GPS_L2), (GPS_L1, GPS_L5), (GLO_G1, GLO_G2)):
        d1, _ = ionosphere.delay(37.5, f1)
        d2, _ = ionosphere.delay(37.5, f2)
        assert (f1 ** 2 * d1 - f2 ** 2 * d2) / (f1 ** 2 - f2 ** 2) == pytest.approx(0.0, abs=1e-12)
        # The geometry-free difference gives the TEC back.
        assert (d2 - d1) / (40.3e16 * (1 / f2 ** 2 - 1 / f1 ** 2)) == pytest.approx(37.5)


@pytest.mark.parametrize("f", [0.0, -1.0, float("nan")])
def test_frequency_must_be_positive(f):
    assert k.delay(1.0, 0.0, f)[0] == k.BAD_FREQUENCY
    with pytest.raises(ValueError):
        ionosphere.delay(1.0, f)


def test_slant_delay_end_to_end():
    # 8.11 in full on a synthetic area: I_P from the polynomial about the
    # centroid, I_G from the grid, the delay on GPS L2.
    a = Area(48.0, -3.0, 3, 3, 1.0, 1.0)            # centroid 47 N, 2 W
    lat, lon = 46.5, -1.25
    dlat, dlon = area.offsets(a, lat, lon)
    ip = area.polynomial(2, C, dlat, dlon)
    nodes, weights = area.grid(a, lat, lon)
    residuals = [0.04 * n for n in range(9)]
    ig = area.interpolate(weights, [residuals[n] for n in nodes])
    code, phase = ionosphere.delay(ip, GPS_L2, ig)
    by_hand_ip = 12.4 - 1.6 * -0.5 + 0.8 * 0.75 + 0.05 * -0.5 * 0.75
    by_hand_ig = 0.04 * (1.5 * 3 + 1.75)
    assert code == pytest.approx(40.3e16 / GPS_L2 ** 2 * (by_hand_ip + by_hand_ig))
    assert math.isclose(phase, -code)


# ---- IGS SSR VTEC (IGS SSR v1.00 4.5.1) ---------------------------------------------
#
# Each step against an independent route: the Legendre functions from Rodrigues'
# formula in exact integers, normalized by their definition; the pierce point as
# the intersection of the line of sight with the layer's sphere in Cartesian
# coordinates, which needs neither the central angle nor the branch conditions of
# the document's longitude formula.

RE = 6370e3


def _legendre(n: int, m: int, x: float) -> float:
    """P_nm(x) = (1 − x²)^(m/2) dᵐ/dxᵐ Pₙ(x), Pₙ = dⁿ/dxⁿ (x² − 1)ⁿ / (2ⁿ n!), no
    Condon–Shortley phase; the polynomial in exact integers."""
    # (x² − 1)ⁿ = Σ_k C(n, k) (−1)^(n−k) x^(2k)
    coefficients = {2 * kk: math.comb(n, kk) * (-1) ** (n - kk) for kk in range(n + 1)}
    for _ in range(n + m):
        coefficients = {p - 1: c * p for p, c in coefficients.items() if p > 0}
    poly = sum(c * x ** p for p, c in coefficients.items())
    return (1 - x * x) ** (m / 2) * poly / (2 ** n * math.factorial(n))


def _normalized(n: int, m: int, x: float) -> float:
    """The 4π normalization: √((2 − δm0)(2n + 1)(n − m)!/(n + m)!) P_nm."""
    return math.sqrt((2 - (m == 0)) * (2 * n + 1) * math.factorial(n - m) / math.factorial(n + m)) \
        * _legendre(n, m, x)


def _vtec(degree, order, cosine, sine, lat, lon, t):
    lam = (lon + (t - 50400.0) * math.pi / 43200.0) % (2 * math.pi)
    total, ic, js = 0.0, 0, 0
    for m in range(min(order, degree) + 1):
        for n in range(m, degree + 1):
            p = _normalized(n, m, math.sin(lat))
            total += cosine[ic] * math.cos(m * lam) * p
            ic += 1
            if m:
                total += sine[js] * math.sin(m * lam) * p
                js += 1
    return total


def _coefficients(degree, order):
    m = min(order, degree)
    nc = (m + 1) * (degree + 1) - m * (m + 1) // 2
    cosine = [((7 * k) % 11 - 5) * 0.37 for k in range(nc)]
    sine = [((5 * k) % 13 - 6) * 0.21 for k in range(nc - degree - 1)]
    return cosine, sine


@pytest.mark.parametrize("degree, order", [(1, 1), (4, 2), (8, 8), (15, 15), (16, 16), (16, 3), (3, 7)])
def test_vtec_is_the_sum_of_fully_normalized_harmonics(degree, order):
    cosine, sine = _coefficients(degree, order)
    for lat, lon, t in ((0.3, 0.2, 50400.0), (-1.2, 2.9, 3600.5), (1.5, -3.0, 86399.0), (0.0, 0.0, 0.0)):
        value = ionosphere.vtec(degree, order, cosine, sine, lat, lon, t)
        assert value == pytest.approx(_vtec(degree, order, cosine, sine, lat, lon, t), rel=1e-11, abs=1e-9)


def test_vtec_sun_fixed_longitude_moves_with_time():
    # IGS 4.5.1: λS = λPP + (t − 50400)·π/43200, so an hour moves the map 15° west.
    cosine, sine = _coefficients(3, 3)
    a = ionosphere.vtec(3, 3, cosine, sine, 0.4, 0.5, 10000.0)
    b = ionosphere.vtec(3, 3, cosine, sine, 0.4, 0.5 - math.radians(15.0), 13600.0)
    assert a == pytest.approx(b, rel=1e-12)


def test_vtec_refuses_a_degree_beyond_the_field_and_a_wrong_count():
    assert k.vtec(17, 1, [0.0] * 35, [0.0] * 17, 0.0, 0.0, 0.0)[0] == k.BAD_DEGREE
    with pytest.raises(ValueError):
        ionosphere.vtec(2, 2, [1.0] * 5, [1.0] * 3, 0.0, 0.0, 0.0)


def _intersection(lat, lon, height, elevation, azimuth, layer):
    """The pierce point by geometry: the point of the ray from the receiver whose
    distance from the centre is RE + layer."""
    up = (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))
    east = (-math.sin(lon), math.cos(lon), 0.0)
    north = (-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat))
    d = [math.cos(elevation) * math.sin(azimuth) * e + math.cos(elevation) * math.cos(azimuth) * nn
         + math.sin(elevation) * u for e, nn, u in zip(east, north, up)]
    r = [(RE + height) * u for u in up]
    b = sum(ri * di for ri, di in zip(r, d))
    c = sum(ri * ri for ri in r) - (RE + layer) ** 2
    s = -b + math.sqrt(b * b - c)
    p = [ri + s * di for ri, di in zip(r, d)]
    return math.atan2(p[2], math.hypot(p[0], p[1])), math.atan2(p[1], p[0])


def _wrap(angle):
    return (angle + math.pi) % (2 * math.pi) - math.pi


@pytest.mark.parametrize("lat_deg", [-89.0, -60.0, -10.0, 0.0, 10.0, 52.0, 85.0, 89.5])
@pytest.mark.parametrize("az_deg", [0.0, 10.0, 90.0, 170.0, 180.0, 200.0, 300.0, 359.0])
@pytest.mark.parametrize("el_deg", [3.0, 20.0, 60.0, 89.0])
def test_pierce_point_lies_on_the_line_of_sight(lat_deg, az_deg, el_deg):
    lat, lon, el, az = map(math.radians, (lat_deg, -100.0, el_deg, az_deg))
    look = ionosphere.SphericalLook(lat, lon, 250.0, el, az)
    p = ionosphere.pierce(look, 450e3)
    lat_g, lon_g = _intersection(lat, lon, 250.0, el, az, 450e3)
    assert p.latitude == pytest.approx(lat_g, abs=1e-9)
    assert _wrap(p.longitude - lon_g) == pytest.approx(0.0, abs=1e-8)
    # STEC = VTEC / sin(E + ψ): the incidence at the layer (IGS 4.5.1)
    assert p.slant_factor == pytest.approx(1.0 / math.sin(el + p.psi), rel=1e-12)
    assert p.slant_factor >= 1.0


def test_pierce_point_at_zenith_is_overhead_and_unmapped():
    look = ionosphere.SphericalLook(0.7, 0.3, 0.0, math.pi / 2, 0.0)
    p = ionosphere.pierce(look, 450e3)
    assert (p.latitude, p.longitude) == pytest.approx((0.7, 0.3), abs=1e-12)
    assert p.psi == pytest.approx(0.0, abs=1e-12) and p.slant_factor == pytest.approx(1.0)


def test_a_layer_below_the_receiver_has_no_pierce_point():
    look = ionosphere.SphericalLook(0.7, 0.3, 500e3, 0.1, 0.0)
    with pytest.raises(ionosphere.KernelError):
        ionosphere.pierce(look, 450e3)


def test_spherical_look_uses_the_geocentric_vertical():
    # A receiver on the equator: the geocentric and geodetic verticals agree, so the
    # ellipsoidal line of sight comes back unchanged; at 45° they differ by ~0.19°.
    from ssrcore.toolkit import geodesy
    for lat_deg, tolerance in ((0.0, 1e-12), (45.0, math.radians(0.2))):
        lat, lon = math.radians(lat_deg), math.radians(30.0)
        receiver = geodesy.geodetic_to_ecef((lat, lon, 100.0))
        el, az = math.radians(40.0), math.radians(120.0)
        enu = (math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el))
        s = ionosphere.spherical(receiver, geodesy.enu_to_ecef(enu, lat, lon))
        assert s.elevation == pytest.approx(el, abs=tolerance)
        assert s.azimuth == pytest.approx(az, abs=tolerance)
        assert s.longitude == pytest.approx(lon, abs=1e-12)
        assert s.latitude == pytest.approx(math.atan2(receiver[2], math.hypot(receiver[0], receiver[1])))
        assert s.height == pytest.approx(math.dist(receiver, (0, 0, 0)) - RE)
    assert k.spherical((0.0, 0.0, 0.0), (0.0, 0.0, 1.0))[0] == k.BAD_POSITION
