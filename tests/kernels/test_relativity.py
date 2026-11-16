"""The relativity kernel and toolkit module (DESIGN 2.7).

Every check takes an independent route where one exists:

- the clock term -2 r.v / c^2 (SPARTN 8.7) against the broadcast-element form
  F e sqrt(A) sin E (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1), with r, v and E from
  the `kepler` kernel on synthetic elements (IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV);
- F against the values the ICDs print;
- the Shapiro delay against a numerical integration of Ashby 2003 Eq. (95)
  along the path, and the magnitudes Ashby states in section 10;
- the Sagnac term against the range to the satellite position rotated by
  the Earth's rotation during the travel time, against the same range
  solved in the inertial frame with the receiver turning, and against
  Ashby's equatorial example, Eq. (8); the remainder of Ashby's
  first-order Eq. (7) is measured over the sky.
"""

import math

import pytest

from ssrcore._kernels import relativity as kernel
from ssrcore.toolkit import KernelError, geodesy, orbits, relativity, ssr

C = 299792458.0

# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1 and Table 20-IV] GPS constants.
MU_GPS = 3.986005e14
OMEGA_E_GPS = 7.2921151467e-5
# [spec: Galileo OS SIS ICD 2.2 5.1.1 Table 68] and [spec: BDS-SIS-ICD-B1I-3.0 5.2.4.9]
MU_GALILEO_BEIDOU = 3.986004418e14

# Synthetic broadcast elements of GPS magnitudes (not from any service).
ELEMENTS = dict(
    sqrt_a=5153.65, e=0.0123, m0=0.71, delta_n=4.5e-9, omega=-1.83,
    omega0=2.42, omega_dot=-8.1e-9, i0=0.958, idot=2.1e-10,
    cuc=-1.9e-6, cus=8.6e-6, crc=221.0, crs=-38.0, cic=1.1e-7, cis=-9.0e-8,
)
UNPERTURBED = dict(ELEMENTS, delta_n=0.0, cuc=0.0, cus=0.0, crc=0.0, crs=0.0, cic=0.0, cis=0.0)


# toe in the middle of a week, so that every tk of TIMES is in the same week.
WEEK, TOE = 2300, 302400.0


def state(el, tk):
    """Position, velocity and E_k at tk seconds from toe, from the ``kepler``
    kernel (IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV, the velocity of sheet 3), whose
    own tests check it against the table and a public ephemeris."""
    s = orbits.kepler_state(orbits.KeplerOrbit(orbits.GPS, WEEK, TOE, **el), WEEK, TOE + tk)
    return s.position, s.velocity, s.eccentric_anomaly


TIMES = [-7200.0 + 900.0 * k for k in range(17)]


# ---- clock term ---------------------------------------------------------------

def test_f_equals_the_icd_values():
    # [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] F = -4.442807633e-10 s/m^1/2
    # for mu = 3.986005e14.
    assert relativity.clock_kepler(MU_GPS, 1.0, 1.0, math.pi / 2) == pytest.approx(-4.442807633e-10, abs=5e-20)
    # [spec: Galileo OS SIS ICD 2.2 5.1.4 Eq. 16] F = -4.442807309e-10 for mu = 3.986004418e14.
    assert relativity.clock_kepler(MU_GALILEO_BEIDOU, 1.0, 1.0, math.pi / 2) == pytest.approx(
        -4.442807309e-10, abs=5e-20)


@pytest.mark.parametrize("tk", TIMES)
def test_vector_form_equals_icd_form_on_a_keplerian_orbit(tk):
    # Without harmonic terms the two forms are identical (IS-GPS-200N with IRN-IS-200N-003 and -004:
    # "alternative but equivalent"); r.v is the same in ECEF and ECI.
    p, v, ek = state(UNPERTURBED, tk)
    vector = relativity.clock(p, v)
    icd = relativity.clock_kepler(MU_GPS, UNPERTURBED["e"], UNPERTURBED["sqrt_a"], ek)
    assert abs(vector - icd) * C < 1e-6  # metres


@pytest.mark.parametrize("tk", TIMES)
def test_vector_form_equals_icd_form_on_a_broadcast_orbit(tk):
    # With every term of Table 20-IV the forms differ by r times the rate of
    # the radial harmonic crc cos 2phi + crs sin 2phi (centimetre level), and
    # by the factor delta_n / n0 (3e-5): the along-track and inclination
    # harmonics and the node rates move the satellite across r, not along it.
    # Bound: 2 r_max * 2 |C_r| * phi_dot_max / c, in metres, plus 0.1 mm.
    p, v, ek = state(ELEMENTS, tk)
    vector = relativity.clock(p, v)
    icd = relativity.clock_kepler(MU_GPS, ELEMENTS["e"], ELEMENTS["sqrt_a"], ek)
    a, e = ELEMENTS["sqrt_a"] ** 2, ELEMENTS["e"]
    n = math.sqrt(MU_GPS / a**3)
    phi_dot_max = n * (1 + e) ** 2 / (1 - e * e) ** 1.5
    bound = 2 * a * (1 + e) * 2 * math.hypot(ELEMENTS["crc"], ELEMENTS["crs"]) * phi_dot_max / C
    assert abs(vector) * C > 0.1  # not trivially small at these times
    assert abs(vector - icd) * C < bound * 1.01 + 1e-4
    # Without the radial harmonic, r.v = a^2 e sin E (n0 + delta_n), so the
    # vector form is the ICD form times (n0 + delta_n) / n0 exactly; the
    # along-track and inclination harmonics leave r.v unchanged.
    no_radial = dict(ELEMENTS, crc=0.0, crs=0.0)
    p, v, ek = state(no_radial, tk)
    vector = relativity.clock(p, v)
    icd = relativity.clock_kepler(MU_GPS, no_radial["e"], no_radial["sqrt_a"], ek)
    assert abs(vector - icd * (n + ELEMENTS["delta_n"]) / n) * C < 1e-6


def test_clock_sign_and_ashby_magnitude():
    # [spec: Ashby 2003 section 5, Eq. (39)] 4.4428e-10 e sqrt(a) sin E; about
    # 23 ns at most for e = 0.01. Ashby's sign is that of a correction to the
    # transmitted time; the tr of IS-GPS-200N with IRN-IS-200N-003 and -004 has the opposite sign.
    peak = relativity.clock_kepler(geodesy.WGS84_GM, 0.01, math.sqrt(26_561_750.0), math.pi / 2)
    assert -23.5e-9 < peak < -22.5e-9
    # Moving away from the Earth (r.v > 0, after perigee) the term is negative.
    assert relativity.clock((2.6e7, 0.0, 0.0), (10.0, 3870.0, 0.0)) == pytest.approx(-2 * 2.6e8 / C**2)


def test_clock_is_frame_invariant():
    # IS-GPS-200N with IRN-IS-200N-003 and -004: "immaterial whether ... earth fixed, rotating
    # coordinates or earth-centered, inertial": omega x r is orthogonal to r.
    p, v, _ = state(ELEMENTS, 1234.0)
    w = OMEGA_E_GPS
    v_inertial = (v[0] - w * p[1], v[1] + w * p[0], v[2])
    assert relativity.clock(p, v_inertial) == pytest.approx(relativity.clock(p, v), rel=1e-12, abs=1e-22)


# ---- Shapiro ------------------------------------------------------------------

RECEIVER = (4_027_894.0, 307_045.0, 4_919_474.0)  # a point near the surface


def satellite_at(elevation_deg, radius=26_560_000.0):
    """A satellite at the given elevation above the receiver's geocentric horizon."""
    r2 = math.sqrt(sum(c * c for c in RECEIVER))
    up = tuple(c / r2 for c in RECEIVER)
    east = (-up[1], up[0], 0.0)
    ne = math.hypot(*east)
    east = tuple(c / ne for c in east)
    el = math.radians(elevation_deg)
    # distance l along the direction (cos el east + sin el up): |rcv + l d| = radius
    b = r2 * math.sin(el)
    l = -b + math.sqrt(b * b - r2 * r2 + radius * radius)
    d = tuple(math.cos(el) * e + math.sin(el) * u for e, u in zip(east, up))
    return tuple(r + l * di for r, di in zip(RECEIVER, d))


def shapiro_by_integration(sat, rcv, steps=20000):
    """Ashby 2003 Eq. (95): c dt = (1 - 2V/c^2) dsigma, V = -GM/r; the delay is
    the integral of 2 GM / (c^2 r) along the straight path (Simpson's rule)."""
    def f(s):
        x = [a + s * (b - a) for a, b in zip(rcv, sat)]
        return 1.0 / math.sqrt(sum(c * c for c in x))
    total = f(0.0) + f(1.0)
    for k in range(1, steps):
        total += (4 if k % 2 else 2) * f(k / steps)
    length = math.dist(sat, rcv)
    return 2 * geodesy.WGS84_GM / C**2 * total * length / (3 * steps)


@pytest.mark.parametrize("elevation", [0.0, 5.0, 15.0, 30.0, 60.0, 90.0])
def test_shapiro_equals_integral_of_the_metric(elevation):
    sat = satellite_at(elevation)
    assert relativity.shapiro(sat, RECEIVER) == pytest.approx(shapiro_by_integration(sat, RECEIVER), abs=1e-7)


def test_shapiro_magnitudes_of_ashby_section_10():
    # GM/c^2 is the 4.43 mm of Ashby Eq. (97); the coefficient of Eq. (96) is twice that.
    assert geodesy.WGS84_GM / C**2 == pytest.approx(4.43e-3, abs=1e-5)
    zenith = relativity.shapiro(satellite_at(90.0), RECEIVER)
    horizon = relativity.shapiro(satellite_at(0.0), RECEIVER)
    assert 0.010 < zenith < horizon < 0.020
    # "The net effect for a satellite to earth link is less than 2 cm": the
    # whole of Eq. (96), with Phi_0 / c^2 of Eq. (18).
    phi0 = -6.96927e-10
    for elevation in (0.0, 10.0, 30.0, 90.0):
        sat = satellite_at(elevation)
        net = phi0 * math.dist(sat, RECEIVER) + relativity.shapiro(sat, RECEIVER)
        assert abs(net) < 0.02


def test_shapiro_through_the_centre_is_refused():
    with pytest.raises(ValueError):
        relativity.shapiro((2.6e7, 0.0, 0.0), (-6.4e6, 0.0, 0.0))
    status, delay = kernel.shapiro(2.6e7, 0.0, 0.0, -6.4e6, 0.0, 0.0)
    assert status == kernel.THROUGH_CENTRE and math.isnan(delay)


# ---- Sagnac -------------------------------------------------------------------

def rotated_range(sat, rcv, w=geodesy.WGS84_OMEGA_E):
    """Exact light-time range: the satellite position of the transmit epoch's
    ECEF frame rotated into the frame of the reception epoch, iterated."""
    tau = math.dist(sat, rcv) / C
    for _ in range(10):
        a = w * tau
        s = (math.cos(a) * sat[0] + math.sin(a) * sat[1],
             -math.sin(a) * sat[0] + math.cos(a) * sat[1], sat[2])
        tau = math.dist(s, rcv) / C
    return C * tau


def inertial_range(sat, rcv, w=geodesy.WGS84_OMEGA_E):
    """The same range by a second route: in the inertial frame that coincides
    with ECEF at transmission the satellite stays put and the receiver turns
    by +w tau (Ashby 2003 section 2); c tau = |r_s - Rz(w tau) r_r| solved by
    bisection on tau."""
    def gap(tau):
        a = w * tau
        r = (math.cos(a) * rcv[0] - math.sin(a) * rcv[1], math.sin(a) * rcv[0] + math.cos(a) * rcv[1], rcv[2])
        return C * tau - math.dist(sat, r)
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if gap(mid) < 0.0 else (lo, mid)
    return C * 0.5 * (lo + hi)


def eq7(sat, rcv):
    """Ashby 2003 Eq. (7), first order in w tau."""
    return geodesy.WGS84_OMEGA_E * (sat[0] * rcv[1] - sat[1] * rcv[0]) / C


def turned(sat, azimuth_turn):
    """The satellite moved about the receiver's vertical by a fraction of a turn."""
    up = tuple(c / math.sqrt(sum(x * x for x in RECEIVER)) for c in RECEIVER)
    a = 2 * math.pi * azimuth_turn
    rel = tuple(s - r for s, r in zip(sat, RECEIVER))
    dot = sum(p * q for p, q in zip(rel, up))
    cross = (up[1] * rel[2] - up[2] * rel[1], up[2] * rel[0] - up[0] * rel[2], up[0] * rel[1] - up[1] * rel[0])
    rel = tuple(p * math.cos(a) + c * math.sin(a) + u * dot * (1 - math.cos(a)) for p, c, u in zip(rel, cross, up))
    return tuple(r + p for r, p in zip(RECEIVER, rel))


@pytest.mark.parametrize("radius", [26_560e3, 42_164e3], ids=["GPS", "GEO"])
@pytest.mark.parametrize("elevation", [5.0, 20.0, 45.0, 90.0])
@pytest.mark.parametrize("azimuth_turn", [0.0, 0.25, 0.5, 0.75])
def test_sagnac_equals_rotation_during_travel(radius, elevation, azimuth_turn):
    # The satellite east, north, west and south of the receiver.
    sat = turned(satellite_at(elevation, radius), azimuth_turn)
    exact = rotated_range(sat, RECEIVER) - math.dist(sat, RECEIVER)
    sagnac = relativity.sagnac(sat, RECEIVER)
    assert abs(sagnac) < 60.0
    # The two routes in double precision differ by the rounding of ranges of
    # 4e7 m (some 1e-8 m); the kernel's form avoids that subtraction.
    assert sagnac == pytest.approx(exact, abs=1e-7)
    assert sagnac == pytest.approx(inertial_range(sat, RECEIVER) - math.dist(sat, RECEIVER), abs=1e-7)


def sky(radius, inclination):
    """Satellites of one orbit radius over receivers from the equator to 80
    degrees, everywhere above 5 degrees elevation."""
    for lat_r in range(0, 81, 10):
        rcv = (6_378_137.0 * math.cos(math.radians(lat_r)), 0.0, 6_378_137.0 * math.sin(math.radians(lat_r)))
        for lat_s in range(-inclination, inclination + 1, 5):
            for lon_s in range(-90, 91, 3):
                p, q = math.radians(lat_s), math.radians(lon_s)
                sat = (radius * math.cos(p) * math.cos(q), radius * math.cos(p) * math.sin(q), radius * math.sin(p))
                d = [s - r for s, r in zip(sat, rcv)]
                if sum(a * b / 6_378_137.0 for a, b in zip(d, rcv)) > math.sin(math.radians(5)) * math.dist(sat, rcv):
                    yield sat, rcv


@pytest.mark.parametrize("radius, inclination, low, high", [
    (26_560e3, 55, 0.08e-3, 0.12e-3),   # GPS
    (29_600e3, 56, 0.11e-3, 0.15e-3),   # Galileo
    (42_164e3, 0, 0.25e-3, 0.32e-3),    # GEO
    (42_164e3, 55, 0.25e-3, 0.32e-3),   # IGSO
], ids=["GPS", "Galileo", "GEO", "IGSO"])
def test_remainder_of_the_first_order_form(radius, inclination, low, high):
    # DESIGN 2.4 row E: the exact form, because the remainder of Eq. (7)
    # reaches the 0.1 mm at which the engine is compared with a reference.
    # Measured here over the sky: its largest value lies in [low, high].
    largest = max(abs(relativity.sagnac(s, r) - eq7(s, r)) for s, r in sky(radius, inclination))
    assert low < largest < high


def test_sagnac_around_the_equator_is_ashby_eq_8():
    # [spec: Ashby 2003 section 2, Eq. (8)] 2 omega_E / c^2 = 1.6227e-21 s/m^2;
    # once round the equator eastward, 207.4 ns.
    assert 2 * geodesy.WGS84_OMEGA_E / C**2 == pytest.approx(1.6227e-21, abs=5e-26)
    a1, n = 6_378_137.0, 3600
    points = [(a1 * math.cos(2 * math.pi * k / n), a1 * math.sin(2 * math.pi * k / n), 0.0) for k in range(n + 1)]
    east = sum(relativity.sagnac(points[k], points[k + 1]) for k in range(n)) / C
    west = sum(relativity.sagnac(points[k + 1], points[k]) for k in range(n)) / C
    assert east == pytest.approx(207.4e-9, abs=0.05e-9)
    assert west == pytest.approx(-207.4e-9, abs=0.05e-9)


def test_sagnac_of_coincident_points_is_refused():
    with pytest.raises(KernelError) as raised:
        relativity.sagnac(RECEIVER, RECEIVER)
    assert raised.value.status == kernel.ZERO_RANGE
    status, delta = kernel.sagnac(*RECEIVER, *RECEIVER)
    assert status == kernel.ZERO_RANGE and math.isnan(delta)


# ---- the constant part of Eq. (96), not applied (DESIGN 2.4, row S) ------------

@pytest.mark.parametrize("radius, low, high", [
    (26_560e3, 0.0140, 0.0181),   # GPS: ranges 20 200 to 25 800 km
    (42_164e3, 0.0249, 0.0291),   # GEO: ranges 35 800 to 41 700 km
], ids=["GPS", "GEO"])
def test_phi0_term_magnitude(radius, low, high):
    # Phi_0 l / c^2 with Phi_0 / c^2 = -6.96927e-10 (Ashby 2003 Eq. 18): the
    # range-proportional part of Eq. (96) the record does not carry; DESIGN
    # 2.4 states its size so that a reference applying it is recognised.
    phi0 = -6.96927e-10
    for elevation in (5.0, 90.0):
        size = -phi0 * math.dist(satellite_at(elevation, radius), RECEIVER)
        assert low < size < high


def test_constants():
    # The values Ashby 2003 evaluates with; one definition each (DESIGN 4).
    assert ssr.SPEED_OF_LIGHT == C
    assert geodesy.WGS84_GM == 3.986004418e14  # Ashby 2003 Eq. (13)
    assert geodesy.WGS84_OMEGA_E == 7.2921151467e-5  # Ashby 2003 section 2
