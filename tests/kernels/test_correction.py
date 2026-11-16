"""The ``correction`` kernel (toolkit ``ssr``) against SPARTN v2.0.3 sections 8.7 and 8.8, on synthetic states.

Independent routes: on a circular orbit the radial, along and cross axes are
known in closed form from the orbital elements (position direction, velocity
direction, orbit normal Rz(RAAN) Rx(i) z); on an eccentric orbit the frame
must be orthonormal and right-handed, and e_radial must equal the part of r
orthogonal to the velocity, normalised (a second formulation of
e_along x e_cross).
"""

import math

import pytest

from ssrcore._kernels import correction as kernel
from ssrcore.toolkit import KernelError, ssr

# m^3/s^2, IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3; any value serves a synthetic orbit
GM = 3.986005e14


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(k, a):
    return tuple(k * x for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def rotate(raan, inc, v):
    """Rz(raan) Rx(inc) v."""
    x, y, z = v
    y, z = y * math.cos(inc) - z * math.sin(inc), y * math.sin(inc) + z * math.cos(inc)
    x, y = x * math.cos(raan) - y * math.sin(raan), x * math.sin(raan) + y * math.cos(raan)
    return (x, y, z)


def kepler_state(a, e, inc, raan, argp, nu):
    """Position and velocity of a two-body orbit at true anomaly nu (perifocal, rotated)."""
    p = a * (1 - e * e)
    r = p / (1 + e * math.cos(nu))
    pos = (r * math.cos(nu), r * math.sin(nu), 0.0)
    k = math.sqrt(GM / p)
    vel = (-k * math.sin(nu), k * (e + math.cos(nu)), 0.0)
    rot = lambda v: rotate(raan, inc, rotate(argp, 0.0, v))  # noqa: E731
    return rot(pos), rot(vel)


CIRCULAR = [  # (radius m, inclination, RAAN, argument of latitude), GPS-, Galileo-, BeiDou-like
    (26_559_800.0, math.radians(55.0), math.radians(30.0), math.radians(10.0)),
    (29_600_000.0, math.radians(56.0), math.radians(200.0), math.radians(250.0)),
    (42_164_000.0, math.radians(0.5), math.radians(120.0), math.radians(95.0)),
    (42_164_000.0, math.radians(55.0), math.radians(300.0), math.radians(180.0)),
]

ECCENTRIC = [  # (a, e, i, RAAN, argument of perigee, true anomaly)
    (26_559_800.0, 0.02, math.radians(55.0), 1.0, 2.0, 0.3),
    (42_164_000.0, 0.08, math.radians(57.0), 4.0, 0.5, 2.9),
    (27_906_100.0, 0.16, math.radians(41.0), 5.5, 4.7, 4.0),  # QZSS-like eccentricity
]


def assert_vec(a, b, tol):
    assert all(abs(x - y) <= tol for x, y in zip(a, b)), (a, b)


@pytest.mark.parametrize("radius, inc, raan, u", CIRCULAR)
def test_circular_orbit_frame_is_position_velocity_normal(radius, inc, raan, u):
    speed = math.sqrt(GM / radius)
    r = rotate(raan, inc, (radius * math.cos(u), radius * math.sin(u), 0.0))
    v = rotate(raan, inc, (-speed * math.sin(u), speed * math.cos(u), 0.0))
    er, ea, ec = ssr.orbit_frame(r, v)
    assert_vec(er, (math.cos(raan) * math.cos(u) - math.sin(raan) * math.cos(inc) * math.sin(u),
                    math.sin(raan) * math.cos(u) + math.cos(raan) * math.cos(inc) * math.sin(u),
                    math.sin(inc) * math.sin(u)), 1e-14)
    assert_vec(ea, scale(1 / speed, v), 1e-14)
    normal = (math.sin(inc) * math.sin(raan), -math.sin(inc) * math.cos(raan), math.cos(inc))
    assert_vec(ec, normal, 1e-14)


@pytest.mark.parametrize("a, e, inc, raan, argp, nu", ECCENTRIC)
def test_frame_orthonormal_right_handed_and_signs(a, e, inc, raan, argp, nu):
    r, v = kepler_state(a, e, inc, raan, argp, nu)
    er, ea, ec = ssr.orbit_frame(r, v)
    for x in (er, ea, ec):
        assert norm(x) == pytest.approx(1.0, abs=1e-15)
    for x, y in ((er, ea), (er, ec), (ea, ec)):
        assert abs(dot(x, y)) < 1e-15
    assert dot(cross(er, ea), ec) == pytest.approx(1.0, abs=1e-15)  # radial, along, cross right-handed
    assert dot(er, r) > 0.0  # radial points away from the Earth
    assert dot(ea, v) == pytest.approx(norm(v), rel=1e-15)  # along is the direction of motion
    assert dot(ec, cross(r, v)) > 0.0  # cross is the angular momentum direction
    # Second formulation of e_radial: the part of r orthogonal to v, normalised.
    vh = scale(1 / norm(v), v)
    rp = add(r, scale(-dot(r, vh), vh))
    assert_vec(er, scale(1 / norm(rp), rp), 1e-15)
    # Eccentric orbits: radial differs from r's direction by the flight-path angle.
    gamma = math.atan2(e * math.sin(nu), 1 + e * math.cos(nu))
    assert dot(er, scale(1 / norm(r), r)) == pytest.approx(math.cos(gamma), abs=1e-14)


@pytest.mark.parametrize("a, e, inc, raan, argp, nu", ECCENTRIC)
def test_correction_is_the_frame_times_do(a, e, inc, raan, argp, nu):
    r, v = kepler_state(a, e, inc, raan, argp, nu)
    er, ea, ec = ssr.orbit_frame(r, v)
    d_o = (0.731, -1.284, 0.052)
    d_x = ssr.orbit_correction(r, v, d_o)
    assert_vec(d_x, add(add(scale(d_o[0], er), scale(d_o[1], ea)), scale(d_o[2], ec)), 1e-15)
    # The frame is orthonormal, so projecting back recovers dO and the length is kept.
    assert_vec((dot(d_x, er), dot(d_x, ea), dot(d_x, ec)), d_o, 1e-15)
    assert norm(d_x) == pytest.approx(norm(d_o), rel=1e-15)
    # X = X_b - dX.
    x = ssr.corrected_position(r, v, d_o)
    assert_vec(x, add(r, scale(-1.0, d_x)), 1e-8)


def test_sign_of_each_component_on_a_circular_orbit():
    """[spec: SPARTN 8.8, X = X_b - dX] A positive radial correction lowers the
    satellite, a positive along correction moves it back along its track, a
    positive cross correction moves it against the orbit normal."""
    radius, inc, raan, u = CIRCULAR[0]
    speed = math.sqrt(GM / radius)
    r = rotate(raan, inc, (radius * math.cos(u), radius * math.sin(u), 0.0))
    v = rotate(raan, inc, (-speed * math.sin(u), speed * math.cos(u), 0.0))
    h = cross(r, v)

    x = ssr.corrected_position(r, v, (1.0, 0.0, 0.0))
    assert norm(x) == pytest.approx(radius - 1.0, abs=1e-8)
    x = ssr.corrected_position(r, v, (0.0, 1.0, 0.0))
    assert dot(add(x, scale(-1, r)), v) / speed == pytest.approx(-1.0, abs=1e-8)
    x = ssr.corrected_position(r, v, (0.0, 0.0, 1.0))
    assert dot(add(x, scale(-1, r)), h) / norm(h) == pytest.approx(-1.0, abs=1e-8)
    assert ssr.corrected_position(r, v, (0.0, 0.0, 0.0)) == r


def test_ecef_velocity_with_earth_rotation():
    """The frame uses the velocity it is given (8.8: the broadcast velocity);
    with an ECEF velocity, which includes -omega x r, it is still orthonormal."""
    omega = 7.2921151467e-5
    r, v = kepler_state(*ECCENTRIC[0])
    v_ecef = add(v, scale(-1.0, cross((0.0, 0.0, omega), r)))
    er, ea, ec = ssr.orbit_frame(r, v_ecef)
    assert_vec(ea, scale(1 / norm(v_ecef), v_ecef), 1e-15)
    assert abs(dot(er, ea)) < 1e-15 and abs(dot(er, ec)) < 1e-15 and abs(dot(ea, ec)) < 1e-15


def test_degenerate_states():
    r = (26_559_800.0, 0.0, 0.0)
    calls = [
        lambda: ssr.orbit_frame(r, (0.0, 0.0, 0.0)),
        lambda: ssr.orbit_frame(r, (3.0, 0.0, 0.0)),  # r parallel to r_dot
        lambda: ssr.orbit_correction(r, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        lambda: ssr.corrected_position(r, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
    ]
    for call in calls:
        with pytest.raises(KernelError) as raised:
            call()
        assert raised.value.status == kernel.DEGENERATE
    assert kernel.orbit_frame(r, (0.0, 0.0, 0.0))[0] == kernel.DEGENERATE


def test_clock_correction():
    """[spec: SPARTN 8.7, t_s = t_b - dC/c; Table 8.0 c] A positive dC of one
    light-nanosecond makes the satellite clock one nanosecond earlier."""
    assert ssr.SPEED_OF_LIGHT == 299792458.0
    assert ssr.corrected_clock(1.0e-4, 0.299792458) == pytest.approx(1.0e-4 - 1.0e-9, abs=1e-20)
    assert ssr.corrected_clock(-2.5e-5, -16.382) == pytest.approx(-2.5e-5 + 16.382 / 299792458.0, abs=1e-20)
    assert ssr.corrected_clock(3.0e-4, 0.0) == 3.0e-4
    # The SF020 range, ±16.382 m, is ±54.6 ns of clock.
    assert (ssr.corrected_clock(0.0, -16.382) - ssr.corrected_clock(0.0, 16.382)) == pytest.approx(109.29e-9, rel=1e-4)


# ---- rates and the clock polynomial (rule C11) ---------------------------------------------
#
# RTCM 10403.3 3.5.13.5 eq. (3.12-6) and 3.5.13.6 eq. (3.12-8), the same in IGS
# SSR v1.00 4.1 and 4.2; values at the ends of the fields' ranges (DF365 to
# DF370, DF376 to DF378; IDF013 to IDF021) and in between.


@pytest.mark.parametrize("d_o, rate, dt", [
    ((0.1234, -0.5, 0.25), (0.001, -0.004, 0.002), 7.5),
    ((209.7151, -209.7148, 209.7148), (1.048575, -1.048572, 1.048572), 2.5),
    ((-1.0, 2.0, -3.0), (0.5, 0.25, -0.125), -2.5),       # t before t0: t0 is mid-interval
])
def test_orbit_at_is_rtcm_eq_3_12_6(d_o, rate, dt):
    out = ssr.orbit_at(d_o, rate, dt)
    assert out == tuple(o + r * dt for o, r in zip(d_o, rate))
    assert kernel.orbit_at(d_o, rate, dt) == out


@pytest.mark.parametrize("c0, c1, c2, dt", [
    (1.234, 0.001, 0.00002, 5.0),
    (-209.7151, 1.048575, -1.34217726, 30.0),
    (0.5, -0.25, 0.125, -15.0),
])
def test_clock_at_is_rtcm_eq_3_12_8_and_its_rate_the_derivative(c0, c1, c2, dt):
    assert ssr.clock_at(c0, c1, c2, dt) == c0 + c1 * dt + c2 * dt * dt
    assert ssr.clock_rate_at(c1, c2, dt) == c1 + 2.0 * c2 * dt
    h = 1e-3
    numeric = (ssr.clock_at(c0, c1, c2, dt + h) - ssr.clock_at(c0, c1, c2, dt - h)) / (2 * h)
    assert numeric == pytest.approx(ssr.clock_rate_at(c1, c2, dt), abs=1e-9)


def test_at_the_reference_time_the_values_are_those_sent():
    d_o = (0.5, -0.25, 0.124)
    assert ssr.orbit_at(d_o, (0.9, -0.8, 0.7), 0.0) == d_o
    assert ssr.clock_at(1.234, 0.9, -0.8, 0.0) == 1.234
    assert ssr.clock_rate_at(0.9, -0.8, 0.0) == 0.9
    with pytest.raises(ValueError):
        ssr.orbit_at((0.0, 0.0), (0.0, 0.0, 0.0), 1.0)
