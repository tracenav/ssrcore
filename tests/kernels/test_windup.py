"""Wind-up kernel and toolkit (DESIGN 2.7).

Checks, each by a route other than the kernel's own arithmetic:

- the satellite frame of SPARTN v2.0.3 section 8.14 against its figure
  (e_x the direction of motion, e_y perpendicular to the orbit plane, e_z to
  the Earth's centre), with the ECI velocity obtained independently by
  differencing positions in a rotating frame;
- the nominal Sun frame of TM-23/1 5.6.3 against a second formulation;
- the wind-up of TM-23/1 5.5 against the literal sign/arccos form of
  eq. (5.72), and against its defining property: turning the satellite or
  the receiver antenna about the line of sight by an angle changes the phase
  by that angle;
- continuity across +-1/2 cycle with eq. (5.73).

All inputs are synthetic.
"""

import math
import random

import pytest

from ssrcore._kernels import windup as k
from ssrcore.toolkit import windup as w

OMEGA = 7.2921151467e-5   # IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV
R_GPS = 26_560_000.0      # m, a GPS-like orbit radius
MU = 3.986005e14          # IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV


# ---- vector helpers, independent of the kernel --------------------------------

def add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def scale(s, u):
    return tuple(s * a for a in u)


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def norm(u):
    return math.sqrt(dot(u, u))


def unit(u):
    return scale(1.0 / norm(u), u)


def rotate(v, axis, angle):
    """Rodrigues: v turned by ``angle`` right-handed about the unit ``axis``."""
    c, s = math.cos(angle), math.sin(angle)
    return add(add(scale(c, v), scale(s, cross(axis, v))), scale(dot(axis, v) * (1 - c), axis))


def rz(theta, v):
    return rotate(v, (0.0, 0.0, 1.0), theta)


def close(u, v, tol):
    return all(abs(a - b) <= tol for a, b in zip(u, v))


def enu(lat, lon):
    """East and north unit vectors (TM-23/1 3.1.4), written here for the tests."""
    e = (-math.sin(lon), math.cos(lon), 0.0)
    n = (-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat))
    return e, n


def wrap(c):
    """A value in cycles into (-1/2, 1/2]."""
    c = c - math.floor(c + 0.5)
    return 0.5 if c == -0.5 else c


def random_state(rng):
    """A satellite ECEF state on a circular inertial orbit of random orientation."""
    r = unit((rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1)))
    t = unit(cross(r, unit((rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1)))))
    pos = scale(R_GPS, r)
    v_eci = scale(math.sqrt(MU / R_GPS), t)
    v_ecef = sub(v_eci, cross((0.0, 0.0, OMEGA), pos))
    return pos, v_ecef


# ---- SPARTN 8.14 satellite frame ----------------------------------------------

def test_spartn_frame_axes_of_the_figure():
    # Equatorial prograde orbit at (R, 0, 0), moving towards +y in ECI.
    pos = (R_GPS, 0.0, 0.0)
    v_eci = (0.0, math.sqrt(MU / R_GPS), 0.0)
    v_ecef = sub(v_eci, cross((0.0, 0.0, OMEGA), pos))
    x, y, z = w.satellite_frame(pos, v_ecef)
    assert close(z, (-1.0, 0.0, 0.0), 1e-15)   # e_z: direction to Earth centre
    assert close(x, (0.0, 1.0, 0.0), 1e-15)    # e_x: direction of satellite motion
    orbit_normal = unit(cross(pos, v_eci))     # (0, 0, 1)
    assert close(y, scale(-1.0, orbit_normal), 1e-15)   # e_y: perpendicular to the orbit plane


def test_spartn_frame_uses_the_inertial_velocity():
    # Independent route: the ECI velocity as a central difference of ECEF
    # positions carried into an inertial frame that coincides with ECEF at t = 0.
    pos = (12_000_000.0, -15_000_000.0, 18_000_000.0)
    vel = (1500.0, 2100.0, -800.0)
    h = 1e-3

    def eci(t):
        return rz(OMEGA * t, add(pos, scale(t, vel)))

    v_eci = scale(1.0 / (2 * h), sub(eci(h), eci(-h)))
    ez = scale(-1.0, unit(pos))
    ey = unit(cross(ez, v_eci))
    ex = cross(ey, ez)
    x, y, z = w.satellite_frame(pos, vel)
    assert close(y, ey, 1e-9) and close(x, ex, 1e-9) and close(z, ez, 1e-15)
    # The Earth-fixed velocity alone would turn this frame by 0.86 degree.
    x0 = cross(unit(cross(ez, vel)), ez)
    assert math.degrees(math.acos(dot(x0, x))) > 0.5


def test_spartn_frame_orthonormal_and_right_handed():
    rng = random.Random(8_14)
    for _ in range(200):
        pos, vel = random_state(rng)
        yaw = rng.uniform(-2 * math.pi, 2 * math.pi)
        x, y, z = w.satellite_frame(pos, vel, yaw)
        for u in (x, y, z):
            assert abs(norm(u) - 1.0) < 1e-15
        assert abs(dot(x, y)) < 1e-15 and abs(dot(y, z)) < 1e-15 and abs(dot(x, z)) < 1e-15
        assert close(cross(x, y), z, 1e-15)


def test_spartn_yaw_turns_the_frame_about_e_z():
    rng = random.Random(21)
    pos, vel = random_state(rng)
    x0, y0, z0 = w.satellite_frame(pos, vel, 0.0)
    for yaw in (math.radians(6.0), math.radians(90.0), math.radians(354.0), 2.0):
        x, y, z = w.satellite_frame(pos, vel, yaw)
        assert close(z, z0, 0.0)
        assert close(x, rotate(x0, z0, yaw), 1e-15)
        assert close(y, rotate(y0, z0, yaw), 1e-15)
    x, y, _ = w.satellite_frame(pos, vel, math.pi / 2)
    assert close(x, y0, 1e-15) and close(y, scale(-1.0, x0), 1e-15)


def test_spartn_frame_degenerate():
    with pytest.raises(ValueError):
        w.satellite_frame((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    # Inertial velocity along the radius: no orbit plane.
    pos = (0.0, 0.0, R_GPS)
    with pytest.raises(ValueError):
        w.satellite_frame(pos, (0.0, 0.0, 3000.0))
    status, *_ = k.frame_spartn((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), 0.0)
    assert status == k.DEGENERATE


# ---- TM-23 5.6.3 Sun frame ----------------------------------------------------

def test_sun_frame_against_projection_form():
    # Second formulation: i is the Sun direction with its k component removed.
    rng = random.Random(5_76)
    for _ in range(200):
        pos, _ = random_state(rng)
        sun = scale(1.496e11, unit((rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1))))
        i, j, kk = w.satellite_frame_sun(pos, sun)
        e = unit(sub(sun, pos))
        assert close(kk, scale(-1.0, unit(pos)), 1e-15)
        assert close(i, unit(sub(e, scale(dot(kk, e), kk))), 1e-12)
        assert abs(dot(j, e)) < 1e-15            # the Sun lies in the i-k plane
        assert dot(i, e) > 0.0                   # on the +i side
        assert close(cross(i, j), kk, 1e-15)


def test_sun_frame_degenerate():
    pos = (R_GPS, 0.0, 0.0)
    with pytest.raises(ValueError):
        w.satellite_frame_sun(pos, (1.496e11, 0.0, 0.0))   # Sun along the radius
    with pytest.raises(ValueError):
        w.satellite_frame_sun(pos, pos)


# ---- TM-23 5.5 wind-up --------------------------------------------------------

def literal_windup(r_sat, r_rcv, a, b, a_s, b_s):
    """Eqs. (5.72) and (5.74) of TM-23/1 as printed, in cycles."""
    rho = unit(sub(r_rcv, r_sat))
    d = add(sub(a, scale(dot(rho, a), rho)), cross(rho, b))
    ds = sub(sub(a_s, scale(dot(rho, a_s), rho)), cross(rho, b_s))
    zeta = dot(rho, cross(ds, d))
    c = max(-1.0, min(1.0, dot(ds, d) / (norm(ds) * norm(d))))
    return math.copysign(1.0, zeta) * math.acos(c) / (2 * math.pi)


def random_geometry(rng):
    lat, lon = rng.uniform(-1.4, 1.4), rng.uniform(-math.pi, math.pi)
    up = (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))
    rcv = scale(6_371_000.0, up)
    while True:
        pos, vel = random_state(rng)
        if dot(unit(sub(pos, rcv)), up) > 0.1:   # above the horizon
            break
    east, north = enu(lat, lon)
    x, y, _ = w.satellite_frame(pos, vel, rng.uniform(0, 2 * math.pi))
    return pos, rcv, east, north, x, y


def test_windup_equals_the_printed_form():
    rng = random.Random(5_72)
    n = 0
    while n < 300:
        g = random_geometry(rng)
        expected = literal_windup(*g)
        if abs(abs(expected) - 0.5) < 1e-6:
            continue
        assert abs(w.windup(*g) - expected) < 1e-12
        n += 1


def test_turning_an_antenna_about_the_line_of_sight():
    # Turning the satellite dipoles right-handed about rho (satellite to
    # receiver) by alpha lowers the phase by alpha; turning the receiver's
    # raises it by alpha. This is the property the formula expresses.
    rng = random.Random(1993)
    for _ in range(100):
        pos, rcv, a, b, x, y = random_geometry(rng)
        rho = unit(sub(rcv, pos))
        w0 = w.windup(pos, rcv, a, b, x, y)
        alpha = rng.uniform(-math.pi, math.pi)
        ws = w.windup(pos, rcv, a, b, rotate(x, rho, alpha), rotate(y, rho, alpha))
        wr = w.windup(pos, rcv, rotate(a, rho, alpha), rotate(b, rho, alpha), x, y)
        assert abs(wrap(ws - (w0 - alpha / (2 * math.pi)))) < 1e-12
        assert abs(wrap(wr - (w0 + alpha / (2 * math.pi)))) < 1e-12


def test_zenith_aligned_and_opposed():
    # Receiver on the equator at longitude 0: east = +y, north = +z, up = +x.
    rcv = (6_378_137.0, 0.0, 0.0)
    sat = (R_GPS, 0.0, 0.0)
    east, north = (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    # Satellite frame with k = -up: i = east and j = -north is right-handed.
    assert w.windup(sat, rcv, east, north, east, scale(-1.0, north)) == 0.0
    # Turned by half a turn: exactly +1/2, the closed end of (-1/2, 1/2].
    assert w.windup(sat, rcv, east, north, scale(-1.0, east), north) == 0.5
    # A quarter turn of the satellite about its boresight (= rho here).
    rho = (-1.0, 0.0, 0.0)
    q = math.pi / 2
    v = w.windup(sat, rcv, east, north, rotate(east, rho, q), rotate(scale(-1.0, north), rho, q))
    assert abs(v + 0.25) < 1e-15


def test_yaw_at_zenith_is_a_phase_rotation():
    # At zenith e_z is the line of sight, so SF021 yaw psi changes the
    # wind-up by -psi/(2 pi) exactly.
    rcv = (6_378_137.0, 0.0, 0.0)
    pos = (R_GPS, 0.0, 0.0)
    v_eci = (0.0, math.sqrt(MU / R_GPS) * 0.6, math.sqrt(MU / R_GPS) * 0.8)
    vel = sub(v_eci, cross((0.0, 0.0, OMEGA), pos))
    east, north = (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    x, y, _ = w.satellite_frame(pos, vel)
    w0 = w.windup(pos, rcv, east, north, x, y)
    for deg in range(0, 360, 6):
        psi = math.radians(deg)
        x, y, _ = w.satellite_frame(pos, vel, psi)
        assert abs(wrap(w.windup(pos, rcv, east, north, x, y) - w0 + deg / 360.0)) < 1e-12


def test_windup_degenerate():
    p = (R_GPS, 0.0, 0.0)
    e, n = (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    with pytest.raises(ValueError):
        w.windup(p, p, e, n, e, n)                       # no line of sight
    # Receiver dipoles with b = rho x a: d = a + rho x (rho x a) = 0.
    rcv = (6_378_137.0, 0.0, 0.0)
    rho = (-1.0, 0.0, 0.0)
    status, _ = k.cycles(p, rcv, e, cross(rho, e), e, n)
    assert status == k.DEGENERATE
    with pytest.raises(ValueError):
        w.windup(p, rcv, e, cross(rho, e), e, n)


# ---- TM-23 5.5 continuity -----------------------------------------------------

def test_unwrap_worked_values():
    assert w.unwrap(0.49, -0.49) == pytest.approx(0.51, abs=1e-15)
    assert w.unwrap(-0.49, 0.49) == pytest.approx(-0.51, abs=1e-15)
    assert w.unwrap(10.2, 0.25) == 10.25
    assert w.unwrap(-3.9, 0.3) == pytest.approx(-3.7, abs=1e-15)
    assert w.unwrap(0.0, 0.3) == 0.3        # N = 0 at the start of an arc
    assert w.unwrap(0.5, 0.0) == 1.0        # a tie rounds away from zero
    assert w.unwrap(-0.5, 0.0) == -1.0


def test_continuity_across_half_cycle():
    # The satellite antenna turns about the line of sight by 0.07 cycle a step
    # for 3 turns: the raw values wrap at +-1/2 six times; unwrapped they are
    # the straight line the rotation gives.
    rcv = (6_378_137.0, 0.0, 0.0)
    sat = (R_GPS, 0.0, 0.0)
    east, north = (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    rho = (-1.0, 0.0, 0.0)
    prev = prev_raw = None
    wraps = 0
    for n in range(43):
        a = 2 * math.pi * 0.07 * n
        raw = w.windup(sat, rcv, east, north, rotate(east, rho, a), rotate(scale(-1.0, north), rho, a))
        assert -0.5 < raw <= 0.5
        cur = raw if prev is None else w.unwrap(prev, raw)
        if prev_raw is not None and abs(raw - prev_raw) > 0.5:
            wraps += 1
        assert abs(cur + 0.07 * n) < 1e-12
        prev, prev_raw = cur, raw
    assert wraps == 3   # 2.94 cycles downwards cross -1/2 at 0.5, 1.5 and 2.5


def test_continuity_along_a_pass():
    # A satellite on a circular orbit passing near a receiver, sampled every
    # 30 s with the SPARTN frame: the unwrapped wind-up never steps by more
    # than a few hundredths of a cycle.
    lat, lon = math.radians(47.0), math.radians(8.0)
    up = (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))
    rcv = scale(6_371_000.0, up)
    east, north = enu(lat, lon)
    # At t = 0 the satellite is 3 degrees east of the zenith direction,
    # moving north in inertial space.
    s = unit(add(up, scale(0.05, east)))
    n_hat = unit(cross(s, north))
    n_mean = math.sqrt(MU / R_GPS ** 3)
    prev = None
    seen = 0
    for i in range(-480, 481):
        t = 30.0 * i
        # Inertial circular orbit through s at t = 0, carried to ECEF.
        u = rotate(s, n_hat, n_mean * t)
        pos_i = scale(R_GPS, u)
        vel_i = scale(R_GPS * n_mean, cross(n_hat, u))
        pos = rz(-OMEGA * t, pos_i)
        vel = sub(rz(-OMEGA * t, vel_i), cross((0.0, 0.0, OMEGA), pos))
        if dot(unit(sub(pos, rcv)), up) < 0.05:
            prev = None
            continue
        x, y, _ = w.satellite_frame(pos, vel)
        cur = w.windup(pos, rcv, east, north, x, y)
        if prev is not None:
            cur = w.unwrap(prev, cur)
            assert abs(cur - prev) < 0.06
        prev = cur
        seen += 1
    assert seen > 100
