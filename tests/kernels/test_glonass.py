"""Kernel ``glonass`` and the GLONASS part of ``toolkit.orbits``.

Checked against GLONASS ICD 5.1 A.3.1.2 with the two corrections of RTCM
10403.3 3.5.13 and its RK4 note, by independent routes (the acceleration
as the gradient of the PZ-90 potential written here; conservation of the
Jacobi integral, which the uncorrected equations break; convergence in
the step; the ICD's own estimate of the luni-solar effect), and against
public precise orbits at the metre level.

Public data used (values copied into this file, nothing else committed):

- Broadcast ephemerides: IGS/BKG merged navigation file
  ``BRDC00WRD_S_20240010000_01D_MN.rnx`` (RINEX 3.05, 2024-01-01), from
  https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2024/001/. RINEX gives tb as a
  UTC epoch, km, km/s, km/s^2, and -tau_n.
- Precise orbits: GFZ multi-GNSS rapid product
  ``GBM0MGXRAP_20240010000_01D_05M_ORB.SP3`` (IGS20, GPS time), from
  ftp://ftp.gfz-potsdam.de/pub/GNSS/products/mgex/2295_IGS20/, openly
  available under the IGS data policy (https://igs.org/data-access/),
  with attribution to the IGS and GFZ.

tb is brought to GPS time with GPS - UTC = 18 s (IS-GPS-200N with IRN-IS-200N-003 and -004
20.3.3.5.2.4; the value in force since 2017). A precise orbit is the centre of mass, the
broadcast one the antenna phase centre (GLONASS ICD 3.3.4), some 2 m away
on GLONASS-M; PZ-90.11 and ITRF agree within centimetres. "Metre level"
here is a bound of 5 m.
"""

import math

import pytest

from ssrcore._kernels import glonass as g
from ssrcore.toolkit import orbits
from ssrcore.toolkit import KernelError
from ssrcore.toolkit.orbits import GlonassOrbit

# GLONASS ICD 5.1 A.3.1.2, written out here.
MU = 398600.44e9
AE = 6378136.0
J2 = 1082625.7e-9
OMEGA = 7.292115e-5

STATE = ((-11175834.47266, -9508662.597656, -20866544.92188),
         (2857.489585876, -39.36195373535, -1511.75403595))
ZERO = (0.0, 0.0, 0.0)


def potential(x, y, z):
    """PZ-90 potential to J2 plus the centrifugal potential of the rotating
    frame: mu/r - mu J2 ae^2 (3 z^2/r^2 - 1) / (2 r^3) + omega^2 (x^2 + y^2) / 2."""
    r = math.sqrt(x * x + y * y + z * z)
    return (MU / r - MU * J2 * AE * AE * (3 * z * z / (r * r) - 1) / (2 * r ** 3)
            + OMEGA * OMEGA * (x * x + y * y) / 2)


def gradient(x, y, z, h=1.0):
    return ((potential(x + h, y, z) - potential(x - h, y, z)) / (2 * h),
            (potential(x, y + h, z) - potential(x, y - h, z)) / (2 * h),
            (potential(x, y, z + h) - potential(x, y, z - h)) / (2 * h))


def jacobi(p, v):
    """Jacobi integral of motion in a uniformly rotating frame, constant for
    a conservative field without the luni-solar term."""
    return (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) / 2 - potential(*p)


# ---- equations of motion ------------------------------------------------------

@pytest.mark.parametrize("p", [STATE[0], (7283031.8, 20786628.5, -12921561.1),
                               (0.0, 15e6, 20e6), (19e6, 0.0, 19e6)])
def test_acceleration_is_gradient_of_potential_plus_coriolis(p):
    """The corrected equations are grad(potential) - 2 omega x v + a_ls:
    the J2 parts must be a gradient (RTCM's (3 - 5 z^2/r^2) correction) and
    the Coriolis term -2 omega x v = (2 omega vy, -2 omega vx, 0) (RTCM's
    sign correction of dVy/dt)."""
    v = (2857.5, -39.4, -1511.8)
    a_ls = (1.0e-6, -2.0e-6, 3.0e-6)
    rate_p, rate_v = g.derivative(p, v, a_ls)
    assert rate_p == v
    grad = gradient(*p)
    expected = (grad[0] + 2 * OMEGA * v[1] + a_ls[0],
                grad[1] - 2 * OMEGA * v[0] + a_ls[1],
                grad[2] + a_ls[2])
    assert math.dist(rate_v, expected) < 2e-8


def test_jacobi_integral_is_conserved_over_the_interval():
    """Integrating +-15 min (the ICD interval) keeps the Jacobi integral.
    Coriolis does no work only with its corrected sign, and the J2 force is
    conservative only with the corrected z term. Over 15 minutes from this
    state the ICD's printed sign changes the integral by about 2e5 m^2/s^2
    and its printed z term by about 1e2; the corrected equations keep it
    to round-off (1e-5 of 1.1e7)."""
    c0 = jacobi(*STATE)
    for span in (-900.0, -450.0, 450.0, 900.0, 2700.0):
        status, p, v = g.state(STATE[0], STATE[1], ZERO, 2295, 87318.0, 2295, 87318.0 + span)
        assert status == g.OK
        assert abs(jacobi(p, v) - c0) < 1e-3, span


def test_steps_are_equal_and_at_most_30_s():
    p0, v0 = STATE
    # 30 s: one step
    assert g.state(p0, v0, ZERO, 2295, 0.0, 2295, 30.0)[1:] == g.rk4_step(p0, v0, ZERO, 30.0)
    # 31 s: two steps of 15.5 s
    p1, v1 = g.rk4_step(p0, v0, ZERO, 15.5)
    assert g.state(p0, v0, ZERO, 2295, 0.0, 2295, 31.0)[1:] == g.rk4_step(p1, v1, ZERO, 15.5)
    # backwards alike
    p1, v1 = g.rk4_step(p0, v0, ZERO, -15.5)
    assert g.state(p0, v0, ZERO, 2295, 31.0, 2295, 0.0)[1:] == g.rk4_step(p1, v1, ZERO, -15.5)
    assert orbits.GLONASS_MAX_STEP == 30.0


def test_state_at_tb_is_the_broadcast_state():
    status, p, v = g.state(STATE[0], STATE[1], (1e-6, 1e-6, 1e-6), 2295, 87318.0, 2295, 87318.0)
    assert (status, p, v) == (g.OK, STATE[0], STATE[1])


def test_30_s_steps_reach_the_rtcm_precision_of_0_1_mm():
    """RTCM 3.5.13: RK4 with steps not above 30 s gives 0.1 mm. Compare with
    1 s steps over 15 minutes."""
    a_ls = (2.8e-6, -1.9e-6, 9.3e-7)
    _, p30, v30 = g.state(STATE[0], STATE[1], a_ls, 2295, 0.0, 2295, 900.0)
    p, v = STATE
    for _ in range(900):
        p, v = g.rk4_step(p, v, a_ls, 1.0)
    assert math.dist(p30, p) < 1e-4
    assert math.dist(v30, v) < 1e-7


def test_there_and_back_returns_to_the_start():
    _, p, v = g.state(STATE[0], STATE[1], ZERO, 2295, 0.0, 2295, 900.0)
    _, p, v = g.state(p, v, ZERO, 2295, 900.0, 2295, 0.0)
    assert math.dist(p, STATE[0]) < 1e-6 and math.dist(v, STATE[1]) < 1e-9


def test_lunisolar_effect_is_the_icd_estimate():
    """GLONASS ICD 5.1 A.3.1.1 note: the luni-solar acceleration J adds
    about J tau^2 / 2 to the position and J tau to the velocity, within 10%."""
    a_ls = (2.8e-6, -1.9e-6, 9.3e-7)
    tau = 900.0
    _, p0, v0 = g.state(STATE[0], STATE[1], ZERO, 2295, 0.0, 2295, tau)
    _, p1, v1 = g.state(STATE[0], STATE[1], a_ls, 2295, 0.0, 2295, tau)
    for j in range(3):
        assert p1[j] - p0[j] == pytest.approx(a_ls[j] * tau * tau / 2, rel=0.1)
        assert v1[j] - v0[j] == pytest.approx(a_ls[j] * tau, rel=0.1)


def test_velocity_is_the_derivative_of_position():
    o = GlonassOrbit(2295, 87318.0, STATE[0], STATE[1], (2.8e-6, -1.9e-6, 9.3e-7))
    t = 87318.0 + 600.0
    h = 0.5
    s = orbits.glonass_state(o, 2295, t)
    p1 = orbits.glonass_state(o, 2295, t + h).position
    p0 = orbits.glonass_state(o, 2295, t - h).position
    assert math.dist(s.velocity, [(b - a) / (2 * h) for a, b in zip(p0, p1)]) < 1e-5


def test_span_beyond_the_bound_is_refused():
    o = GlonassOrbit(2295, 0.0, STATE[0], STATE[1], ZERO)
    with pytest.raises(KernelError) as err:
        orbits.glonass_state(o, 2295, orbits.GLONASS_MAX_SPAN + 1.0)
    assert (err.value.kernel, err.value.status) == ("glonass", g.TOO_FAR)
    assert g.state(STATE[0], STATE[1], ZERO, 2295, 0.0, 2295, float("nan"))[0] == g.TOO_FAR


# ---- clock --------------------------------------------------------------------

def test_clock_is_minus_tau_plus_gamma_dt():
    """GLONASS ICD 5.1 3.3.3: t_GLONASS = t + tau_n - gamma_n (t - tb), so
    the satellite clock offset t - t_GLONASS = -tau_n + gamma_n (t - tb)."""
    o = GlonassOrbit(2295, 604500.0, STATE[0], STATE[1], ZERO, tau_n=-7.3e-5, gamma_n=9.1e-13)
    assert orbits.glonass_clock(o, 2295, 604500.0) == 7.3e-5
    dt = 400.0  # across the week end
    assert orbits.glonass_clock(o, 2296, 100.0) == pytest.approx(7.3e-5 + 9.1e-13 * dt, abs=1e-20)
    assert orbits.glonass_clock_drift(o) == 9.1e-13


# ---- public broadcast ephemerides against public precise orbits ---------------
# (tb as GPS week, seconds; position m, velocity m/s, luni-solar m/s^2 from km;
# tau_n is minus the RINEX "-TauN"; precise: GPS week, sow, x, y, z m, clock us)

PUBLIC = {
    "R01": dict(
        tb=(2295, 86400.0 + 900.0 + 18.0),
        position=(-11175834.47266, -9508662.597656, -20866544.92188),
        velocity=(2857.489585876, -39.36195373535, -1511.75403595),
        lunisolar=(-2.793967723846e-06, -1.862645149231e-06, 0.0),
        tau_n=-7.346831262112e-05, gamma_n=9.094947017729e-13,
        sp3=[(2295, 86700.0, -12899555.791, -9528742.824, -19838102.219, 73.473071),
             (2295, 87300.0, -11227238.908, -9507992.587, -20839253.070, 73.473228),
             (2295, 88200.0, -8586799.402, -9639213.241, -22001180.827, 73.474054)]),
    "R09": dict(
        tb=(2295, 86400.0 + 8100.0 + 18.0),
        position=(7285992.1875, 20818277.83203, -12869019.53125),
        velocity=(165.5836105347, 1754.727363586, 2921.821594238),
        lunisolar=(0.0, 0.0, 9.313225746155e-07),
        tau_n=-1.619076356292e-04, gamma_n=2.728484105319e-12,
        sp3=[(2295, 93900.0, 7209046.094, 19657293.045, -14613155.698, 161.898552),
             (2295, 94500.0, 7283031.838, 20786628.534, -12921561.061, 161.899579),
             (2295, 95400.0, 7471381.243, 22200516.190, -10180327.231, 161.901486)]),
}


@pytest.mark.parametrize("name", sorted(PUBLIC))
def test_broadcast_orbit_against_precise_orbit(name):
    r = PUBLIC[name]
    o = GlonassOrbit(*r["tb"], r["position"], r["velocity"], r["lunisolar"],
                     r["tau_n"], r["gamma_n"])
    for week, sow, x, y, z, _ in r["sp3"]:
        s = orbits.glonass_state(o, week, sow)
        assert math.dist(s.position, (x, y, z)) < 5.0, (name, sow)


@pytest.mark.parametrize("name", sorted(PUBLIC))
def test_broadcast_clock_against_precise_clock(name):
    """Precise clocks are against GPS time, the broadcast one against
    GLONASS time; the two time scales differ by tens of nanoseconds at
    most (GLONASS ICD 3.3.3 gives the GPS-GLONASS difference within 30 ns),
    so the sign and size of -tau_n + gamma_n (t - tb) show."""
    r = PUBLIC[name]
    o = GlonassOrbit(*r["tb"], r["position"], r["velocity"], r["lunisolar"],
                     r["tau_n"], r["gamma_n"])
    for week, sow, _, _, _, clock_us in r["sp3"]:
        assert abs(orbits.glonass_clock(o, week, sow) - clock_us * 1e-6) < 50e-9
