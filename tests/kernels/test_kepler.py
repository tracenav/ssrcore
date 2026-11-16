"""Kernel ``kepler`` and the Keplerian part of ``toolkit.orbits``.

Checked against the ICDs (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1 and 20.3.3.4.3
Table 20-IV; Galileo OS SIS ICD 2.2 5.1.1 and 5.1.4; BDS-SIS-ICD-B1I-3.0 5.2.4.9
and 5.2.4.12), by independent routes where there is one (Kepler's
equation, vis-viva, a geostationary orbit, numerical derivatives, the
alternative relativistic form of 20.3.3.3.3.1), and against public precise
orbits at the metre level.

Public data used (values copied into this file, nothing else committed):

- Broadcast ephemerides: IGS/BKG merged navigation file
  ``BRDC00WRD_S_20240010000_01D_MN.rnx`` (RINEX 3.05, 2024-01-01), from
  https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2024/001/. RINEX gives the
  angles already in radians.
- Precise orbits: GFZ multi-GNSS rapid product
  ``GBM0MGXRAP_20240010000_01D_05M_ORB.SP3`` (IGS20, GPS time), from
  ftp://ftp.gfz-potsdam.de/pub/GNSS/products/mgex/2295_IGS20/.
  IGS products are openly available under the IGS data policy
  (https://igs.org/data-access/), with attribution to the IGS and the
  analysis centre (GFZ).

A precise orbit is the centre of mass; a broadcast orbit is the antenna
phase centre (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.3.1, Galileo OS SIS ICD 2.2
5.1.1, BDS ICD 5.2.4.12), offset by up to about 2.5 m, and broadcast orbits are good to
1 to 2 m (BeiDou GEO a few metres). So "metre level" here is a bound of
3 m for GPS, QZSS, Galileo and BeiDou MEO/IGSO and 5 m for BeiDou GEO.
"""

import math

import pytest

from ssrcore._kernels import kepler as k
from ssrcore.toolkit import orbits
from ssrcore.toolkit import KernelError
from ssrcore.toolkit.orbits import KeplerClock, KeplerOrbit

# Constants as the ICDs print them, written out here so the tests do not
# take them from the code under test.
MU = {orbits.GPS: 3.986005e14, orbits.QZSS: 3.986005e14, orbits.GALILEO: 3.986004418e14,
      orbits.BEIDOU: 3.986004418e14}
OMEGA_E = {orbits.GPS: 7.2921151467e-5, orbits.QZSS: 7.2921151467e-5, orbits.GALILEO: 7.2921151467e-5,
           orbits.BEIDOU: 7.2921150e-5}
PI_ICD = 3.1415926535898
C = 299792458.0
WEEK = 604800.0


def kepler_orbit(system, r, week, geo=False):
    """A KeplerOrbit from a RINEX 3 record dict."""
    return KeplerOrbit(
        system=system, toe_week=week, toe=r["toe"], sqrt_a=r["sqrt_a"], e=r["e"],
        m0=r["m0"], delta_n=r["dn"], omega=r["omega"], omega0=r["omega0"],
        omega_dot=r["omega_dot"], i0=r["i0"], idot=r["idot"], cuc=r["cuc"], cus=r["cus"],
        crc=r["crc"], crs=r["crs"], cic=r["cic"], cis=r["cis"], geo=geo)


# ---- public data: 2024-01-01, GPS week 2295 -----------------------------------
# Each record: the RINEX fields (s, m, rad), and three precise positions as
# (GPS week, GPS seconds of week, x, y, z in m, clock in microseconds).
RECORDS = {
    "G05": dict(
        toc=93600.0, af0=-1.562642864883e-04, af1=-1.477928890381e-12, af2=0.0,
        toe=93600.0, sqrt_a=5153.783403397, e=0.005790898343548, m0=-2.494229256397,
        dn=4.453756945425e-09, omega=1.235180510432, omega0=-0.6299845609285,
        omega_dot=-8.155696860368e-09, i0=0.9693533956084, idot=-3.928735076157e-11,
        cuc=-1.747161149979e-06, cus=3.31737101078e-06, crc=319.34375, crs=-34.125,
        cic=1.098960638046e-07, cis=4.097819328308e-08,
        sp3=[(2295, 93600.0, -10165892.577, -12970215.216, -20987562.731, -156.263473),
             (2295, 95400.0, -8514441.729, -17133498.617, -18560936.773, -156.265716),
             (2295, 97200.0, -7459642.749, -20779562.463, -14875907.434, -156.268086)]),
    "J02": dict(
        toc=93600.0, af0=-8.819624781609e-07, af1=-2.273736754432e-13, af2=0.0,
        toe=93600.0, sqrt_a=6493.216625214, e=0.07529062137473, m0=1.846571600548,
        dn=1.298268363803e-09, omega=-1.578717840941, omega0=2.742180015536,
        omega_dot=-8.357490980188e-10, i0=0.7130722967005, idot=-1.160762636137e-09,
        cuc=-4.878267645836e-06, cus=1.765787601471e-06, crc=108.90625, crs=-149.8125,
        cic=1.126900315285e-06, cis=3.281980752945e-06,
        sp3=[(2295, 93600.0, -33857438.328, 24434445.127, 11240184.620, -0.885167),
             (2295, 95400.0, -32650060.004, 25028628.443, 14472252.508, -0.885499),
             (2295, 97200.0, -31378027.151, 25336893.322, 17479393.464, -0.885852)]),
    "E03": dict(
        toc=90000.0, af0=-8.395168697461e-05, af1=-2.515321284591e-12, af2=0.0,
        toe=90000.0, sqrt_a=5440.616941452, e=0.0002315054880455, m0=-1.93780374359,
        dn=2.966552140233e-09, omega=1.011677624228, omega0=0.4691320002023,
        omega_dot=-5.541659403786e-09, i0=0.9604110396719, idot=-9.821837690392e-11,
        cuc=-6.573274731636e-06, cus=9.343028068542e-06, crc=141.59375, crs=-140.4375,
        cic=6.519258022308e-08, cis=-5.587935447693e-09,
        sp3=[(2295, 90000.0, 20015430.240, -9976729.429, -19395290.465, -83.94886),
             (2295, 91800.0, 23174365.929, -9644619.461, -15690011.550, -83.953456),
             (2295, 93600.0, 25608234.689, -9735867.506, -11206664.661, -83.958124)]),
    # BeiDou: toe and toc in BDT, week 939 of BDT (GPS week 2295).
    "C20": dict(
        toc=93600.0, af0=-2.318993210793e-05, af1=-3.502442780245e-11, af2=0.0,
        toe=93600.0, sqrt_a=5282.634473801, e=0.0005596516421065, m0=-1.819140909332,
        dn=3.604078695773e-09, omega=-0.5213434498987, omega0=-0.1089002567271,
        omega_dot=-6.85278544602e-09, i0=0.9765811735201, idot=-1.542921411727e-10,
        cuc=-4.664994776249e-06, cus=5.309004336596e-06, crc=262.5625, crs=-95.25,
        cic=7.543712854385e-08, cis=-1.629814505577e-08,
        sp3=[(2295, 93600.0, -22280635.475, 2885791.893, -16559399.331, -23.272058),
             (2295, 95400.0, -19496431.701, 340852.500, -19963087.048, -23.335675),
             (2295, 97200.0, -16659653.910, -2960547.514, -22186234.761, -23.399255)]),
    "C06": dict(  # IGSO
        toc=93600.0, af0=2.286629751325e-05, af1=2.05409023124e-11, af2=0.0,
        toe=93600.0, sqrt_a=6493.350938797, e=0.003801974235103, m0=-2.111332808389,
        dn=9.207526387575e-10, omega=-2.987112057847, omega0=1.255370911397,
        omega_dot=-1.66149777948e-09, i0=0.9460036153163, idot=6.932431620746e-10,
        cuc=1.187436282635e-06, cus=2.15214677155e-05, crc=-423.90625, crs=41.375,
        cic=-3.417953848839e-07, cis=1.355074346066e-07,
        sp3=[(2295, 93600.0, -2706008.002, 27856841.445, 31646234.563, 22.778363),
             (2295, 95400.0, -3972011.717, 25951758.289, 33075236.180, 22.815501),
             (2295, 97200.0, -5689048.736, 24442572.009, 33937785.848, 22.852595)]),
    "C02": dict(  # GEO (BDS-2)
        toc=93600.0, af0=3.015995025635e-04, af1=3.297273565295e-11, af2=0.0,
        toe=93600.0, sqrt_a=6493.454561234, e=0.0008897357620299, m0=-2.692316991158,
        dn=-1.870435053985e-09, omega=0.5946301070622, omega0=-2.175657386554,
        omega_dot=2.83690388272e-09, i0=0.06884128507013, idot=-4.071598169835e-10,
        cuc=7.162801921368e-06, cus=1.789582893252e-05, crc=-550.21875, crs=224.421875,
        cic=1.047737896442e-07, cis=-1.541338860989e-07,
        sp3=[(2295, 93600.0, 4414137.888, 41951868.889, 1153880.113, 301.532826),
             (2295, 95400.0, 4423317.968, 41952426.705, 1003205.224, 301.591952),
             (2295, 97200.0, 4432723.262, 41952183.945, 835309.697, 301.651154)]),
    "C59": dict(  # GEO (BDS-3)
        toc=93600.0, af0=-7.369089871645e-07, af1=1.767475055203e-13, af2=0.0,
        toe=93600.0, sqrt_a=6493.454147339, e=0.0001097865169868, m0=2.559401345431,
        dn=-1.592923494515e-09, omega=2.898151811953, omega0=-2.470151406647,
        omega_dot=2.769758228691e-09, i0=0.05742849059572, idot=-3.907305612105e-10,
        cuc=-1.462176442146e-06, cus=1.59814953804e-05, crc=-487.453125, crs=-37.40625,
        cic=-7.962808012962e-08, cis=6.658956408501e-08,
        sp3=[(2295, 93600.0, -32345709.692, 27048687.513, 582129.268, -0.901325),
             (2295, 95400.0, -32345357.964, 27052337.299, 422746.623, -0.901109),
             (2295, 97200.0, -32344236.344, 27055931.745, 256084.721, -0.900894)]),
}

# System, week of toe in the system's count, offset (s) of the system's time
# from GPS time, bound (m). BDT = GPS - 14 s and BDT week = GPS week - 1356:
# BDT starts at 2006-01-01 00:00:00 UTC (BDS-SIS-ICD-B1I-3.0 3.3), the start
# of GPS week 1356, when GPS - UTC was 14 s. GST and GPS weeks coincide in RINEX 3 (GST - GPS is a few
# ns, the GGTO, negligible for the orbit).
SATELLITES = {
    "G05": (orbits.GPS, 2295, 0, 0.0, 3.0),
    "J02": (orbits.QZSS, 2295, 0, 0.0, 3.0),
    "E03": (orbits.GALILEO, 2295, 0, 0.0, 3.0),
    "C20": (orbits.BEIDOU, 939, 1356, -14.0, 3.0),
    "C06": (orbits.BEIDOU, 939, 1356, -14.0, 3.0),
    "C02": (orbits.BEIDOU, 939, 1356, -14.0, 5.0),
    "C59": (orbits.BEIDOU, 939, 1356, -14.0, 5.0),
}


def public(name):
    system, week, week_offset, offset, bound = SATELLITES[name]
    geo = name[0] == "C" and orbits.beidou_geo(int(name[1:]))
    return kepler_orbit(system, RECORDS[name], week, geo), week_offset, offset, bound


# ---- constants and helpers ----------------------------------------------------

def test_semicircles_use_the_icd_pi():
    assert orbits.semicircles(1.0) == PI_ICD
    assert orbits.semicircles(-0.5) == -0.5 * PI_ICD
    assert PI_ICD != math.pi


def test_beidou_geo_satellites_of_table_4_1():
    geo = {1, 2, 3, 4, 5, 59, 60, 61, 62, 63}
    for prn in range(1, 64):
        assert orbits.beidou_geo(prn) == (prn in geo), prn
    assert not orbits.beidou_geo(0) and not orbits.beidou_geo(64)


def test_qzss_uses_the_gps_constants():
    gps = kepler_orbit(orbits.GPS, dict(RECORDS["J02"]), 2295)
    fields = {f: getattr(gps, f) for f in gps.__slots__}
    fields["system"] = orbits.QZSS
    assert orbits.QZSS != orbits.GPS
    assert orbits.kepler_state(KeplerOrbit(**fields), 2295, 93600.0) == orbits.kepler_state(gps, 2295, 93600.0)


# ---- Kepler's equation --------------------------------------------------------

@pytest.mark.parametrize("e", [0.0, 0.0002, 0.01, 0.075, 0.3, 0.5, 0.9])
def test_eccentric_anomaly_solves_keplers_equation(e):
    for j in range(-40, 41):
        m = j * 0.157
        ek = orbits.eccentric_anomaly(m, e)
        assert abs(ek - e * math.sin(ek) - m) <= 4e-15 * max(1.0, abs(m))


def test_eccentric_anomaly_of_a_circle_is_the_mean_anomaly():
    assert orbits.eccentric_anomaly(1.234, 0.0) == 1.234


def test_eccentric_anomaly_status_when_it_cannot_converge():
    status, _ = k.eccentric_anomaly(1.0, float("nan"))
    assert status == k.NOT_CONVERGED
    with pytest.raises(KernelError) as err:
        orbits.eccentric_anomaly(1.0, float("nan"))
    assert (err.value.kernel, err.value.status) == ("kepler", k.NOT_CONVERGED)


@pytest.mark.parametrize("change", [dict(e=1.0), dict(e=-0.1), dict(sqrt_a=0.0), dict(system=7)])
def test_bad_elements_are_refused(change):
    r = dict(RECORDS["G05"])
    base = kepler_orbit(orbits.GPS, r, 2295)
    fields = {f: getattr(base, f) for f in base.__slots__}
    fields.update(change)
    with pytest.raises(KernelError) as err:
        orbits.kepler_state(KeplerOrbit(**fields), 2295, 93600.0)
    assert (err.value.kernel, err.value.status) == ("kepler", k.BAD_ELEMENTS)


# ---- independent checks on synthetic orbits -----------------------------------

def pure_kepler(system, sqrt_a, e, i0, omega0, omega, m0, geo=False, toe=93600.0, week=2295):
    """Elements with no perturbation: omega_dot = 0 fixes the node in
    inertial space, so this is a fixed ellipse in inertial space."""
    return KeplerOrbit(system=system, toe_week=week, toe=toe, sqrt_a=sqrt_a, e=e, m0=m0,
                       delta_n=0.0, omega=omega, omega0=omega0, omega_dot=0.0,
                       i0=i0, idot=0.0, cuc=0.0, cus=0.0, crc=0.0, crs=0.0, cic=0.0,
                       cis=0.0, geo=geo)


@pytest.mark.parametrize("system,sqrt_a,e", [
    (orbits.GPS, 5153.7, 0.0058), (orbits.GPS, 6493.2, 0.075),
    (orbits.GALILEO, 5440.6, 0.0002), (orbits.BEIDOU, 5282.6, 0.15)])
def test_pure_kepler_orbit_obeys_vis_viva_and_keeps_angular_momentum(system, sqrt_a, e):
    """With omega_dot = 0 the node is fixed in inertial space and the
    orbit is a Kepler ellipse: r = a (1 - e cos E), |v_inertial|^2 =
    mu (2/r - 1/a), and r x v_inertial is constant. v_inertial =
    v_ecef + omega_e z x r (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1: r.v is the same
    in both frames)."""
    mu, we = MU[system], OMEGA_E[system]
    o = pure_kepler(system, sqrt_a, e, 0.96, 0.47, 1.0, -1.9)
    a = sqrt_a * sqrt_a
    h0 = None
    for dt in range(-7200, 7201, 900):
        s = orbits.kepler_state(o, 2295, 93600.0 + dt)
        x, y, z = s.position
        vx, vy, vz = s.velocity
        vix, viy, viz = vx - we * y, vy + we * x, vz
        r = math.sqrt(x * x + y * y + z * z)
        assert r == pytest.approx(a * (1 - e * math.cos(s.eccentric_anomaly)), abs=1e-6)
        assert vix * vix + viy * viy + viz * viz == pytest.approx(mu * (2 / r - 1 / a), rel=1e-12)
        # r x v_inertial, in rotating axes: its length and its z component
        # (along the rotation axis) are constant.
        h = (y * viz - z * viy, z * vix - x * viz, x * viy - y * vix)
        if h0 is None:
            h0 = h
        assert math.hypot(*h) == pytest.approx(math.hypot(*h0), rel=1e-12)
        assert h[2] == pytest.approx(h0[2], rel=1e-12)
        # the relativistic term in both forms of IS-GPS-200N with IRN-IS-200N-003 and -004
        # 20.3.3.3.3.1
        f = -2 * math.sqrt(mu) / C ** 2
        rel_e = f * e * sqrt_a * math.sin(s.eccentric_anomaly)
        rel_rv = -2 * (x * vx + y * vy + z * vz) / C ** 2
        assert rel_rv == pytest.approx(rel_e, abs=1e-16)


def test_geostationary_beidou_geo_stays_fixed():
    """BDS ICD Table 5-11 GEO algorithm on a circular orbit of mean motion
    omega_e, inclined by 5 degrees in the GEO frame with its node at 180
    degrees so that Rx(-5 deg) makes it equatorial: the result must be a
    geostationary satellite, fixed in BDCS, on the equator."""
    we = OMEGA_E[orbits.BEIDOU]
    a = (MU[orbits.BEIDOU] / (we * we)) ** (1 / 3)
    toe = 93600.0
    o = pure_kepler(orbits.BEIDOU, math.sqrt(a), 0.0, 5 * PI_ICD / 180,
                    PI_ICD + we * toe, 0.0, 0.3, geo=True, toe=toe, week=939)
    first = orbits.kepler_state(o, 939, toe).position
    assert abs(first[2]) < 1e-6
    assert math.hypot(*first) == pytest.approx(a, abs=1e-6)
    for dt in range(-43200, 43201, 3600):
        s = orbits.kepler_state(o, 939, toe + dt)
        assert math.dist(s.position, first) < 1e-5, dt
        assert math.hypot(*s.velocity) < 1e-9, dt


def test_week_crossover_is_the_actual_time_difference():
    """Table 20-IV note: tk is the actual total time difference, across the
    end of a week."""
    r = dict(RECORDS["G05"], toe=604000.0)
    o = kepler_orbit(orbits.GPS, r, 2295)
    after = orbits.kepler_state(o, 2296, 100.0)
    same = orbits.kepler_state(o, 2295, 604900.0)
    assert after == same
    before = orbits.kepler_state(o, 2295, 604000.0 - 900.0)
    assert math.dist(after.position, before.position) > 1e6


# ---- clock --------------------------------------------------------------------

def test_clock_polynomial_and_drift():
    c = KeplerClock(toc_week=2295, toc=604000.0, af0=1e-4, af1=-2e-11, af2=3e-18)
    dt = 900.0  # week 2296, 100 s: across the week end
    assert orbits.kepler_clock(c, 2296, 100.0) == pytest.approx(
        1e-4 - 2e-11 * dt + 3e-18 * dt * dt, abs=1e-20)
    assert orbits.kepler_clock_drift(c, 2296, 100.0) == pytest.approx(
        -2e-11 + 2 * 3e-18 * dt, abs=1e-24)
    h = 10.0
    numerical = (orbits.kepler_clock(c, 2296, 100.0 + h)
                 - orbits.kepler_clock(c, 2296, 100.0 - h)) / (2 * h)
    assert orbits.kepler_clock_drift(c, 2296, 100.0) == pytest.approx(numerical, abs=1e-19)
    assert orbits.kepler_clock(c, 2295, 604000.0) == 1e-4


# ---- public broadcast ephemerides against public precise orbits ---------------

@pytest.mark.parametrize("name", sorted(RECORDS))
def test_broadcast_orbit_against_precise_orbit(name):
    o, week_offset, offset, bound = public(name)
    for week, sow, x, y, z, _ in RECORDS[name]["sp3"]:
        s = orbits.kepler_state(o, week - week_offset, sow + offset)
        assert math.dist(s.position, (x, y, z)) < bound, (name, sow)


@pytest.mark.parametrize("name", sorted(RECORDS))
def test_velocity_is_the_derivative_of_position(name):
    """Table 20-IV sheet 3 velocity (and the GEO derivative) against a
    central difference of the positions, on the real orbits."""
    o, week_offset, offset, _ = public(name)
    week, sow = 2295 - week_offset, RECORDS[name]["toe"] + 1234.5
    h = 0.5
    s = orbits.kepler_state(o, week, sow)
    p1 = orbits.kepler_state(o, week, sow + h).position
    p0 = orbits.kepler_state(o, week, sow - h).position
    numerical = [(b - a) / (2 * h) for a, b in zip(p0, p1)]
    assert math.dist(s.velocity, numerical) < 1e-4


@pytest.mark.parametrize("name", sorted(RECORDS))
def test_relativistic_forms_agree_on_broadcast_orbits(name):
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1: F e sqrt(A) sin E and -2 r.v/c^2 are
    equivalent for a Kepler ellipse. A broadcast orbit adds the harmonic
    radius terms, whose radial rate is at most 2 n (|Crs| + |Crc|) (Table
    20-IV sheet 3), so the two forms differ by at most 2 r 2 n (|Crs| +
    |Crc|) / c^2, plus a little for delta n and the other terms."""
    o, week_offset, offset, _ = public(name)
    s = orbits.kepler_state(o, 2295 - week_offset, RECORDS[name]["toe"] + 1800.0)
    mu = MU[o.system]
    rel_e = -2 * math.sqrt(mu) / C ** 2 * o.e * o.sqrt_a * math.sin(s.eccentric_anomaly)
    rel_rv = -2 * sum(p * v for p, v in zip(s.position, s.velocity)) / C ** 2
    r = math.hypot(*s.position)
    n = math.sqrt(mu) / o.sqrt_a ** 3
    bound = 2 * r * 2 * n * (abs(o.crs) + abs(o.crc)) / C ** 2
    assert abs(rel_e - rel_rv) < 1.1 * bound + 1e-12


@pytest.mark.parametrize("name", ["G05", "E03"])
def test_broadcast_clock_against_precise_clock(name):
    """The precise clocks exclude the periodic relativistic term, as the
    broadcast polynomial does; both refer to the dual-frequency
    combination of the system's clock definition. GPS and Galileo agree
    within a few nanoseconds."""
    r = RECORDS[name]
    c = KeplerClock(2295, r["toc"], r["af0"], r["af1"], r["af2"])
    for week, sow, _, _, _, clock_us in r["sp3"]:
        assert abs(orbits.kepler_clock(c, week, sow) - clock_us * 1e-6) < 10e-9
