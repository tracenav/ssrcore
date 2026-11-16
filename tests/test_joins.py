"""Joins between items: real code on both sides of each interface.

Each test feeds one layer's real output to the next layer's real code, so a
stand-in on either side cannot hide a disagreement of units, time scale or
convention (the toolkit's time scales, orbits and astronomy with the model's stores).
"""

from __future__ import annotations

import dataclasses
import math

import pytest
from gnsswire import rtcm3

from ssrcore.adapters import rtcm3 as adapter
from ssrcore.model.atmosphere import AreaChoice, AtmosphereStore, TroposphereZenith
from ssrcore.model.corrections import Corrections, Selection
from ssrcore.model.ephemeris import EphemerisStore, broadcast_state
from ssrcore.model.stream import Streams
from ssrcore.model.timebase import TimeBase
from ssrcore.model.types import (
    AreaDefinition,
    ClockCorrection,
    FeedOutcome,
    GadInput,
    GnssTime,
    HpacArea,
    HpacInput,
    MaxAge,
    OcbInput,
    OcbSatellite,
    OrbitCorrection,
    Sat,
    SpartnTag,
    StreamId,
    System,
    TroposphereBlock,
)
from ssrcore.toolkit import astro, orbits, relativity, tides, windup

PI_ICD = 3.1415926535898   # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2
LEAP = 18                  # GPS − UTC on 2024-01-01 (USNO / IERS Bulletin C)


def wire(cls, **fields):
    """A message with every field 0 but those given, through bytes and back."""
    values = {f.name: 0 for f in dataclasses.fields(cls) if f.name != "extension"}
    values.update(number=int(cls.__name__[3:]), **fields)
    frames = list(rtcm3.scan(rtcm3.encode(cls(**values))))
    assert len(frames) == 1
    return rtcm3.decode(frames[0])


def q(value, exponent, semicircle=False):
    """The broadcast integer of a value at resolution 2^exponent (angles from rad)."""
    return round((value / PI_ICD if semicircle else value) / 2.0**exponent)


# ---- 2c with 3b: RTCM ephemeris -> adapter -> store -> orbit kernels -------------
# Public data of 2024-01-01 (IGS/BKG BRDC00WRD_S_20240010000_01D_MN.rnx and GFZ
# GBM0MGXRAP_20240010000_01D_05M_ORB.SP3), the same records as
# tests/kernels/test_kepler.py and test_glonass.py: broadcast values in RINEX
# units (BeiDou toe, toc in BDT, week 939) and precise positions as (GPS week,
# GPS seconds, x, y, z in m, clock in microseconds).

G05 = dict(
    toc=93600.0, af0=-1.562642864883e-04, af1=-1.477928890381e-12, af2=0.0,
    toe=93600.0, sqrt_a=5153.783403397, e=0.005790898343548, m0=-2.494229256397,
    dn=4.453756945425e-09, omega=1.235180510432, omega0=-0.6299845609285,
    omega_dot=-8.155696860368e-09, i0=0.9693533956084, idot=-3.928735076157e-11,
    cuc=-1.747161149979e-06, cus=3.31737101078e-06, crc=319.34375, crs=-34.125,
    cic=1.098960638046e-07, cis=4.097819328308e-08,
    sp3=[(2295, 93600.0, -10165892.577, -12970215.216, -20987562.731, -156.263473),
         (2295, 97200.0, -7459642.749, -20779562.463, -14875907.434, -156.268086)])
C59 = dict(  # GEO (BDS-3)
    toc=93600.0, af0=-7.369089871645e-07, af1=1.767475055203e-13, af2=0.0,
    toe=93600.0, sqrt_a=6493.454147339, e=0.0001097865169868, m0=2.559401345431,
    dn=-1.592923494515e-09, omega=2.898151811953, omega0=-2.470151406647,
    omega_dot=2.769758228691e-09, i0=0.05742849059572, idot=-3.907305612105e-10,
    cuc=-1.462176442146e-06, cus=1.59814953804e-05, crc=-487.453125, crs=-37.40625,
    cic=-7.962808012962e-08, cis=6.658956408501e-08,
    sp3=[(2295, 93600.0, -32345709.692, 27048687.513, 582129.268, -0.901325),
         (2295, 97200.0, -32344236.344, 27055931.745, 256084.721, -0.900894)])
R01 = dict(  # tb 2024-01-02 00:15 UTC = 03:15 Moscow time, tb index 13
    position=(-11175834.47266, -9508662.597656, -20866544.92188),
    velocity=(2857.489585876, -39.36195373535, -1511.75403595),
    lunisolar=(-2.793967723846e-06, -1.862645149231e-06, 0.0),
    tau_n=-7.346831262112e-05, gamma_n=9.094947017729e-13,
    sp3=[(2295, 86700.0, -12899555.791, -9528742.824, -19838102.219, 73.473071),
         (2295, 87300.0, -11227238.908, -9507992.587, -20839253.070, 73.473228)])


def angles(r, harmonic, exp_c, exp_r):
    """The orbit fields common to RTCM 1019 and 1042 (resolutions of Tables 3.5-21, 3.5-113)."""
    return dict(
        sqrt_a=q(r["sqrt_a"], -19), e=q(r["e"], -33), m0=q(r["m0"], -31, True),
        delta_n=q(r["dn"], -43, True), omega=q(r["omega"], -31, True),
        omega0=q(r["omega0"], -31, True), omega_dot=q(r["omega_dot"], -43, True),
        i0=q(r["i0"], -31, True), idot=q(r["idot"], -43, True),
        cuc=q(r["cuc"], harmonic), cus=q(r["cus"], harmonic), cic=q(r["cic"], harmonic),
        cis=q(r["cis"], harmonic), crc=q(r["crc"], exp_r), crs=q(r["crs"], exp_r))


def through_store(message, near, sat, iod, ephemeris_type):
    store = EphemerisStore()
    assert store.add(adapter.ephemeris(message, near, LEAP), received=near).outcome == FeedOutcome.STORED
    eph = store.find(sat, iod, ephemeris_type, near).ephemeris
    assert eph is not None
    return eph


def check(eph, sp3, bound, clock_bound):
    for week, sow, x, y, z, clock_us in sp3:
        state, clock = broadcast_state(eph, GnssTime(week, sow))
        assert math.dist(state.position, (x, y, z)) < bound, sow
        if clock_bound is None:   # BeiDou: broadcast and precise clocks have different references
            dt = GnssTime(week, sow) - eph.toc
            assert clock == pytest.approx(eph.af0 + eph.af1 * dt, abs=1e-15)
        else:
            assert abs(clock - clock_us * 1e-6) < clock_bound, sow


def test_gps_1019_through_adapter_and_store_matches_precise_orbit():
    r = G05
    m = wire(rtcm3.Msg1019, satellite=5, week=2295 % 1024, iode=33, iodc=33,
             toe=int(r["toe"] // 16), toc=int(r["toc"] // 16),
             af0=q(r["af0"], -31), af1=q(r["af1"], -43), **angles(r, -29, -29, -5))
    eph = through_store(m, GnssTime(2295, 93600.0), Sat.gps(5), 33, 0)
    check(eph, r["sp3"], 3.0, 10e-9)


def test_beidou_geo_1042_in_bdt_through_adapter_and_store_matches_precise_orbit():
    """toe goes BDT -> GPS in the adapter and back to BDT for the kernel: a
    node term omega_e x toe on the GPS scale (14 s later) would be ~40 km off."""
    r = C59
    m = wire(rtcm3.Msg1042, satellite=59, week=939, toe=int(r["toe"] // 8),
             toc=int(r["toc"] // 8), a0=q(r["af0"], -33), a1=q(r["af1"], -50),
             **angles(r, -31, -31, -6))
    near = GnssTime(2295, 93600.0 + 14.0)
    eph = through_store(m, near, Sat.beidou(59), (93600 // 720) % 240, 1)
    assert eph.toe == near
    assert orbits.beidou_geo(59)
    check(eph, r["sp3"], 5.0, None)


def test_glonass_1020_placed_by_timescale_matches_precise_orbit():
    """tb = 13 (03:15 Moscow time) placed through ``timescale`` on GPS time."""
    r = R01
    km = 1000.0
    m = wire(rtcm3.Msg1020, satellite=1, channel=1 + 7, tb=13,
             x=q(r["position"][0] / km, -11), y=q(r["position"][1] / km, -11),
             z=q(r["position"][2] / km, -11),
             x_dot=q(r["velocity"][0] / km, -20), y_dot=q(r["velocity"][1] / km, -20),
             z_dot=q(r["velocity"][2] / km, -20),
             x_ddot=q(r["lunisolar"][0] / km, -30), y_ddot=q(r["lunisolar"][1] / km, -30),
             z_ddot=q(r["lunisolar"][2] / km, -30),
             tau=q(r["tau_n"], -30), gamma=q(r["gamma_n"], -40))
    near = GnssTime(2295, 86400.0 + 600.0)
    eph = through_store(m, near, Sat.glonass(1), 13, 0)
    assert eph.tb == GnssTime(2295, 86400.0 + 900.0 + LEAP)
    check(eph, r["sp3"], 5.0, 50e-9)


@pytest.mark.parametrize("near_offset", [-40000.0, 40000.0])
def test_glonass_tb_day_is_the_one_nearest_the_receipt_time(near_offset):
    """Across a GLONASS day boundary the adapter chooses the tb nearest ``near``."""
    tb_gps = GnssTime(2295, 86400.0 + 900.0 + LEAP)
    m = wire(rtcm3.Msg1020, satellite=1, channel=8, tb=13)
    eph = adapter.from_1020(m, tb_gps + near_offset, LEAP)
    assert eph.tb == tb_gps


# ---- 3a with 3c and 3d: time base -> one stream store -> corrections and atmosphere ----

STREAM = StreamId(1, 0)
TAG = 500_000_000


def ocb(tag, siou=3):
    block = OcbSatellite(Sat.gps(5), False, 320.0, OrbitCorrection(33, 0.1, -0.2, 0.05),
                         ClockCorrection(320.0, 0.3, 2))
    return OcbInput(System.GPS, STREAM, SpartnTag(tag, True), siou, True, False, 0, 0, (block,))


def gad(tag, siou=3, aiou=1):
    return GadInput(STREAM, SpartnTag(tag, True), siou, aiou,
                    (AreaDefinition(1, 50.0, 10.0, 3, 4, 1.0, 1.5),))


def hpac(tag, siou=3, aiou=1):
    tropo = TroposphereBlock(equation_type=0, quality=2, hydrostatic=2.3, coefficients=(0.1,))
    return HpacInput(System.GPS, STREAM, SpartnTag(tag, True), siou, aiou, (HpacArea(1, 0, tropo),))


def test_one_stream_store_serves_corrections_and_atmosphere():
    """The SIOU issue a selection names is the issue the atmosphere was stored
    under, because both stores note their messages in one ``Streams``; an OCB
    gap of 510 s starts a new SIOU era for the atmosphere as well (S3, S4)."""
    time, streams = TimeBase(), Streams()
    corrections, atmosphere = Corrections(streams), AtmosphereStore(streams, 120.0)

    def at(tag):
        t = time.resolve(SpartnTag(tag, True), System.GPS)
        assert isinstance(t, GnssTime)
        return t

    atmosphere.feed_gad(gad(TAG), at(TAG))
    atmosphere.feed_hpac(hpac(TAG + 1), at(TAG + 1))
    corrections.feed(ocb(TAG + 2), at(TAG + 2))
    t = at(TAG + 3)
    selection = corrections.select(STREAM, Sat.gps(5), t, MaxAge())
    assert isinstance(selection, Selection)
    choice = atmosphere.area(STREAM, selection.issue, 49.0, 12.0, t)
    assert isinstance(choice, AreaChoice) and choice.issue.siou == selection.issue
    assert isinstance(atmosphere.troposphere(choice, t), TroposphereZenith)

    # 520 s without any message: the next OCB starts a new SIOU era, and the
    # atmosphere stored before it no longer combines with the new selection.
    corrections.feed(ocb(TAG + 523), at(TAG + 523))
    t = at(TAG + 524)
    later = corrections.select(STREAM, Sat.gps(5), t, MaxAge())
    assert isinstance(later, Selection) and later.issue.era == selection.issue.era + 1
    assert not isinstance(atmosphere.area(STREAM, later.issue, 49.0, 12.0, t), AreaChoice)



# ---- 2g with 2h: one evaluation of the Sun and Moon for the tide and the wind-up ------


def test_sky_feeds_the_tide_and_the_sun_steered_frame():
    """``astro.sky`` is computed once per epoch; the tide and the TM-23 wind-up
    frame take it as an argument and agree with the kernels fed directly."""
    week, sow, ut1_minus_gps = 2310, 345600.0, -18.0
    s = astro.sky(week, sow, ut1_minus_gps)
    site = (4075578.385, 931852.890, 4801570.154)
    assert tides.solid_earth_tide(site, s) == tides.displacement(site, s.sun, s.moon, s.tt_centuries, s.gmst)
    sat = (15600e3, 7540e3, 20140e3)
    i, j, k = windup.satellite_frame_sun(sat, s.sun)
    to_sun = [a - b for a, b in zip(s.sun, sat)]
    # TM-23 (5.77): j is normal to the plane of the Earth's centre, satellite and Sun.
    assert abs(sum(a * b for a, b in zip(j, to_sun))) / math.dist(s.sun, sat) < 1e-12
    assert abs(sum(a * b for a, b in zip(j, k))) < 1e-12


# ---- 2c with 2i: the relativistic clock term from the orbit kernel's state ------------


def test_relativity_from_the_kepler_state_both_forms():
    """The two forms of IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1 on the state
    ``broadcast_state`` gives: -2 r.v/c^2 and F e sqrt(A) sin E agree within what the harmonic
    radius terms allow."""
    r = G05
    m = wire(rtcm3.Msg1019, satellite=5, week=2295 % 1024, iode=33, iodc=33,
             toe=int(r["toe"] // 16), toc=int(r["toc"] // 16), **angles(r, -29, -29, -5))
    eph = through_store(m, GnssTime(2295, 93600.0), Sat.gps(5), 33, 0)
    for dt in (0.0, 3600.0, 7200.0):
        state, _ = broadcast_state(eph, GnssTime(2295, 93600.0 + dt))
        a = relativity.clock(state.position, state.velocity)
        b = relativity.clock_kepler(3.986005e14, eph.e, eph.sqrt_a, state.eccentric_anomaly)
        # radial rate of the harmonic terms at most 2 n (|Crs| + |Crc|) (Table 20-IV
        # sheet 3), as tests/kernels/test_kepler.py bounds it
        n = math.sqrt(3.986005e14) / eph.sqrt_a ** 3
        bound = 2 * math.hypot(*state.position) * 2 * n * (abs(eph.crs) + abs(eph.crc)) / 299792458.0 ** 2
        assert abs(a - b) < 1.1 * bound + 1e-12 and abs(b) > 1e-9
