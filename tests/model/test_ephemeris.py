"""The ephemeris store (rules E1 to E5, DESIGN 2.5 validity) and the RTCM adapter.

Messages are synthetic: chosen integers encoded with ``gnsswire.rtcm3.encode``
and decoded back from the bytes, so every test starts from a frame as a
receiver or caster would send it. Expected values are the integers times the
resolutions of RTCM 10403.3 (and the ICDs they repeat), and the time-scale
definitions of the ICDs, computed here independently of the adapter.
"""

from __future__ import annotations

import dataclasses

import pytest
from gnsswire import rtcm3

from ssrcore.adapters import rtcm3 as adapter
from ssrcore.model.ephemeris import (
    DEPTH,
    EphemerisStore,
    beidou_iodc,
    ephemeris_message,
    gps_fit_interval,
    healthy,
    is_valid,
    validity,
)
from ssrcore.model.types import FeedOutcome, FeedReason, GnssTime, Missing, Sat, System
from ssrcore.toolkit.orbits import beidou_geo

WEEK = 2400                       # a GPS week of 2026
NEAR = GnssTime(WEEK, 302400.0)   # mid-week receipt time
HOUR = 3600.0


def wire(cls, **fields):
    """A message with every field 0 but those given, through bytes and back."""
    values = {f.name: 0 for f in dataclasses.fields(cls) if f.name != "extension"}
    values.update(number=int(cls.__name__[3:]), **fields)
    frame = rtcm3.encode(cls(**values))
    frames = list(rtcm3.scan(frame))
    assert len(frames) == 1
    return rtcm3.decode(frames[0])


def gps(satellite=5, iode=40, toe=302400.0, week=WEEK, iodc=None, **fields):
    """A 1019 with toe = toc given in GPS seconds of week (multiple of 16 s); IODC = IODE unless given."""
    return wire(rtcm3.Msg1019, satellite=satellite, week=week % 1024, iode=iode,
                iodc=iode if iodc is None else iodc,
                toe=int(toe // 16), toc=int(toe // 16), sqrt_a=2_702_000_000, **fields)


def galileo(cls=rtcm3.Msg1045, satellite=11, iodnav=100, toe=302400.0, **fields):
    return wire(cls, satellite=satellite, week=(WEEK - 1024) % 4096, iodnav=iodnav,
                toe=int(toe // 60), toc=int(toe // 60), **fields)


def beidou(satellite=20, toe_bdt=302400.0, week=WEEK - 1356, **fields):
    return wire(rtcm3.Msg1042, satellite=satellite, week=week % 8192,
                toe=int(toe_bdt // 8), toc=int(toe_bdt // 8), **fields)


def qzss(satellite=2, iode=40, toe=302400.0, week=WEEK, iodc=None, **fields):
    """A 1044 with toe = toc in GPS seconds of week; IODC = IODE unless given."""
    return wire(rtcm3.Msg1044, satellite=satellite, week=week % 1024, iode=iode,
                iodc=iode if iodc is None else iodc,
                toe=int(toe // 16), toc=int(toe // 16), sqrt_a=2_702_000_000, **fields)


def glonass(satellite=3, tb=60, channel=8, **fields):
    """A 1020; tb = 60 is 15:00 Moscow time, near NEAR (Wednesday 12:00 GPS = 14:59:42 Moscow)."""
    return wire(rtcm3.Msg1020, satellite=satellite, tb=tb, channel=channel, **fields)


def store_with(*messages, leap=18):
    store = EphemerisStore()
    for m in messages:
        result = store.add(adapter.ephemeris(m, NEAR, leap), received=NEAR)
        assert result.outcome == FeedOutcome.STORED, result
    return store


# ---- adapter: scaling and times (rule E3: the ICD's parameters in the ICD's units) ------


def test_rule_E3_1019_fields_in_icd_units():
    """RTCM Table 3.5-21 and IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-III: integer ×
    resolution, angles in semicircles."""
    m = gps(idot=-3, af2=2, af1=-5, af0=7, crs=-11, delta_n=13, m0=-17, cuc=19, e=23, cus=-29,
            cic=31, omega0=-37, cis=41, i0=43, crc=-47, omega=53, omega_dot=-59, tgd=-6,
            accuracy=2, fit_interval=1)
    e = adapter.from_1019(m, NEAR)
    assert e.sat == Sat.gps(5) and e.message == 1019 and e.iod == 40
    assert e.toe == e.toc == GnssTime(WEEK, 302400.0)
    expected = dict(sqrt_a=2_702_000_000 * 2.0**-19, e=23 * 2.0**-33, m0=-17 * 2.0**-31,
                    delta_n=13 * 2.0**-43, omega=53 * 2.0**-31, omega0=-37 * 2.0**-31,
                    omega_dot=-59 * 2.0**-43, i0=43 * 2.0**-31, idot=-3 * 2.0**-43,
                    cuc=19 * 2.0**-29, cus=-29 * 2.0**-29, crc=-47 * 2.0**-5, crs=-11 * 2.0**-5,
                    cic=31 * 2.0**-29, cis=41 * 2.0**-29, af0=7 * 2.0**-31, af1=-5 * 2.0**-43,
                    af2=2 * 2.0**-55)
    for name, value in expected.items():
        assert getattr(e, name) == value, name
    assert e.group_delays == (-6 * 2.0**-31,)
    assert e.accuracy == 2 and e.fit_interval == 1 and e.health == 0


def test_rule_E3_1044_fields_in_icd_units():
    """RTCM Table 3.5-112 (DF429 to DF457): the resolutions of GPS, the GPS week, the
    SPARTN QZSS number (DF429 1 is PRN 193) and the LEX bit of the health left out."""
    m = qzss(idot=-3, af2=2, af1=-5, af0=7, crs=-11, delta_n=13, m0=-17, cuc=19, e=23, cus=-29,
             cic=31, omega0=-37, cis=41, i0=43, crc=-47, omega=53, omega_dot=-59, tgd=-6,
             accuracy=2, fit_interval=1, health=1)
    e = adapter.ephemeris(m, NEAR)
    assert e.sat == Sat.qzss(2) and e.message == 1044 and e.iod == 40 and e.iodc == 40
    assert e.toe == e.toc == GnssTime(WEEK, 302400.0)
    assert (e.m0, e.crc, e.af0, e.af2, e.group_delays) == (
        -17 * 2.0**-31, -47 * 2.0**-5, 7 * 2.0**-31, 2 * 2.0**-55, (-6 * 2.0**-31,))
    assert e.health == 0 and e.fit_interval == 1
    assert adapter.from_1044(qzss(health=0b100000), NEAR).health != 0
    with pytest.raises(adapter.EphemerisError):
        adapter.from_1044(qzss(iode=40, iodc=41), NEAR)


def test_rule_E3_galileo_fields_and_gst_week():
    """RTCM Tables 3.5-110, 3.5-111: clock resolutions 2^-34, 2^-46, 2^-59; GST week + 1024 = GPS week."""
    m = galileo(af0=-9, af1=4, af2=-1, bgd_e5a_e1=-3, sisa=107)
    e = adapter.from_1045(m, NEAR)
    assert e.sat == Sat.galileo(11) and e.message == 1045 and e.iod == 100
    assert e.toe == GnssTime(WEEK, 302400.0)
    assert (e.af0, e.af1, e.af2) == (-9 * 2.0**-34, 4 * 2.0**-46, -1 * 2.0**-59)
    assert e.group_delays == (-3 * 2.0**-32,) and e.accuracy == 107
    i = adapter.from_1046(galileo(rtcm3.Msg1046, bgd_e5a_e1=2, bgd_e5b_e1=-2), NEAR)
    assert i.message == 1046 and i.group_delays == (2 * 2.0**-32, -2 * 2.0**-32)


def test_rule_E3_beidou_fields_and_bdt():
    """RTCM Table 3.5-113: harmonic terms 2^-31 rad, 2^-6 m, clock 2^-33, 2^-50, 2^-66, TGD 0.1 ns;
    BDT week 0 is GPS week 1356 and GPS − BDT = 14 s (BDS ICD 3.3, TM-23 3.1.1)."""
    m = beidou(cuc=3, crs=-5, a0=7, a1=-11, a2=13, tgd1=-25, tgd2=40, urai=3)
    e = adapter.from_1042(m, NEAR)
    assert e.sat == Sat.beidou(20) and e.message == 1042
    assert e.toe == e.toc == GnssTime(WEEK, 302414.0)  # BDT 302400 s + 14 s
    assert (e.cuc, e.crs) == (3 * 2.0**-31, -5 * 2.0**-6)
    assert (e.af0, e.af1, e.af2) == (7 * 2.0**-33, -11 * 2.0**-50, 13 * 2.0**-66)
    assert e.group_delays == (-25e-10, 40e-10) and e.accuracy == 3


def test_week_rollover_resolved_nearest_the_receipt_time():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1: the 10-bit week is the full week
    modulo 1024."""
    e = adapter.from_1019(gps(), NEAR)
    assert gps().week == WEEK - 2048 and e.toe.week == WEEK
    assert adapter.from_1019(gps(), GnssTime(1376, 302400.0)).toe.week == 1376  # one rollover earlier


def test_reference_time_across_the_end_of_week():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV, week crossover: a toe of 01:00
    broadcast late on Saturday belongs to the next week, and one of 23:00 Saturday received early
    Sunday to the last."""
    late_saturday = GnssTime(WEEK, 604800.0 - 1800.0)
    e = adapter.from_1019(gps(toe=3600.0), late_saturday)
    assert e.toe == GnssTime(WEEK + 1, 3600.0)
    early_sunday = GnssTime(WEEK + 1, 600.0)
    e = adapter.from_1019(gps(toe=604800.0 - 3600.0, week=WEEK + 1), early_sunday)
    assert e.toe == GnssTime(WEEK, 604800.0 - 3600.0)


def test_gps_cutover_iode_iodc_mismatch_refused():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.1: IODE and the 8 LSBs of IODC differ
    during a cutover."""
    m = wire(rtcm3.Msg1019, satellite=5, iode=40, iodc=256 + 41)
    with pytest.raises(adapter.EphemerisError):
        adapter.from_1019(m, NEAR)
    assert adapter.from_1019(wire(rtcm3.Msg1019, satellite=5, iode=40, iodc=256 + 40), NEAR).iod == 40


def test_satellite_id_zero_refused():
    """RTCM DF038: 0 means the slot is unknown; there is no satellite to key it by."""
    with pytest.raises(adapter.EphemerisError):
        adapter.from_1020(glonass(satellite=0), NEAR, 18)


def test_glonass_tb_on_the_gps_scale():
    """GLONASS ICD 3.3.3: GLONASS time = UTC + 3 h, GPS = UTC + leap seconds; tb in 15 min steps."""
    leap = 18
    midnight = GnssTime(WEEK, 3 * 86400.0)  # Wednesday 00:00:00 GPS = Tuesday 23:59:42 UTC
    # = Wednesday 02:59:42 Moscow; tb = 12 is 03:00 Moscow = 00:00 UTC = 00:00:18 GPS.
    e = adapter.from_1020(glonass(tb=12), midnight, leap)
    assert e.tb == midnight + 18.0
    assert e.iod == 12 and e.sat == Sat.glonass(3) and e.frequency_channel == 1
    # tb = 95 is 23:45 Moscow; received at 02:59 Moscow it belongs to the day before.
    e = adapter.from_1020(glonass(tb=95), midnight, leap)
    assert e.tb == midnight + 18.0 - 3 * HOUR - 15 * 60.0
    # The same message received just before midnight Moscow is of the same day.
    e = adapter.from_1020(glonass(tb=95), midnight - 3 * HOUR + 18.0 - 60.0, leap)
    assert e.tb == midnight + 18.0 - 3 * HOUR - 15 * 60.0


def test_glonass_fields_and_health():
    """RTCM Table 3.5-22: 2^-11 km, 2^-20 km/s, 2^-30 km/s², γ 2^-40, τ 2^-30 s; P1 Table 4.3."""
    m = glonass(x=-1000, y=2000, z=3, x_dot=-7, y_dot=0, z_dot=9, x_ddot=-1, z_ddot=2,
                tau=-123, gamma=45, delta_tau=-3, p1=2, bn_msb=0, m=1, ln_third=1)
    e = adapter.from_1020(m, NEAR, 18)
    assert e.position == (-1000 * 2.0**-11 * 1000, 2000 * 2.0**-11 * 1000, 3 * 2.0**-11 * 1000)
    assert e.velocity == (-7 * 2.0**-20 * 1000, 0.0, 9 * 2.0**-20 * 1000)
    assert e.acceleration == (-1 * 2.0**-30 * 1000, 0.0, 2 * 2.0**-30 * 1000)
    assert (e.tau, e.gamma, e.delta_tau) == (-123 * 2.0**-30, 45 * 2.0**-40, -3 * 2.0**-30)
    assert e.p1 == 45 * 60.0
    assert e.health == 2 and not healthy(e)          # ln = 1 on a GLONASS-M satellite
    assert healthy(adapter.from_1020(glonass(m=0, ln_third=1), NEAR, 18))  # ln invalid when not M
    assert not healthy(adapter.from_1020(glonass(bn_msb=1), NEAR, 18))


def test_glonass_without_leap_seconds_is_held():
    """Rule T4: no GLONASS time on the GPS scale without the leap seconds."""
    with pytest.raises(adapter.LeapSecondsNeeded):
        adapter.ephemeris(glonass(), NEAR)


def test_glonass_channel_outside_the_icd_refused():
    """GLONASS ICD 3.3.1.1: channels −7 to 6 (DF040 can carry up to +13)."""
    with pytest.raises(adapter.EphemerisError):
        adapter.from_1020(glonass(channel=14), NEAR, 18)


def test_not_an_ephemeris_message():
    """DESIGN 2.1: a type the engine does not know is a programming error."""
    with pytest.raises(TypeError):
        adapter.ephemeris(object(), NEAR)


def test_rule_E4_group_delays_kept_apart_from_the_clock():
    """Rule E4: the broadcast group delay is carried separately and not folded into af0,
    so the engine, which applies it in no case, leaves it out of every clock."""
    e = adapter.from_1019(gps(af0=1000, tgd=-20), NEAR)
    assert e.af0 == 1000 * 2.0**-31 and e.group_delays == (-20 * 2.0**-31,)


# ---- rule E1: issue of data --------------------------------------------------------------


def test_rule_E1_beidou_iodc_from_toc():
    """SPARTN SF100: BDS IODC = mod(toc/720, 240), toc in BDT seconds."""
    assert beidou_iodc(0) == 0
    assert beidou_iodc(3600) == 5
    assert beidou_iodc(719) == 0 and beidou_iodc(720) == 1
    assert beidou_iodc(240 * 720) == 0
    assert beidou_iodc(604792) == (604792 // 720) % 240 == 119
    e = adapter.from_1042(beidou(toe_bdt=7200.0), GnssTime(WEEK, 7200.0))
    assert e.iod == 10


def test_rule_E1_lookup_by_issue_of_data():
    """SPARTN 8.3, SF018: the correction applies to the ephemeris with its IODE."""
    store = store_with(gps(iode=40), gps(iode=41, toe=302400.0 + 2 * HOUR))
    found = store.find(Sat.gps(5), 41, 0, NEAR)
    assert found.miss is None and found.ephemeris.iod == 41
    assert store.find(Sat.gps(5), 40, 0, NEAR).ephemeris.iod == 40
    missing = store.find(Sat.gps(5), 42, 0, NEAR)
    assert missing.ephemeris is None and missing.miss == Missing.IODE_MISMATCH
    assert store.find(Sat.gps(6), 40, 0, NEAR).miss == Missing.NO_EPHEMERIS


def test_rule_E1_glonass_iod_is_tb_and_beidou_iod_is_iodc():
    """SF019 (DESIGN Q3, RTCM DF392) and SF100 name the ephemeris by tb and by IODC."""
    store = store_with(glonass(tb=60), beidou(satellite=20))
    assert store.find(Sat.glonass(3), 60, 0, NEAR).ephemeris.iod == 60
    iodc = beidou_iodc(302400)
    assert store.find(Sat.beidou(20), iodc, 0, NEAR).ephemeris.sat == Sat.beidou(20)


# ---- rule E2: ephemeris type -------------------------------------------------------------


def test_rule_E2_galileo_fnav_and_inav():
    """SF096: 0 is F/NAV (1045), 1 is I/NAV (1046); the same IODnav in both is two ephemerides."""
    store = store_with(galileo(rtcm3.Msg1045, af0=1), galileo(rtcm3.Msg1046, af0=2))
    fnav = store.find(Sat.galileo(11), 100, 0, NEAR).ephemeris
    inav = store.find(Sat.galileo(11), 100, 1, NEAR).ephemeris
    assert (fnav.message, inav.message) == (1045, 1046) and fnav.af0 != inav.af0
    only_fnav = store_with(galileo(rtcm3.Msg1045))
    assert only_fnav.find(Sat.galileo(11), 100, 1, NEAR).miss == Missing.EPHEMERIS_TYPE
    assert only_fnav.find(Sat.galileo(11), 100, 2, NEAR).miss == Missing.EPHEMERIS_TYPE  # C/NAV


def test_rule_E2_beidou_d1_and_d2():
    """SF097 0 is D1 and 1 is D2 on B1I; BDS ICD 5.1.1: GEO satellites broadcast D2."""
    assert beidou_geo(5) and beidou_geo(59) and not beidou_geo(20)
    assert ephemeris_message(Sat.beidou(20), 0) == 1042 and ephemeris_message(Sat.beidou(20), 1) is None
    assert ephemeris_message(Sat.beidou(3), 1) == 1042 and ephemeris_message(Sat.beidou(3), 0) is None
    assert ephemeris_message(Sat.beidou(20), 2) is None  # B3I: no RTCM message here
    store = store_with(beidou(satellite=3))
    iodc = beidou_iodc(302400)
    assert store.find(Sat.beidou(3), iodc, 1, NEAR).ephemeris is not None
    assert store.find(Sat.beidou(3), iodc, 0, NEAR).miss == Missing.EPHEMERIS_TYPE


def test_rule_E2_types_of_gps_glonass_qzss():
    """SF016, SF017, SF098: only type 0 is defined (L1C/A, L1C/A, LNAV); 1 to 3 TBD."""
    assert ephemeris_message(Sat.gps(1), 0) == 1019 and ephemeris_message(Sat.gps(1), 1) is None
    assert ephemeris_message(Sat.glonass(1), 0) == 1020 and ephemeris_message(Sat.glonass(1), 3) is None
    assert ephemeris_message(Sat.qzss(1), 0) == 1044 and ephemeris_message(Sat.qzss(1), 2) is None
    store = store_with(gps())
    assert store.find(Sat.gps(5), 40, 1, NEAR).miss == Missing.EPHEMERIS_TYPE


# ---- rule E5: health ----------------------------------------------------------------------


def test_rule_E5_unhealthy_stored_not_used():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.4: a non-zero health; stored with the
    reason, never found."""
    store = EphemerisStore()
    result = store.add(adapter.from_1019(gps(health=0b100000), NEAR), NEAR)
    assert result.outcome == FeedOutcome.STORED and result.reason == FeedReason.UNHEALTHY
    assert len(store.issues(Sat.gps(5))) == 1
    assert store.find(Sat.gps(5), 40, 0, NEAR).miss == Missing.UNHEALTHY
    assert store.broadcast(Sat.gps(5), NEAR) is None


def test_rule_E5_health_of_each_system():
    """Galileo OS SIS ICD 2.2 5.1.9.3 (SHS and DVS), BDS ICD 5.2.4.6 (SatH1): anything but 0 is not nominal."""
    assert healthy(adapter.from_1045(galileo(), NEAR))
    assert not healthy(adapter.from_1045(galileo(e5a_validity=1), NEAR))
    assert not healthy(adapter.from_1045(galileo(e5a_health=2), NEAR))
    assert not healthy(adapter.from_1046(galileo(rtcm3.Msg1046, e1b_health=1), NEAR))
    assert not healthy(adapter.from_1046(galileo(rtcm3.Msg1046, e5b_validity=1), NEAR))
    assert healthy(adapter.from_1042(beidou(), NEAR))
    assert not healthy(adapter.from_1042(beidou(health=1), NEAR))


# ---- issues kept --------------------------------------------------------------------------


def test_three_issues_kept():
    """DESIGN 3.2: current, previous and next; a fourth drops the oldest reference time."""
    store = EphemerisStore()
    assert DEPTH == 3
    for k, iode in enumerate((40, 41, 42)):
        e = adapter.from_1019(gps(iode=iode, toe=302400.0 + 2 * HOUR * k), NEAR)
        assert store.add(e, NEAR).outcome == FeedOutcome.STORED
    newest = adapter.from_1019(gps(iode=43, toe=302400.0 + 6 * HOUR), NEAR + 6 * HOUR)
    assert store.add(newest, NEAR + 6 * HOUR).outcome == FeedOutcome.STORED
    assert [s.ephemeris.iod for s in store.issues(Sat.gps(5))] == [41, 42, 43]
    assert store.find(Sat.gps(5), 40, 0, NEAR).miss == Missing.IODE_MISMATCH
    # An issue older than the three kept is not stored.
    oldest = adapter.from_1019(gps(iode=39, toe=302400.0 - 2 * HOUR), NEAR)
    result = store.add(oldest, NEAR)
    assert result.outcome == FeedOutcome.SUPERSEDED and result.reason == FeedReason.OLDER
    assert [s.ephemeris.iod for s in store.issues(Sat.gps(5))] == [41, 42, 43]


def test_duplicate_keeps_first_receipt():
    store = EphemerisStore()
    e = adapter.from_1019(gps(), NEAR)
    store.add(e, NEAR)
    assert store.add(e, NEAR + 30.0).outcome == FeedOutcome.DUPLICATE
    (only,) = store.issues(Sat.gps(5))
    assert only.received == NEAR


def test_iod_reused_with_a_later_reference_time_replaces():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.4.4: an IODE is not repeated within six hours,
    so a repeated IOD with a later toe is a new data set; with an earlier toe it is an old one."""
    store = EphemerisStore()
    store.add(adapter.from_1019(gps(iode=40), NEAR), NEAR)
    later = adapter.from_1019(gps(iode=40, toe=302400.0 + 8 * HOUR), NEAR + 8 * HOUR)
    assert store.add(later, NEAR + 8 * HOUR).outcome == FeedOutcome.STORED
    (only,) = store.issues(Sat.gps(5))
    assert only.ephemeris.toe == later.toe
    result = store.add(adapter.from_1019(gps(iode=40), NEAR), NEAR)
    assert result.outcome == FeedOutcome.SUPERSEDED and result.reason == FeedReason.OLDER


def test_same_iod_and_time_with_another_body_is_a_conflict():
    store = EphemerisStore()
    store.add(adapter.from_1019(gps(af0=1), NEAR), NEAR)
    result = store.add(adapter.from_1019(gps(af0=2), NEAR), NEAR + 1.0)
    assert result.outcome == FeedOutcome.STORED and result.reason == FeedReason.CONFLICT
    assert store.find(Sat.gps(5), 40, 0, NEAR).ephemeris.af0 == 2 * 2.0**-31


def test_broadcast_without_corrections_takes_the_nearest_valid_issue():
    store = store_with(gps(iode=40), gps(iode=41, toe=302400.0 + 2 * HOUR))
    assert store.broadcast(Sat.gps(5), NEAR + 0.5 * HOUR).ephemeris.iod == 40
    assert store.broadcast(Sat.gps(5), NEAR + 1.5 * HOUR).ephemeris.iod == 41
    assert store.broadcast(Sat.gps(5), NEAR + 1.0 * HOUR).ephemeris.iod == 41   # tie: the later
    assert store.broadcast(Sat.gps(5), NEAR + 5.0 * HOUR) is None


def test_glonass_channel_from_the_ephemeris():
    """DESIGN 2.6.3: the frequency channel comes from the GLONASS ephemerides."""
    store = store_with(glonass(channel=3))
    assert store.glonass_channel(Sat.glonass(3)) == -4
    assert store.glonass_channel(Sat.glonass(4)) is None
    assert store.satellites() == (Sat.glonass(3),)


# ---- validity table of DESIGN 2.5 ---------------------------------------------------------


@pytest.mark.parametrize("message, half", [
    (lambda: adapter.from_1019(gps(), NEAR), 2 * HOUR),                    # fit flag 0: 4 h fit
    (lambda: adapter.from_1019(gps(fit_interval=1), NEAR), 3 * HOUR),      # flag 1: at least 6 h
    (lambda: adapter.from_1045(galileo(), NEAR), 4 * HOUR),
    (lambda: adapter.from_1042(beidou(), NEAR), 2 * HOUR),
    (lambda: adapter.from_1020(glonass(p1=0), NEAR, 18), 15 * 60.0),
    (lambda: adapter.from_1020(glonass(p1=1), NEAR, 18), 30 * 60.0),       # P1 = 30 min
    (lambda: adapter.from_1020(glonass(p1=3), NEAR, 18), 45 * 60.0),       # P1 = 60 min
    (lambda: adapter.from_1044(qzss(), NEAR), 1 * HOUR),                   # DF457 0: 2 h fit
    (lambda: adapter.from_1044(qzss(fit_interval=1), NEAR), 1 * HOUR),     # DF457 1: over 2 h
], ids=["gps-fit0", "gps-fit1", "galileo", "beidou", "glonass-p1-0", "glonass-p1-30", "glonass-p1-60",
        "qzss-fit0", "qzss-fit1"])
def test_validity_table(message, half):
    """DESIGN 2.5: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1 and Table 20-XIII for GPS;
    design values for the others."""
    e = message()
    ref = e.tb if e.sat.system == System.GLONASS else e.toe
    assert validity(e) == half
    assert is_valid(e, ref + half) and is_valid(e, ref - half)
    assert not is_valid(e, ref + half + 1.0) and not is_valid(e, ref - half - 1.0)


@pytest.mark.parametrize("flag, iodc, hours", [
    (0, 250, 4), (1, None, 6), (1, 40, 6), (1, 256 + 40, 6), (1, 239, 6), (1, 240, 8), (1, 247, 8),
    (1, 248, 14), (1, 255, 14), (1, 496, 14), (1, 497, 26), (1, 503, 26), (1, 504, 6),
    (1, 1020, 6), (1, 1021, 26), (1, 1023, 26),
])
def test_gps_fit_interval_by_iodc(flag, iodc, hours):
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1 and Table 20-XII: the fit's length
    from the flag and IODC."""
    assert gps_fit_interval(flag, iodc) == hours * HOUR


def test_gps_validity_uses_the_iodc_of_the_1019():
    """DESIGN 2.5: a 1019 with flag 1 and IODC 497 (a 26 h fit) is valid 13 h either side of toe."""
    e = adapter.from_1019(gps(fit_interval=1, iode=497 & 0xFF, iodc=497), NEAR)
    assert e.iodc == 497 and validity(e) == 13 * HOUR


def test_lookup_outside_validity():
    store = store_with(gps())
    assert store.find(Sat.gps(5), 40, 0, NEAR + 2 * HOUR).miss is None
    assert store.find(Sat.gps(5), 40, 0, NEAR + 2 * HOUR + 1.0).miss == Missing.OUTSIDE_VALIDITY


# ---- live service: issues chosen by number while corrections lag -----------------------


def test_live_issue_every_ten_minutes_with_corrections_90_s_behind():
    """DESIGN 3.2 and rule E1, live-service fact (DESIGN 5.1): one constellation
    changes issue every 10 minutes and the corrections name a new issue about 90 s
    after it first appears. The issue in use and the newer one are both held and
    chosen by the issue number the corrections name, not by recency."""
    store = EphemerisStore()
    sat = Sat.galileo(11)
    for k in range(6):
        arrival = NEAR + 600.0 * k
        e = adapter.from_1045(galileo(iodnav=100 + k, toe=302400.0 + 600.0 * k), arrival)
        assert store.add(e, arrival).outcome == FeedOutcome.STORED
        if k:
            for lag in (0.0, 45.0, 89.0):         # the corrections still name the previous issue
                found = store.find(sat, 100 + k - 1, 0, arrival + lag)
                assert found.miss is None and found.ephemeris.iod == 100 + k - 1
        found = store.find(sat, 100 + k, 0, arrival + 90.0)
        assert found.miss is None and found.ephemeris.iod == 100 + k
        assert len(store.issues(sat)) == min(k + 1, DEPTH)
