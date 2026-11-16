"""The time base (rules T1 to T5 of DESIGN 3.3) and the NMEA and UBX time adapters.

The time base uses ``toolkit.timescale``; the expected values here
are computed independently, from the ICD definitions cited and the proleptic
Gregorian calendar of Python's ``datetime``. Messages are synthetic, encoded
with gnsswire, or the worked examples of the u-blox F9 HPG 1.51 Interface
Description, which the u-blox X20 HPG 2.11 one repeats.
"""

from __future__ import annotations

import math
from datetime import date

import pytest
from gnsswire import nmea, ubx

from ssrcore.adapters import nmea as nmea_adapter
from ssrcore.adapters import ubx as ubx_adapter
from ssrcore.model.timebase import TimeBase
from ssrcore.model.types import (
    FeedReason,
    GnssTime,
    LeapSeconds,
    RoverPosition,
    SpartnTag,
    System,
    UtcTime,
)

WEEK = 604800.0
DAY = 86400.0
# SPARTN Table 7.1 TF009: the 16-bit tag runs from 0 to 43,199 s, one half day.
HALF_DAY = 43200.0
# IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4: GPS time starts at midnight 1980-01-06 (UTC, then
# equal to GPS).
GPS_EPOCH = date(1980, 1, 6).toordinal()
# SPARTN Table 7.1 TF009: the 32-bit tag counts from 2010-01-01 00:00:00.
SPARTN_EPOCH_DAYS = date(2010, 1, 1).toordinal() - GPS_EPOCH


def seconds(t: GnssTime) -> float:
    return t.week * WEEK + t.seconds


def gps(total: float) -> GnssTime:
    return GnssTime.normalised(0, total)


def utc(year, month, day, hour=0, minute=0, second=0.0, leap=0) -> GnssTime:
    """The GPS time of a UTC reading under the GPS − UTC count ``leap``."""
    days = date(year, month, day).toordinal() - GPS_EPOCH
    return gps(days * DAY + hour * 3600 + minute * 60 + second + leap)


def base(**kw) -> TimeBase:
    return TimeBase(**kw)


# A leap second: 2016-12-31 23:59:60 UTC, GPS − UTC 17 before and 18 after (IERS
# Bulletin C 52; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 gives the form of the
# announcement).
LEAP_2017 = LeapSeconds(17, 18, utc(2017, 1, 1, leap=18))


# ---- the join with toolkit.timescale ----------------------------------------------------


def test_spartn_epoch_through_the_toolkit():
    """SPARTN TF009: 2010-01-01 00:00:00 is GPS week 1564, Friday, 432000 s."""
    assert base().spartn_to_gps(System.GPS, 0) == GnssTime(1564, 432000.0)
    assert LEAP_2017.effective == GnssTime(1930, 18.0)


# ---- T1, T3: 32-bit tags in each constellation's time scale ---------------------------


def test_rule_T3_32_bit_tag_counts_from_2010():
    tb = base()
    assert tb.resolve(SpartnTag(0, True), System.GPS) == GnssTime(1564, 432000.0)
    assert tb.resolve(SpartnTag(86400 * 2 + 5, True), System.GPS) == GnssTime(1565, 5.0)   # Sunday


def test_rule_T1_tag_in_the_scale_of_the_constellation():
    tb = base()
    tag = SpartnTag(220_000_000, True)
    t_gps = tb.resolve(tag, System.GPS)
    assert tb.resolve(tag, System.GALILEO) == t_gps
    assert tb.resolve(tag, System.QZSS) == t_gps
    assert tb.resolve(tag, System.BEIDOU) - t_gps == 14.0   # BDT = GPS − 14 s
    # GAD and multi-constellation messages: GPS, the default.
    assert tb.resolve(tag) == t_gps


def test_rule_T3_glonass_epoch_in_glonass_time():
    """Q2: 2010-01-01 00:00:00 GLONASS time is 2009-12-31 21:00:00 UTC; GPS − UTC was 15 s."""
    tb = base()
    tb.set_leap_seconds(LeapSeconds(15))
    assert tb.resolve(SpartnTag(0, True), System.GLONASS) == utc(2009, 12, 31, 21, leap=15)
    assert tb.resolve(SpartnTag(10800, True), System.GLONASS) == utc(2010, 1, 1, leap=15)


# ---- T2: 16-bit tags near half-day edges -------------------------------------------------


@pytest.mark.parametrize("received_offset, tag, expected_offset", [
    (+2.0, 43199, -1.0),      # just after the half-day edge, a tag of the half day before
    (-2.0, 0, 0.0),           # just before the edge, a tag of the next half day
    (+0.5, 0, 0.0),
    (-0.5, 43199, -1.0),
    (+21599.0, 0, 0.0),       # 6 h less 1 s after the edge: still that half day
    (+21601.0, 0, 43200.0),   # past 6 h: the next half day is nearer
    (+21600.0, 0, 0.0),       # exactly 6 h: the earlier half day
])
def test_rule_T2_half_day_nearest_received(received_offset, tag, expected_offset):
    tb = base()
    edge = GnssTime(2100, 3 * DAY + HALF_DAY)   # a half-day edge of GPS time
    received = edge + received_offset
    assert tb.resolve(SpartnTag(tag, False), System.GPS, received) == edge + expected_offset


def test_rule_T2_half_day_across_a_week():
    tb = base()
    received = GnssTime(2100, 1.0)          # just after the start of a GPS week
    assert tb.resolve(SpartnTag(43199, False), System.GPS, received) == GnssTime(2099, WEEK - 1.0)


def test_rule_T2_half_day_of_beidou_time():
    """The half day is that of BDT, 14 s behind GPS time: tag 43190 BDT 11:59:50."""
    tb = base()
    received = utc(2024, 3, 1, 12, 0, 5)               # GPS 12:00:05, BDT 11:59:51
    t = tb.resolve(SpartnTag(43190, False), System.BEIDOU, received)
    assert t == utc(2024, 3, 1, 12, 0, 4)


def test_rule_T2_glonass_half_day_is_moscow():
    """GLONASS half days start at 00:00 and 12:00 GLONASS time, 21:00 and 09:00 UTC."""
    tb = base()
    tb.set_leap_seconds(LeapSeconds(18))
    received = utc(2024, 3, 1, 9, 0, 30, leap=18)
    assert tb.resolve(SpartnTag(0, False), System.GLONASS, received) == utc(2024, 3, 1, 9, leap=18)
    assert tb.resolve(SpartnTag(43170, False), System.GLONASS, received) == \
        utc(2024, 3, 1, 8, 59, 30, leap=18)


def test_rule_T2_uses_current_time_without_receipt_time():
    tb = base()
    tb.advance(GnssTime(2100, HALF_DAY + 3.0))
    assert tb.resolve(SpartnTag(43198, False)) == GnssTime(2100, HALF_DAY - 2.0)


# ---- T4: GLONASS and leap seconds ----------------------------------------------------------


def test_rule_T4_glonass_needs_leap_seconds():
    tb = base()
    assert tb.resolve(SpartnTag(500_000_000, True), System.GLONASS) is FeedReason.NO_LEAP_SECONDS
    assert tb.resolve(SpartnTag(5, False), System.GLONASS, GnssTime(2100, 0.0)) \
        is FeedReason.NO_LEAP_SECONDS
    assert tb.spartn_to_gps(System.GLONASS, 0.0) is FeedReason.NO_LEAP_SECONDS
    assert tb.missing(None, System.GLONASS) is FeedReason.NO_LEAP_SECONDS   # an ephemeris (1020)
    assert tb.missing(None, System.GPS) is None
    # Other constellations do not wait for leap seconds.
    assert isinstance(tb.resolve(SpartnTag(500_000_000, True), System.GALILEO), GnssTime)


def test_rule_T4_glonass_inputs_held_and_released_with_leap_seconds():
    tb = base()
    tb.advance(GnssTime(2100, 1000.0))
    tb.hold("glonass ocb", System.GLONASS, SpartnTag(5, False))
    tb.hold("glonass 1020", System.GLONASS)
    assert tb.release() == []
    assert tb.held == 2
    tb.set_leap_seconds(LeapSeconds(18))
    assert tb.release() == ["glonass ocb", "glonass 1020"]
    assert tb.held == 0


def glonass_tag(year, month, day, hour, minute, second) -> SpartnTag:
    """The 32-bit tag of a GLONASS clock reading (Q2: calendar seconds since 2010-01-01)."""
    days = date(year, month, day).toordinal() - date(2010, 1, 1).toordinal()
    return SpartnTag(int(days * DAY + hour * 3600 + minute * 60 + second), True)


def test_rule_T4_glonass_across_a_leap_second():
    """GLONASS time follows UTC through the leap second at the end of 2016 (GLONASS ICD
    3.3.3), so GPS − GLONASS time steps from 17 s − 3 h to 18 s − 3 h at UTC midnight."""
    tb = base()
    tb.set_leap_seconds(LEAP_2017)
    before = tb.resolve(glonass_tag(2017, 1, 1, 2, 59, 59), System.GLONASS)
    after = tb.resolve(glonass_tag(2017, 1, 1, 3, 0, 0), System.GLONASS)
    assert before == utc(2016, 12, 31, 23, 59, 59, leap=17)
    assert after == utc(2017, 1, 1, leap=18) == LEAP_2017.effective
    assert after - before == 2.0        # the inserted second lies between them


def test_rule_T4_leap_seconds_replaced_by_the_latest_statement():
    tb = base()
    tb.set_leap_seconds(LEAP_2017)
    assert tb.leap_at(LEAP_2017.effective - 1.0) == 17
    assert tb.leap_at(LEAP_2017.effective) == 18
    tb.set_leap_seconds(LeapSeconds(18))                 # after the event, no change announced
    assert tb.leap_at(LEAP_2017.effective - 1.0) == 18
    with pytest.raises(TypeError):
        tb.set_leap_seconds(18)


def test_rule_T4_utc_across_a_leap_second():
    """ZDA readings either side of 23:59:60 and the leap second itself (IS-GPS-200N with
    IRN-IS-200N-003 and -004 20.3.3.5.2.4 b: the second 23:59:60 belongs to the old day)."""
    tb = base()
    tb.set_leap_seconds(LEAP_2017)
    t59 = tb.observe_utc(UtcTime(2016, 12, 31, 23, 59, 59.0))
    t60 = tb.observe_utc(UtcTime(2016, 12, 31, 23, 59, 60.5))
    t00 = tb.observe_utc(UtcTime(2017, 1, 1, 0, 0, 0.0))
    assert t59 == GnssTime(1930, 16.0)
    assert t60 == GnssTime(1930, 17.5)
    assert t00 == GnssTime(1930, 18.0)


def test_rule_T4_negative_leap_second():
    """A negative leap second removes 23:59:59: 23:59:58 is followed by 00:00:00."""
    tb = base()
    effective = utc(2030, 1, 1, leap=17)
    tb.set_leap_seconds(LeapSeconds(18, 17, effective))
    assert tb.observe_utc(UtcTime(2029, 12, 31, 23, 59, 58.0)) == utc(2029, 12, 31, 23, 59, 58, leap=18)
    assert tb.observe_utc(UtcTime(2030, 1, 1, 0, 0, 0.0)) == effective
    assert effective - utc(2029, 12, 31, 23, 59, 58, leap=18) == 1.0


# ---- T5: time from the data alone -------------------------------------------------------


def test_rule_T5_sixteen_bit_tags_held_until_time_is_known():
    tb = base()
    first, second = SpartnTag(100, False), SpartnTag(101, False)
    assert tb.resolve(first) is FeedReason.NO_TIME
    tb.hold("ocb 1", System.GPS, first)
    tb.hold("ocb 2", System.GPS, second)
    # A ZDA without leap seconds gives no time (T4: no guessed count).
    assert tb.observe_utc(UtcTime(2024, 3, 1, 0, 1, 50.0)) is None
    assert tb.current is None and tb.release() == []
    # NAV-TIMELS then ZDA: the time is known, the held messages come out in order.
    tb.set_leap_seconds(LeapSeconds(18))
    t = tb.observe_utc(UtcTime(2024, 3, 1, 0, 1, 50.0))
    assert t == tb.current == utc(2024, 3, 1, 0, 1, 50, leap=18)
    assert tb.release() == ["ocb 1", "ocb 2"]
    assert tb.resolve(first) == utc(2024, 3, 1, 0, 1, 40)


def test_rule_T5_32_bit_tag_gives_the_time():
    tb = base()
    tag32 = SpartnTag(int(seconds(utc(2024, 3, 1, 6)) - SPARTN_EPOCH_DAYS * DAY), True)
    tb.hold("hpac", System.GPS, SpartnTag(21590, False))
    assert tb.observe_tag(SpartnTag(7, False)) is None          # a 16-bit tag is no time
    assert tb.observe_tag(tag32) == utc(2024, 3, 1, 6) == tb.current
    assert tb.release() == ["hpac"]
    assert tb.resolve(SpartnTag(21590, False)) == utc(2024, 3, 1, 5, 59, 50)


def test_rule_T5_glonass_32_bit_tag_without_leap_seconds_gives_no_time():
    tb = base()
    assert tb.observe_tag(glonass_tag(2024, 3, 1, 3, 0, 0), System.GLONASS) is None
    assert tb.current is None


def test_rule_T5_receipt_time_and_never_backwards():
    tb = base()
    assert tb.observe_received(GnssTime(2100, 500.0))
    assert not tb.observe_received(GnssTime(2100, 400.0))
    assert tb.current == GnssTime(2100, 500.0)
    tb.set_leap_seconds(LeapSeconds(18))
    assert tb.observe_utc(UtcTime(1980 + 40, 1, 1, 0, 0, 0.0)) is not None   # an older ZDA
    assert tb.current == GnssTime(2100, 500.0)


def test_rule_T5_receipt_time_resolves_without_current_time():
    tb = base()
    assert tb.missing(SpartnTag(5, False), System.GPS, GnssTime(2100, 0.0)) is None
    tb.hold("with receipt", System.GLONASS, SpartnTag(5, False), GnssTime(2100, 7.0))
    tb.set_leap_seconds(LeapSeconds(18))
    assert tb.release() == ["with receipt"]


def test_rule_T5_gga_time_of_day_placed_in_the_known_day():
    tb = base()
    rover = lambda s: RoverPosition(0.8, 0.1, 400.0, s)   # noqa: E731
    assert tb.observe_rover(rover(10.0)) is None          # no day known yet
    tb.set_leap_seconds(LeapSeconds(18))
    tb.advance(utc(2024, 3, 1, 23, 59, 50, leap=18))
    # Ten seconds after midnight: the next day.
    assert tb.observe_rover(rover(10.0)) == utc(2024, 3, 2, 0, 0, 10, leap=18) == tb.current
    # A late fix from just before midnight: the previous day; the time stays.
    assert tb.observe_rover(rover(86395.0)) == utc(2024, 3, 1, 23, 59, 55, leap=18)
    assert tb.current == utc(2024, 3, 2, 0, 0, 10, leap=18)
    # Across a month end (the date comes from the time scales, not from adding a day).
    tb.advance(utc(2024, 3, 31, 23, 59, 59, leap=18))
    assert tb.observe_rover(rover(1.0)) == utc(2024, 4, 1, 0, 0, 1, leap=18)


def test_rule_T5_hold_limit_drops_the_oldest():
    tb = base(hold_limit=2)
    assert tb.hold("a", System.GPS, SpartnTag(1, False)) is None
    assert tb.hold("b", System.GPS, SpartnTag(2, False)) is None
    assert tb.hold("c", System.GPS, SpartnTag(3, False)) == "a"
    assert tb.dropped == 1 and tb.held == 2
    tb.advance(GnssTime(2100, 0.0))
    assert tb.release() == ["b", "c"]
    with pytest.raises(ValueError):
        base(hold_limit=0)


def test_rule_T5_release_keeps_the_order_of_what_stays():
    tb = base()
    tb.advance(GnssTime(2100, 0.0))
    tb.hold("glo 1", System.GLONASS, SpartnTag(1, False))
    tb.hold("gps", System.GPS, SpartnTag(1, False))
    tb.hold("glo 2", System.GLONASS)
    assert tb.release() == ["gps"]
    tb.set_leap_seconds(LeapSeconds(18))
    assert tb.release() == ["glo 1", "glo 2"]


# ---- the UBX adapter: NAV-TIMELS -------------------------------------------------------


def timels(curr=18, src=2, change=0, src_change=0, to_event=0, wn=0, dn=0, valid=0x01):
    m = ubx.NavTimels(itow=0, version=0, reserved0=(0, 0, 0), src_of_curr_ls=src, curr_ls=curr,
                      src_of_ls_change=src_change, ls_change=change, time_to_ls_event=to_event,
                      date_of_ls_gps_wn=wn, date_of_ls_gps_dn=dn, reserved1=(0, 0, 0), valid=valid)
    frame = next(ubx.scan(ubx.encode(m)))
    return ubx.decode(frame)


def test_nav_timels_current_only():
    assert ubx_adapter.leap_seconds(timels()) == LeapSeconds(18)
    assert ubx_adapter.leap_seconds(timels(change=1, valid=0x01)) == LeapSeconds(18)   # event not valid


def test_nav_timels_untrusted():
    """srcOfCurrLs 0 is the firmware default, "can be outdated"; 255 unknown; valid bit 0 off."""
    assert ubx_adapter.leap_seconds(timels(src=0)) is None
    assert ubx_adapter.leap_seconds(timels(src=255)) is None
    assert ubx_adapter.leap_seconds(timels(valid=0x00)) is None
    assert ubx_adapter.leap_seconds(timels(src=1)) == LeapSeconds(18)   # from GPS − GLONASS


def test_nav_timels_announced_leap_second():
    """The 2016-12-31 leap second: WN 1929, DN 7 (Saturday), 17 → 18."""
    ls = ubx_adapter.leap_seconds(timels(curr=17, change=1, src_change=2, to_event=3600,
                                         wn=1929, dn=7, valid=0x03))
    assert ls == LEAP_2017
    # BeiDou counts days from 0 = Sunday.
    assert ubx_adapter.leap_seconds(timels(curr=17, change=1, src_change=4, wn=1929, dn=6,
                                           valid=0x03)) == LEAP_2017


def test_nav_timels_zero_change_is_no_event():
    """u-blox F9 HPG 1.51 3.15.26.1 (and X20 HPG 2.11 3.14.22.1), lsChange: "If the value is
    0, then the amount of leap seconds did not change and the event should be ignored",
    even with validTimeToLsEvent set and the last event's week and day given."""
    ls = ubx_adapter.leap_seconds(timels(curr=18, change=0, src_change=2, to_event=-10**8,
                                         wn=1929, dn=7, valid=0x03))
    assert ls == LeapSeconds(18)


def test_nav_timels_malformed():
    with pytest.raises(ValueError):
        ubx_adapter.leap_seconds(timels(change=2, wn=1929, dn=7, valid=0x03))
    with pytest.raises(ValueError):
        ubx_adapter.leap_seconds(timels(change=1, src_change=2, wn=1929, dn=0, valid=0x03))


def test_nav_timels_to_time_base():
    tb = base()
    tb.set_leap_seconds(ubx_adapter.leap_seconds(timels(curr=17, change=1, src_change=2,
                                                        wn=1929, dn=7, valid=0x03)))
    assert tb.observe_utc(UtcTime(2017, 1, 1, 0, 0, 0.0)) == LEAP_2017.effective


# ---- the NMEA adapter: ZDA and GGA -------------------------------------------------------


def sentence(text: bytes):
    return nmea.decode(next(nmea.scan(text)))


def test_zda_worked_example():
    """u-blox F9 HPG 1.51 2.7.21.1, X20 HPG 2.11 2.7.20.1: $GPZDA,082710.00,16,09,2002,00,00*64."""
    zda = sentence(b"$GPZDA,082710.00,16,09,2002,00,00*64\r\n")
    assert nmea_adapter.utc_time(zda) == UtcTime(2002, 9, 16, 8, 27, 10.0)


def test_zda_without_time_and_malformed():
    assert nmea_adapter.utc_time(nmea.Zda("GN", "", "", "", "", "00", "00")) is None
    for fields in (("082710.00", "31", "04", "2002"), ("082760.00", "16", "09", "2002"),
                   ("0827.00", "16", "09", "2002"), ("082710.00", "16", "13", "2002"),
                   ("082710.00", "16", "9", "2002")):
        with pytest.raises(ValueError):
            nmea_adapter.utc_time(nmea.Zda("GN", *fields, "00", "00"))
    # The leap second is a valid reading in the last minute of the day only.
    assert nmea_adapter.utc_time(nmea.Zda("GN", "235960.50", "31", "12", "2016", "00", "00")) \
        == UtcTime(2016, 12, 31, 23, 59, 60.5)
    with pytest.raises(ValueError):
        nmea_adapter.utc_time(nmea.Zda("GN", "235860.00", "31", "12", "2016", "00", "00"))


def test_gga_worked_example():
    """u-blox F9 HPG 1.51 2.7.5.1: $GPGGA,092725.00,4717.11399,N,00833.91590,E,1,08,1.01,
    499.6,M,48.0,M,,*5B; height on the ellipsoid is alt + sep."""
    gga = sentence(b"$GPGGA,092725.00,4717.11399,N,00833.91590,E,1,08,1.01,499.6,M,48.0,M,,*5B\r\n")
    p = nmea_adapter.rover_position(gga)
    assert math.degrees(p.latitude) == pytest.approx(47 + 17.11399 / 60, abs=1e-12)
    assert math.degrees(p.longitude) == pytest.approx(8 + 33.91590 / 60, abs=1e-12)
    assert p.height == pytest.approx(547.6, abs=1e-9)
    assert p.utc_seconds == 9 * 3600 + 27 * 60 + 25.0


def test_latitude_longitude_worked_example():
    """u-blox F9 HPG 1.51 2.5.3: 4722.80340 N is 47.38005667°, 00831.68218 E 8.52803633°."""
    assert nmea_adapter.degrees("4722.80340", "N", longitude=False) == pytest.approx(47.38005667, abs=5e-9)
    assert nmea_adapter.degrees("00831.68218", "E", longitude=True) == pytest.approx(8.52803633, abs=5e-9)
    assert nmea_adapter.degrees("4722.80340", "S", longitude=False) < 0
    assert nmea_adapter.degrees("00831.68218", "W", longitude=True) < 0
    for text, hemi, lon in (("4760.0", "N", False), ("9100.0", "N", False), ("18100.0", "E", True),
                            ("4722.8", "E", False), ("472.8", "N", False)):
        with pytest.raises(ValueError):
            nmea_adapter.degrees(text, hemi, lon)


def gga(quality="4", lat="4717.11399", alt="499.6", unit="M"):
    return nmea.Gga("GN", "092725.00", lat, "N", "00833.91590", "E", quality, "12", "0.5",
                    alt, unit, "48.0", "M", "", "")


def test_sentences_with_fields_missing_at_the_end():
    """gnsswire decodes a sentence with fewer fields than its table, the missing
    ones None: a ZDA without its date is no time, a GGA without its position is
    malformed, and a GGA cut after its quality of no fix has nothing to store."""
    assert nmea_adapter.utc_time(nmea.decode(next(nmea.scan(nmea.encode_frame(b"GPZDA,082710.00,16,09"))))) is None
    no_fix = nmea.decode(next(nmea.scan(nmea.encode_frame(b"GPGGA,092725.00,,,,,0"))))
    assert no_fix.alt is None and nmea_adapter.rover_position(no_fix) is None
    with pytest.raises(ValueError):
        nmea_adapter.rover_position(nmea.decode(next(nmea.scan(nmea.encode_frame(b"GPGGA,092725.00,,,,,1")))))


def test_gga_fix_qualities():
    """2.5.5: 1, 2, 4, 5 are fixes; 0 (no fix or limits exceeded) and 6 (dead reckoning,
    also with limits exceeded) are not."""
    for q in "1245":
        assert nmea_adapter.rover_position(gga(q)) is not None
    for q in "0678":
        assert nmea_adapter.rover_position(gga(q)) is None
    assert nmea_adapter.rover_position(gga("0", lat="")) is None   # no fix, null fields
    for bad in (gga(""), gga(lat=""), gga(alt="4x"), gga(unit="F")):
        with pytest.raises(ValueError):
            nmea_adapter.rover_position(bad)
