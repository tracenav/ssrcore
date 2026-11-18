"""The ``timescale`` kernel through its toolkit module (DESIGN 2.7).

Independent routes used here, not the kernel's own arithmetic:
- Python's ``datetime`` (proleptic Gregorian calendar) for day counts, dates
  and days of year over the whole valid range;
- published epochs: JD 2 444 244.5 for 1980-01-06 and J2000.0 = JD 2 451 545.0
  (TM-23 3.1.1); the GPS week rollovers of 1999-08-22 and 2019-04-07 (IS-GPS-200N with
  IRN-IS-200N-003 and -004 6.2.4 numbering); the GST and BDT start epochs stated in UTC (Galileo OS
  SIS ICD 2.2 5.1.2, BDS ICD 3.3) converted with GPS - UTC as published by the IERS in
  Bulletin C (13 s in 1999, 14 s in 2006, 15 s in 2010, 17 s to 18 s at the end
  of 2016); the worked example of GLONASS ICD A.3.1.3 (N4 = 2, NT = 839: 2002);
- the leap second at the end of 2016-12-31 (IERS Bulletin C 52), checked
  against the labels UTC must take around it, 23:59:59, 23:59:60, 00:00:00.
"""

import math
import random
from datetime import date, timedelta

import pytest

from ssrcore.toolkit import timescale as ts
from ssrcore.toolkit.timescale import System as Scale
from ssrcore.toolkit.timescale import UtcParameters

EPOCH = date(1980, 1, 6)
WEEK = 604800.0
DAY = 86400.0


def days(d: date) -> int:
    """Independent day count from 1980-01-06, by datetime."""
    return (d - EPOCH).days


def gps(d: date, seconds: float = 0.0) -> tuple[int, float]:
    """GPS (week, sow) of a GPS calendar reading, by datetime."""
    n = days(d)
    return n // 7, (n % 7) * DAY + seconds


# ---- calendar -------------------------------------------------------------------

def test_calendar_matches_datetime_over_the_whole_valid_range():
    first, last = date(1900, 3, 1), date(2100, 2, 28)
    d = first
    while d <= last:
        n = days(d)
        assert ts.days_from_date(d.year, d.month, d.day) == n
        assert ts.date_from_days(n) == (d.year, d.month, d.day)
        if 1901 <= d.year <= 2099:
            assert ts.day_of_year(n) == d.timetuple().tm_yday
        d += timedelta(days=1)


def test_calendar_refuses_dates_outside_the_range_and_dates_that_do_not_exist():
    for y, m, d in [(1900, 2, 28), (2100, 3, 1), (2001, 2, 29), (2024, 4, 31),
                    (2024, 13, 1), (2024, 0, 1), (2024, 1, 0), (2024, 1, 32)]:
        with pytest.raises(ValueError):
            ts.days_from_date(y, m, d)
    for n in (days(date(1900, 3, 1)) - 1, days(date(2100, 2, 28)) + 1):
        with pytest.raises(ValueError):
            ts.date_from_days(n)
    assert ts.days_from_date(2000, 2, 29) == days(date(2000, 2, 29))  # 2000 is a leap year


def test_julian_dates_of_tm23():
    # TM-23 3.1.1: 6 January 1980 0h = JD 2 444 244.5; J2000.0 = JD 2 451 545.0 (12h).
    assert ts.days_from_date(1980, 1, 6) == 0
    assert ts.days_from_date(2000, 1, 1) == int(2451545.0 - 0.5 - 2444244.5)


# ---- weeks ------------------------------------------------------------------------

def test_gps_week_rollovers():
    # Week 1024 began 1999-08-22 and week 2048 2019-04-07 (the two published rollovers).
    assert ts.days_from_date(1999, 8, 22) == 1024 * 7
    assert ts.days_from_date(2019, 4, 7) == 2048 * 7
    # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1: ten bits of the week, resolved near a
    # known week.
    assert ts.resolve_week(0, 10, 2047) == 2048
    assert ts.resolve_week(1023, 10, 2048) == 2047
    assert ts.resolve_week(2048 % 1024, 10, 2048 + 511) == 2048
    assert ts.resolve_week(2048 % 1024, 10, 2048 - 511) == 2048
    assert ts.resolve_week(2048 % 1024, 10, 2048 - 512) == 1024  # [ref - 512, ref + 512)
    # Galileo OS SIS ICD 2.2 5.1.2 (12 bits, 4096 weeks) and BDS ICD 5.2.4.4 (13 bits, 8192).
    assert ts.resolve_week(5, 12, 4095) == 4101
    assert ts.resolve_week(8191, 13, 8193) == 8191
    # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4: WN_LSF is 8 bits, within 127 of WN.
    for wn in range(1900, 2300):
        for k in range(-127, 128):
            assert ts.resolve_week((wn + k) % 256, 8, wn) == wn + k
    for bad in [(256, 8, 0), (-1, 8, 0), (0, 0, 0), (0, 31, 0)]:
        with pytest.raises(ValueError):
            ts.resolve_week(*bad)


def test_end_of_week_crossover():
    # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1: times on either side of the end/start
    # of week.
    assert ts.add(1929, 604799.5, 1.0) == (1930, 0.5)
    assert ts.add(1930, 0.5, -1.0) == (1929, 604799.5)
    assert ts.normalize(1930, -3 * WEEK - 1.0) == (1926, WEEK - 1.0)
    assert ts.normalize(1930, 2 * WEEK) == (1932, 0.0)
    assert ts.diff(1930, 0.5, 1929, 604799.5) == 1.0
    # Weeks are differenced exactly: the result keeps the resolution of the seconds
    # of week (2**-33 s, the spacing of doubles just below 604 800) at any week.
    tiny = 2.0**-33
    assert ts.diff(4000, 302400.0 + tiny, 3999, 302400.0) == WEEK + tiny
    assert ts.diff(4000, tiny, 4000, 0.0) == tiny


# ---- scales -----------------------------------------------------------------------

def test_gst_epoch_is_gps_week_1024():
    # Galileo OS SIS ICD 2.2 5.1.2: GST = 13 s at 1999-08-22 00:00:00 UTC. With GPS - UTC = 13 s
    # then (IERS Bulletin C), GST 0 is GPS week 1024 starting.
    gst_zero_in_gps = ts.add(*ts.utc_to_gps(days(date(1999, 8, 22)), 0.0, 13), -13.0)
    assert ts.gst_to_gps(0, 0.0) == gst_zero_in_gps == (1024, 0.0)
    assert ts.gps_to_gst(2300, 1234.5) == (2300 - 1024, 1234.5)
    assert ts.gst_to_gps(*ts.gps_to_gst(2300, 1234.5)) == (2300, 1234.5)


def test_bdt_epoch_and_offset():
    # BDS ICD 3.3: BDT starts at 2006-01-01 00:00:00 UTC; GPS - UTC was 14 s then.
    assert ts.bdt_to_gps(0, 0.0) == ts.utc_to_gps(days(date(2006, 1, 1)), 0.0, 14)
    # RTCM DF427: BDT time of week is 14 s less than GPS time of week.
    assert ts.gps_to_bdt(2300, 100.0) == (2300 - 1356, 86.0)
    # ... across the start of a BDT week.
    assert ts.gps_to_bdt(2300, 10.0) == (2300 - 1356 - 1, WEEK - 4.0)
    assert ts.bdt_to_gps(*ts.gps_to_bdt(2300, 10.0)) == (2300, 10.0)
    # BDS ICD 5.2.4.19: A0GPS + A1GPS tE.
    assert ts.bdt_offset(2e-9, 1e-12, 1000.0) == pytest.approx(2e-9 + 1e-9, rel=1e-15)


def test_ggto_equation_24():
    # Galileo OS SIS ICD 2.2 5.1.8 Eq. 24: A0G + A1G (TOW - t0G + 604800 (WN - WN0G)), WN0G of 6 bits.
    wn, tow, a0g, a1g, wn0g, t0g = 1000, 3600.0, 5e-9, 2e-13, 999, 7200.0
    expected = a0g + a1g * (tow - t0g + WEEK * (wn - wn0g))
    assert ts.ggto(wn, tow, a0g, a1g, wn0g % 64, t0g) == pytest.approx(expected, rel=1e-15)


def test_glonass_time_is_utc_plus_three_hours():
    # GLONASS ICD 3.3.3: T_GLONASS = T_UTC(SU) + 3 h.
    d = days(date(2024, 5, 1))
    for sod in (0.0, 10799.0, 10800.0, 50000.25):
        assert ts.glonass_to_gps(d, sod, 18) == ts.utc_to_gps(d, sod - 10800.0, 18)
    # GLONASS midnight is 21:00 UTC of the previous day.
    assert ts.glonass_to_gps(d, 0.0, 18) == gps(date(2024, 4, 30), 21 * 3600.0 + 18.0)
    assert ts.gps_to_glonass(*ts.glonass_to_gps(d, 1234.5, 18), 18) == (d, 1234.5)


def test_glonass_n4_nt_worked_example_and_calendar():
    # GLONASS ICD A.3.1.3: N4 = 2, NT = 839 gives J = 3 and the year 2002.
    assert ts.date_from_days(ts.glonass_day(2, 839))[0] == 2002
    # Day and month from the civil calendar: NT counts from 1 January of the leap year.
    for n4, nt in [(1, 1), (2, 839), (7, 366), (7, 367), (8, 1461), (8, 1097)]:
        first = date(1996 + 4 * (n4 - 1), 1, 1)
        assert ts.glonass_day(n4, nt) == days(first + timedelta(days=nt - 1))
        assert ts.glonass_n4_nt(ts.glonass_day(n4, nt)) == (n4, nt)
    for bad in [(0, 1), (32, 1), (1, 0), (1, 1462)]:
        with pytest.raises(ValueError):
            ts.glonass_day(*bad)


# ---- UTC and leap seconds ----------------------------------------------------------

# IERS Bulletin C 52: a positive leap second at the end of 2016-12-31; GPS - UTC
# went from 17 s to 18 s. 2016-12-31 is a Saturday of GPS week 1929: DN = 7.
LEAP_2016 = UtcParameters(dt_ls=17, dt_lsf=18, wn_lsf=1929 % 256, dn=7,
                          wn_t=1929 % 256, t_ot=0.0)
DEC31, JAN1 = days(date(2016, 12, 31)), days(date(2017, 1, 1))


def test_week_1929_holds_the_2016_leap_second():
    assert days(date(2016, 12, 25)) == 1929 * 7
    assert date(2016, 12, 31).isoweekday() == 6  # Saturday, day 7 of the GPS week


def test_gps_to_utc_through_the_2016_leap_second():
    # Around the event: UTC 23:59:59, 23:59:60, then 00:00:00 of 1 January.
    assert ts.gps_to_utc(1930, 16.5, LEAP_2016) == (DEC31, 86399.5, 2)
    assert ts.gps_to_utc(1930, 17.5, LEAP_2016) == (DEC31, 86400.5, 2)
    assert ts.gps_to_utc(1930, 18.5, LEAP_2016) == (JAN1, 0.5, 2)
    # Outside the six-hour window: case a with 17 s, case c with 18 s.
    assert ts.gps_to_utc(1929, 0.0, LEAP_2016) == (1929 * 7 - 1, DAY - 17.0, 1)
    assert ts.gps_to_utc(1931, 0.0, LEAP_2016) == (1931 * 7 - 1, DAY - 18.0, 3)
    # UTC is continuous across both window edges and rises by 1 s per second, except
    # for the leap second itself, which repeats nothing and skips nothing.
    previous = None
    for k in range(-7 * 3600, 7 * 3600):
        t = ts.add(1930, 0.0, float(k))
        d, s, _ = ts.gps_to_utc(*t, LEAP_2016)
        if previous is not None:
            pd, ps = previous
            if d == pd:
                assert s - ps == 1.0
            else:
                assert d == pd + 1 and s == 0.0 and ps in (86399.0, 86400.0)
        previous = (d, s)
    # Round trip through utc_to_gps, with the count in force at that UTC.
    assert ts.utc_to_gps(DEC31, 86400.5, 17) == (1930, 17.5)
    assert ts.utc_to_gps(JAN1, 0.5, 18) == (1930, 18.5)


def test_gps_to_utc_negative_leap_second_skips_23_59_59():
    p = UtcParameters(dt_ls=18, dt_lsf=17, wn_lsf=2400 % 256, dn=3)
    end = ts.add(2400, 3 * DAY, 18.0)  # 00:00:00 UTC that would have been
    d, s, case = ts.gps_to_utc(*ts.add(*end, -1.5), p)
    assert (s, case) == (86398.5, 2)
    d2, s2, _ = ts.gps_to_utc(*ts.add(*end, -0.5), p)
    assert (d2, s2) == (d + 1, 0.5)


def test_gps_to_utc_polynomial():
    # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a: A0 + A1 (tE - tot + 604800 (WN -
    # WNt)), WNt of 8 bits.
    p = UtcParameters(dt_ls=18, dt_lsf=18, wn_lsf=0, dn=1, a0=3e-9, a1=-4e-14,
                      wn_t=2399 % 256, t_ot=405504.0)
    week, sow = 2400, 1000.0
    poly = 3e-9 + -4e-14 * (sow - 405504.0 + WEEK)
    d, s, case = ts.gps_to_utc(week, sow, p)
    assert case == 1 and d == 2400 * 7
    assert s == pytest.approx(sow - 18.0 - poly, abs=1e-12)


def test_gst_to_utc_same_event():
    # Galileo OS SIS ICD 2.2 5.1.7, GST weeks; the 2016 leap second as Galileo broadcast it.
    p = UtcParameters(dt_ls=17, dt_lsf=18, wn_lsf=(1929 - 1024) % 256, dn=7,
                      wn_t=(1929 - 1024) % 256)
    assert ts.gst_to_utc(1930 - 1024, 17.5, p) == (DEC31, 86400.5, 2)
    # Case a with delta t_LS before the window: Eq. 21 of issue 2.2, which corrects
    # the delta t_LSF of issue 2.1's Eq. 20; case c with delta t_LSF after it (Eq. 23).
    assert ts.gst_to_utc(1929 - 1024, 0.0, p) == (1929 * 7 - 1, DAY - 17.0, 1)
    assert ts.gst_to_utc(1931 - 1024, 0.0, p) == (1931 * 7 - 1, DAY - 18.0, 3)


def test_bdt_to_utc_same_event():
    # BDS ICD 5.2.4.18: BDT - UTC went from 3 s to 4 s; DN 6 is the Saturday.
    p = UtcParameters(dt_ls=3, dt_lsf=4, wn_lsf=(1929 - 1356) % 256, dn=6)
    assert ts.bdt_to_utc(1930 - 1356, 3.5, p) == (DEC31, 86400.5, 2)
    assert ts.bdt_to_utc(1930 - 1356, 4.5, p) == (JAN1, 0.5, 2)
    # Window DN + 2/3 to DN + 5/4: 8 h before to 6 h after the end of day DN.
    assert ts.bdt_to_utc(1930 - 1356, -8 * 3600.0 - 1.0, p)[2] == 1
    assert ts.bdt_to_utc(1930 - 1356, -8 * 3600.0 + 1.0, p)[2] == 2
    assert ts.bdt_to_utc(1930 - 1356, 6 * 3600.0 - 1.0, p)[2] == 2
    assert ts.bdt_to_utc(1930 - 1356, 6 * 3600.0 + 1.0, p)[2] == 3
    # The same instant through GPS time agrees.
    assert ts.gps_to_utc(*ts.bdt_to_gps(1930 - 1356, 3.5), LEAP_2016)[:2] == (DEC31, 86400.5)
    with pytest.raises(ValueError):
        ts.bdt_to_utc(500, 0.0, UtcParameters(3, 4, 0, 7))


# ---- SPARTN time tags -----------------------------------------------------------------

def test_spartn_32_bit_tag_epoch():
    # SPARTN TF009: seconds since 2010-01-01 00:00:00 in the message's system time.
    assert ts.spartn_tag_to_gps(0, Scale.GPS) == gps(date(2010, 1, 1)) == (1564, 432000.0)
    assert ts.spartn_tag_to_gps(0, Scale.GALILEO) == (1564, 432000.0)
    # BDT reads 00:00:00 when GPS time reads 00:00:14.
    assert ts.spartn_tag_to_gps(0, Scale.BEIDOU) == gps(date(2010, 1, 1), 14.0)
    # GLONASS midnight is 21:00 UTC; GPS - UTC was 15 s in 2010 (IERS Bulletin C).
    assert ts.spartn_tag_to_gps(0, Scale.GLONASS, 15) == gps(date(2009, 12, 31), 21 * 3600.0 + 15.0)
    with pytest.raises(ValueError):
        ts.spartn_tag_to_gps(0, Scale.GLONASS)
    # One day of tag is one day of GPS time; the largest tag is in 2146.
    assert ts.spartn_tag_to_gps(86400, Scale.GPS) == gps(date(2010, 1, 2))
    w, s = ts.spartn_tag_to_gps(2**32 - 1, Scale.GPS)
    assert ts.diff(w, s, 1564, 432000.0) == 2**32 - 1
    assert (EPOCH + timedelta(days=w * 7 + s // DAY)).year == 2146


def test_spartn_tag_round_trip_and_range():
    rng = random.Random(2010)
    for scale, leap in [(Scale.GPS, None), (Scale.GALILEO, None), (Scale.BEIDOU, None),
                        (Scale.GLONASS, 18)]:
        for _ in range(2000):
            tag = rng.randrange(0, 2**32)
            w, s = ts.spartn_tag_to_gps(tag, scale, leap)
            assert ts.spartn_tag_from_gps(w, s, scale, leap) == tag
            assert ts.spartn_tag_from_gps(w, s + 0.999, scale, leap) == tag
            # Just below the next second: still this tag, never the next (the count from
            # 2010 as one double cannot hold the fraction).
            assert ts.spartn_tag_from_gps(w, math.nextafter(s + 1.0, 0.0), scale, leap) == tag
            assert ts.spartn_tag_from_gps(w, s + 0.99999998, scale, leap) == tag
    with pytest.raises(ValueError):
        ts.spartn_tag_from_gps(1564, 431999.0, Scale.GPS)
    # QZSS time keeps the GPS week and has no whole-second offset (RTCM DF428, unlike
    # DF427 for BeiDou); a code that names no system is refused.
    assert ts.spartn_tag_to_gps(123456, Scale.QZSS) == ts.spartn_tag_to_gps(123456, Scale.GPS)
    with pytest.raises(ValueError):
        ts.spartn_tag_to_gps(0, 5)


def test_spartn_half_day_resolution_rule_T2():
    # A 16-bit tag is the 32-bit tag modulo 43 200; any reference within 6 h of the
    # time recovers the 32-bit tag (DESIGN T2).
    rng = random.Random(43200)
    for _ in range(20000):
        tag = rng.randrange(43200, 2**32 - 43200)
        ref = tag + rng.randrange(-21599, 21600)
        assert ts.spartn_resolve_half_day(tag % 43200, ref) == tag
    # Across a half-day edge in both directions.
    edge = 43200 * 1000
    assert ts.spartn_resolve_half_day(43199, edge + 5) == edge - 1
    assert ts.spartn_resolve_half_day(3, edge - 5) == edge + 3
    # At an exact tie (6 h either way) the earlier half day.
    assert ts.spartn_resolve_half_day(0, edge + 21600) == edge
    # Near the start of the count there is no earlier half day.
    assert ts.spartn_resolve_half_day(43000, 100) == 43000
    with pytest.raises(ValueError):
        ts.spartn_resolve_half_day(43200, edge)


def test_spartn_received_time_to_tag_then_resolve():
    # The receipt time as a 32-bit reference in the tag's own scale (BeiDou here).
    received = gps(date(2024, 5, 1), 12 * 3600.0 + 30.0)
    tag32 = ts.spartn_tag_from_gps(*ts.add(*received, -2.0), Scale.BEIDOU)
    ref = ts.spartn_tag_from_gps(*received, Scale.BEIDOU)
    assert ts.spartn_resolve_half_day(tag32 % 43200, ref) == tag32
    assert ts.spartn_tag_to_gps(tag32, Scale.BEIDOU) == ts.add(*received, -2.0)
