"""The MSM arithmetic of vrs/msm.py against the tables of RTCM 10403.3 3.5.16."""

from __future__ import annotations

import pytest

from ssrcore.model.types import GnssTime, Sat, Signal, System
from ssrcore.vrs import msm


def test_lock_indicator_follows_table_3_5_74():
    # (lock time ms, DF402): the lower edges of each row and values inside
    for t, i in ((0, 0), (31, 0), (32, 1), (63, 1), (64, 2), (1000, 5), (1024, 6), (262143, 13),
                 (262144, 14), (524287, 14), (524288, 15), (10 ** 9, 15)):
        assert msm.lock_indicator(t) == i, t
        assert msm.lock_minimum(i, 4) <= t


def test_extended_lock_indicator_follows_tables_3_5_75_and_3_5_76():
    # GetIndex of Table 3.5-76 written out row by row, against the closed form
    def get_index(t):
        if t < 64:
            return t
        lo = 64
        for k, base in enumerate(range(64, 704, 32)):
            if lo <= t < 2 * lo:
                return base + (t - lo) // (lo // 32)
            lo *= 2
        return 704

    for t in list(range(0, 300)) + [511, 512, 1000, 1023, 1024, 4095, 65536, 67108863, 67108864, 10 ** 9]:
        i = msm.lock_indicator_extended(t)
        assert i == get_index(t), t
        assert msm.lock_minimum(i, 7) <= t, t
    # Table 3.5-75: minimum lock times of chosen indicators (k·i − c)
    for i, minimum in ((63, 63), (64, 64), (95, 126), (96, 128), (190, 992), (222, 1984), (704, 67108864)):
        assert msm.lock_minimum(i, 7) == minimum


def test_rough_and_fine_ranges_restore_the_range_within_the_resolution():
    for level in (4, 5, 6, 7):
        res_p, res_l = msm.resolution(level)
        for p in (19_123_456.789, 25_987_654.3217, 21_000_000.0):
            rough = msm.rough_range(p)
            r = msm.rough_metres(*rough)
            fp = msm.fine_pseudorange(p, r, level)
            fl = msm.fine_phaserange(p + 12.3456, r, level)
            assert abs(msm.pseudorange_metres(*rough, fp, level) - p) <= res_p / 2 + 1e-8
            assert abs(msm.phaserange_metres(*rough, fl, level) - p - 12.3456) <= res_l / 2 + 1e-8
    assert msm.resolution(7)[1] == pytest.approx(1.396e-4, rel=1e-3)        # DF406, about 0.14 mm
    assert msm.resolution(4)[1] == pytest.approx(5.58e-4, rel=1e-3)         # DF401, about 0.56 mm


def test_fields_out_of_range_and_invalid_patterns():
    rough = msm.rough_metres(*msm.rough_range(2.0e7))
    assert msm.fine_pseudorange(rough + 293.0, rough, 7) is None            # beyond ±2⁻¹⁰ ms (292.8 m)
    assert msm.fine_pseudorange(rough - 292.0, rough, 7) is not None
    assert msm.fine_phaserange(rough + 1172.0, rough, 4) is None            # beyond ±2⁻⁸ ms (1171 m)
    assert msm.fine_phaserange(rough - 1170.0, rough, 4) is not None
    assert msm.rough_range(255 * msm.LIGHT_MS) is None                      # DF397 0 to 254 ms
    assert (msm.invalid_pseudorange(4), msm.invalid_pseudorange(7)) == (-0x4000, -0x80000)
    assert (msm.invalid_phaserange(4), msm.invalid_phaserange(7)) == (-0x200000, -0x800000)


def test_rates():
    assert msm.rates(-523.45678) == (-523, -4568)
    assert msm.rate_metres(*msm.rates(812.00004)) == pytest.approx(812.0, abs=1e-4)
    assert msm.rates(9000.0) == (-8192, -16384)                             # invalid patterns
    assert msm.rates(float("nan")) == (-8192, -16384)


def test_identifiers():
    assert msm.signal_id(Signal.GPS_L1C) == 2 and msm.signal_id(Signal.GPS_L2W) == 10
    assert msm.signal_id(Signal.GPS_L2L) == 16 and msm.signal_id(Signal.GPS_L5Q) == 23
    assert msm.signal_id(Signal.GLONASS_L2C) == 8 and msm.signal_id(Signal.GALILEO_L7Q) == 15
    assert msm.signal_id(Signal.QZSS_L2L) == 16 and msm.signal_id(Signal.BEIDOU_L6I) == 8
    for code in ("1P", "5P", "7P", "8P"):                                   # not in Table 3.5-108
        assert msm.signal_id(Signal(System.BEIDOU, code)) is None
    assert msm.satellite_id(Sat.gps(63)) == 63 and msm.satellite_id(Sat.gps(64)) is None
    assert msm.satellite_id(Sat.glonass(24)) == 24 and msm.satellite_id(Sat.glonass(25)) is None
    assert msm.satellite_id(Sat.galileo(51)) is None and msm.satellite_id(Sat.beidou(38)) is None
    assert msm.satellite_id(Sat.qzss(1)) == 1 and msm.satellite_id(Sat.qzss(11)) is None
    assert msm.frequency_channel(-7) == 0 and msm.frequency_channel(6) == 13


def test_epoch_times_of_each_system():
    t = GnssTime(2295, 87400.0)                       # Monday 00:16:40 GPS time
    ms = 87_400_000
    assert msm.epoch_fields(System.GPS, t, None)["epoch"] == ms
    assert msm.epoch_fields(System.QZSS, t, None)["epoch"] == ms
    assert msm.epoch_fields(System.GALILEO, t, None)["epoch"] == ms        # GST weeks start with GPS's
    assert msm.epoch_fields(System.BEIDOU, t, None)["epoch"] == ms - 14_000
    assert msm.epoch_fields(System.GLONASS, t, None) is None
    g = msm.epoch_fields(System.GLONASS, t, 18)
    # GPS − 18 s + 3 h: Monday 03:16:22 Moscow time
    assert (g["epoch"], g["glonass_day"], g["glonass_epoch"]) == (None, 1, (3 * 3600 + 16 * 60 + 22) * 1000)
    # near the end of the week, BeiDou and GLONASS wrap into the next day or week
    late = GnssTime(2295, 604_790.0)
    assert msm.epoch_fields(System.BEIDOU, GnssTime(2295, 5.0), None)["epoch"] == 604_791_000
    g = msm.epoch_fields(System.GLONASS, late, 18)
    assert g["glonass_day"] == 0 and g["glonass_epoch"] == (3 * 3600 - 18 - 10) * 1000
