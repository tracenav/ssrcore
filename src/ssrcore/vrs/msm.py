"""The arithmetic of a Multiple Signal Message: identifiers, rough and fine
ranges, lock time indicators, rates and epoch times.

RTCM 10403.3 section 3.5.16. Ranges are carried in light-milliseconds, with
c = 299 792 458 m/s (3.5.16.1), as a satellite's rough range (DF397 whole
milliseconds and DF398 in 2⁻¹⁰ ms) and per signal a fine pseudorange and
fine phase range against it, at the standard resolution of MSM1 to MSM5 or
the extended one of MSM6 and MSM7 (3.5.16.3.1). Every function takes and
returns numbers; the observations module assembles them into
``gnsswire.rtcm3.Msm`` dataclasses, whose bit layout is gnsswire's.

[design: DESIGN 2.6.3, the encoding of a wire format written here rather
than as a kernel]
"""

from __future__ import annotations

import math

from ssrcore.model.types import SPEED_OF_LIGHT, WEEK_SECONDS, GnssTime, Sat, Signal, System
from ssrcore.toolkit import timescale

# [spec: RTCM 10403.3 3.5.16.1, "In MSM it is assumed that the speed of light is
# c = 299,792,458 meters per second"; the same value as SPARTN Table 8.0] metres per
# light-millisecond.
LIGHT_MS = SPEED_OF_LIGHT / 1000.0

# [spec: RTCM 10403.3 Table 3.5-71] the first three digits of each system's MSM numbers.
PREFIX = {System.GPS: 107, System.GLONASS: 108, System.GALILEO: 109, System.QZSS: 111, System.BEIDOU: 112}

# [spec: RTCM 10403.3 Tables 3.5-90, 3.5-95, 3.5-98, 3.5-104, 3.5-107] the highest
# satellite ID with a satellite: GPS PRN 1 to 63 (64 reserved), GLONASS slots 1
# to 24, Galileo PRN 1 to 50 (51 and 52 are GIOVE-A and -B, not PRNs), QZSS 1 to
# 10 (PRN 193 to 202, which are SPARTN's QZSS numbers 1 to 10, ssrcore.model.types.Sat),
# BeiDou PRN 1 to 37. The ID is the satellite's number in every system.
LAST_SATELLITE = {System.GPS: 63, System.GLONASS: 24, System.GALILEO: 50, System.QZSS: 10, System.BEIDOU: 37}

# [spec: RTCM 10403.3 Tables 3.5-91 (GPS), 3.5-96 (GLONASS), 3.5-99 (Galileo),
# 3.5-105 (QZSS), 3.5-108 (BeiDou)] signal ID in DF395 of each RINEX code.
SIGNAL_ID = {
    System.GPS: {"1C": 2, "1P": 3, "1W": 4, "2C": 8, "2P": 9, "2W": 10, "2S": 15, "2L": 16, "2X": 17,
                 "5I": 22, "5Q": 23, "5X": 24, "1S": 30, "1L": 31, "1X": 32},
    System.GLONASS: {"1C": 2, "1P": 3, "2C": 8, "2P": 9},
    System.GALILEO: {"1C": 2, "1A": 3, "1B": 4, "1X": 5, "1Z": 6, "6C": 8, "6A": 9, "6B": 10, "6X": 11,
                     "6Z": 12, "7I": 14, "7Q": 15, "7X": 16, "8I": 18, "8Q": 19, "8X": 20, "5I": 22,
                     "5Q": 23, "5X": 24},
    System.QZSS: {"1C": 2, "6S": 9, "6L": 10, "6X": 11, "2S": 15, "2L": 16, "2X": 17, "5I": 22, "5Q": 23,
                  "5X": 24, "1S": 30, "1L": 31, "1X": 32},
    System.BEIDOU: {"2I": 2, "2Q": 3, "2X": 4, "6I": 8, "6Q": 9, "6X": 10, "7I": 14, "7Q": 15, "7X": 16},
}

# [spec: RTCM 10403.3 Table 3.4-1 DF397, 0 to 254 ms; 255 is the invalid value]
_ROUGH_MS_LIMIT = 254
# [spec: RTCM 10403.3 DF398, 2⁻¹⁰ ms] rough range units per millisecond.
_ROUGH_UNITS = 1024
# [spec: RTCM 10403.3 Table 3.5-72 and DF400, DF401, DF405, DF406] per resolution class
# (standard MSM1 to MSM5, extended MSM6 and MSM7): the binary exponent of the fine
# pseudorange and fine phase range resolutions in ms, and their widths in bits.
_FINE = {False: (24, 15, 29, 22), True: (29, 20, 31, 24)}
# [spec: RTCM 10403.3 DF399, ±8191 m/s, 1 m/s; DF404, ±1.6383 m/s, 0.0001 m/s]
_ROUGH_RATE_LIMIT = 8191
_FINE_RATE_UNIT = 0.0001
_FINE_RATE_LIMIT = 16383
# [spec: RTCM 10403.3 Table 3.4-6 DF419, the frequency channel number plus 7]
_CHANNEL_OFFSET = 7
# [spec: RTCM 10403.3 DF004, DF248, DF427, DF428, ms of the week; DF034, ms of the day]
_MS_PER_WEEK = int(WEEK_SECONDS) * 1000
_MS_PER_DAY = 86400 * 1000


def extended(level: int) -> bool:
    """[spec: RTCM 10403.3 Table 3.5-72] whether an MSM level carries the extended
    resolution fields (MSM6, MSM7)."""
    return level >= 6


def has_rates(level: int) -> bool:
    """[spec: RTCM 10403.3 Table 3.5-72] whether an MSM level carries phase range
    rates and the extended satellite information (MSM5, MSM7)."""
    return level in (5, 7)


def round_half_away(value: float) -> int:
    """[design: a value quantised to a field's resolution is rounded to the nearest
    unit, halves away from zero, so both signs are treated alike]"""
    n = math.floor(abs(value) + 0.5)
    return int(n) if value >= 0.0 else -int(n)


def satellite_id(sat: Sat) -> int | None:
    """[spec: RTCM 10403.3 Tables 3.5-90, 3.5-95, 3.5-98, 3.5-104, 3.5-107] the
    satellite's ID in DF394, or None when RTCM defines none for it."""
    last = LAST_SATELLITE.get(sat.system)
    return sat.number if last is not None and sat.number <= last else None


def signal_id(signal: Signal) -> int | None:
    """[spec: RTCM 10403.3 Tables 3.5-91, 3.5-96, 3.5-99, 3.5-105, 3.5-108] the signal's
    ID in DF395, or None when RTCM defines none for its RINEX code."""
    return SIGNAL_ID.get(signal.system, {}).get(signal.code)


def rough_range(pseudorange: float) -> tuple[int, int] | None:
    """(DF397, DF398) of a range in metres: rounded to 2⁻¹⁰ ms and split into whole
    milliseconds and the remainder; None beyond the 254 ms of DF397.

    [spec: RTCM 10403.3 DF397, DF398]
    """
    units = round_half_away(pseudorange / LIGHT_MS * _ROUGH_UNITS)
    whole, rest = divmod(units, _ROUGH_UNITS)
    if not 0 <= whole <= _ROUGH_MS_LIMIT:
        return None
    return whole, rest


def rough_metres(rough_ms: int, rough_mod: int) -> float:
    """[spec: RTCM 10403.3 3.5.16.3.1] the rough range in metres, c/1000 × (Nms + Rough_range/1024)."""
    return LIGHT_MS * (rough_ms + rough_mod / _ROUGH_UNITS)


def _fine(value: float, rough: float, exponent: int, bits: int) -> int | None:
    """``value − rough`` (m) in units of 2^−exponent ms, or None when it does not fit
    ``bits`` signed bits without the most negative pattern, which is the invalid value."""
    n = round_half_away((value - rough) / LIGHT_MS * 2.0 ** exponent)
    limit = (1 << (bits - 1)) - 1
    return n if -limit <= n <= limit else None


def fine_pseudorange(pseudorange: float, rough: float, level: int) -> int | None:
    """DF400 (MSM1 to MSM5) or DF405 (MSM6, MSM7) of a pseudorange against the rough
    range, both in metres; None when it does not fit.

    [spec: RTCM 10403.3 DF400, 2⁻²⁴ ms, int15; DF405, 2⁻²⁹ ms, int20; 3.5.16.3.1]
    """
    exponent, bits, _, _ = _FINE[extended(level)]
    return _fine(pseudorange, rough, exponent, bits)


def fine_phaserange(phaserange: float, rough: float, level: int) -> int | None:
    """DF401 (MSM1 to MSM5) or DF406 (MSM6, MSM7) of a phase range against the rough
    range, both in metres; None when it does not fit.

    [spec: RTCM 10403.3 DF401, 2⁻²⁹ ms, int22; DF406, 2⁻³¹ ms, int24; 3.5.16.3.1]
    """
    _, _, exponent, bits = _FINE[extended(level)]
    return _fine(phaserange, rough, exponent, bits)


def invalid_pseudorange(level: int) -> int:
    """[spec: RTCM 10403.3 DF400, 4000h (−2⁻¹⁰ ms); DF405, 80000h] the invalid fine pseudorange."""
    return -(1 << (_FINE[extended(level)][1] - 1))


def invalid_phaserange(level: int) -> int:
    """[spec: RTCM 10403.3 DF401, 200000h (−2⁻⁸ ms); DF406, 800000h] the invalid fine phase range."""
    return -(1 << (_FINE[extended(level)][3] - 1))


def pseudorange_metres(rough_ms: int, rough_mod: int, fine: int, level: int) -> float:
    """[spec: RTCM 10403.3 3.5.16.3.1] Pseudorange(i) = c/1000 × (Nms + Rough_range/1024
    + 2⁻²⁴ (2⁻²⁹ for MSM6, MSM7) × Fine_Pseudorange(i)), metres."""
    exponent = _FINE[extended(level)][0]
    return LIGHT_MS * (rough_ms + rough_mod / _ROUGH_UNITS + fine * 2.0 ** -exponent)


def phaserange_metres(rough_ms: int, rough_mod: int, fine: int, level: int) -> float:
    """[spec: RTCM 10403.3 3.5.16.3.1] Phaserange(i) = c/1000 × (Nms + Rough_range/1024
    + 2⁻²⁹ (2⁻³¹ for MSM6, MSM7) × Fine_Phaserange(i)), metres."""
    exponent = _FINE[extended(level)][2]
    return LIGHT_MS * (rough_ms + rough_mod / _ROUGH_UNITS + fine * 2.0 ** -exponent)


def resolution(level: int) -> tuple[float, float]:
    """[spec: RTCM 10403.3 DF400, DF401, DF405, DF406] the resolutions of the fine
    pseudorange and phase range in metres."""
    exponent_p, _, exponent_l, _ = _FINE[extended(level)]
    return LIGHT_MS * 2.0 ** -exponent_p, LIGHT_MS * 2.0 ** -exponent_l


def rates(rate: float) -> tuple[int, int]:
    """(DF399, DF404) of a phase range rate in m/s: the whole metres per second
    nearest and the rest in 0.0001 m/s; the invalid patterns (−8192, −16384) beyond
    their ranges.

    [spec: RTCM 10403.3 DF399, ±8191 m/s, 2000h invalid; DF404, ±1.6383 m/s, 4000h
    invalid; full rate = DF399 + 0.0001 × DF404 (3.5.16.3.1)]
    """
    if not math.isfinite(rate):
        return -_ROUGH_RATE_LIMIT - 1, -_FINE_RATE_LIMIT - 1
    rough = round_half_away(rate)
    if abs(rough) > _ROUGH_RATE_LIMIT:
        return -_ROUGH_RATE_LIMIT - 1, -_FINE_RATE_LIMIT - 1
    return rough, round_half_away((rate - rough) / _FINE_RATE_UNIT)


def rate_metres(rough: int, fine: int) -> float:
    """[spec: RTCM 10403.3 3.5.16.3.1] PhaseRangeRate(i) = Rough_PhaseRangeRate + 0.0001 ×
    Fine_PhaseRangeRate(i), m/s."""
    return rough + _FINE_RATE_UNIT * fine


def lock_indicator(milliseconds: int) -> int:
    """DF402 of a lock time in whole milliseconds: 0 below 32 ms, i for 2^(i+4) ≤ t <
    2^(i+5), 15 from 524 288 ms.

    [spec: RTCM 10403.3 Table 3.5-74]
    """
    if milliseconds < 32:
        return 0
    return min(15, milliseconds.bit_length() - 5)


def lock_indicator_extended(milliseconds: int) -> int:
    """DF407 of a lock time in whole milliseconds.

    [spec: RTCM 10403.3 Table 3.5-76 and its GetIndex algorithm: t below 64 ms is
    itself; in [2^n, 2^(n+1)) for n from 6 to 25 it is 64 + 32 (n − 6) + (t − 2^n) /
    2^(n−5); from 67 108 864 ms it is 704]
    """
    t = milliseconds
    if t < 64:
        return t
    if t >= 67108864:
        return 704
    n = t.bit_length() - 1
    return 64 + 32 * (n - 6) + ((t - (1 << n)) >> (n - 5))


def lock_minimum(indicator: int, level: int) -> int:
    """The minimum lock time in ms an indicator stands for, as a decoder reconstructs it.

    [spec: RTCM 10403.3 Table 3.5-74 (DF402: 0, then 2^(i+4)); Table 3.5-75 (DF407:
    k·i − c for the row of i)]
    """
    if not extended(level):
        return 0 if indicator == 0 else 1 << (indicator + 4)
    if indicator < 64:
        return indicator
    row = (indicator - 64) // 32          # 64–95 is row 0, k = 2
    k = 2 << row
    return k * indicator - k * (64 + 32 * row) + (64 << row)


def lock(milliseconds: int, level: int) -> int:
    """[spec: RTCM 10403.3 DF402 (MSM1 to MSM5), DF407 (MSM6, MSM7)] the lock time
    indicator of the level."""
    return lock_indicator_extended(milliseconds) if extended(level) else lock_indicator(milliseconds)


def frequency_channel(channel: int) -> int:
    """[spec: RTCM 10403.3 3.5.16.5.2, Table 3.4-6] DF419 of a GLONASS frequency channel K: K + 7."""
    return channel + _CHANNEL_OFFSET


def _milliseconds(week: int, seconds: float) -> int:
    """[design: an epoch lies on the grid within 1 µs (DESIGN 2.6.2), so its time in ms is
    the nearest whole millisecond]"""
    return round_half_away(seconds * 1000.0) % _MS_PER_WEEK


def epoch_fields(system: System, t: GnssTime, leap_seconds: int | None) -> dict:
    """The epoch time of an MSM header of ``system`` at GPS time ``t``, as the
    ``epoch``, ``glonass_day`` and ``glonass_epoch`` fields of ``gnsswire.rtcm3.Msm``;
    None for GLONASS without the leap seconds.

    [spec: RTCM 10403.3 3.5.16.4.2 DF004 (GPS ms of week), 3.5.16.8.2 DF428 (QZSS, the
    GPS week: rule T1), 3.5.16.6.2 DF248 (Galileo, ms of the GST week, Galileo OS SIS
    ICD 2.2 5.1.2), 3.5.15.8.2 DF427 (BeiDou, ms of the BDT week, BDS ICD 3.3), 3.5.16.5.2
    Table 3.5-94 DF416 (GLONASS day of week from Sunday) and DF034 (ms of the GLONASS
    day), GLONASS time UTC(SU) + 3 h (GLONASS ICD 3.3.3)]
    """
    week, sow = t.week, t.second + t.fraction
    if system == System.GLONASS:
        if leap_seconds is None:
            return None
        days, sod = timescale.gps_to_glonass(week, sow, leap_seconds)
        ms = round_half_away(sod * 1000.0)
        day = days % 7                    # day 0 of the count is Sunday 1980-01-06
        if ms >= _MS_PER_DAY:
            ms -= _MS_PER_DAY
            day = (day + 1) % 7
        return dict(epoch=None, glonass_day=day, glonass_epoch=ms)
    if system == System.GALILEO:
        week, sow = timescale.gps_to_gst(week, sow)
    elif system == System.BEIDOU:
        week, sow = timescale.gps_to_bdt(week, sow)
    return dict(epoch=_milliseconds(week, sow), glonass_day=None, glonass_epoch=None)
