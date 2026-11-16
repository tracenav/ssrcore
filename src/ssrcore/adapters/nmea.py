"""gnsswire.nmea messages: ZDA to UTC date and time, GGA to a rover position.

DESIGN.md section 2.2. gnsswire keeps every NMEA field as the string on the
wire; this module gives each its meaning. NMEA 0183 itself is not public, so
the meaning is that of the u-blox F9 HPG 1.51 Interface Description, which
gnsswire's sentence tables also follow: sections 2.5.3 (latitude and
longitude), 2.5.5 (position fix flags), 2.7.5.1 (GGA) and 2.7.21.1 (ZDA). The
u-blox X20 HPG 2.11 Interface Description defines the same in the same
sections, but for ZDA, which is its 2.7.20.1.

A sentence whose fields are malformed raises ``ValueError``; a well-formed
sentence that carries nothing usable (no fix, no time yet) gives ``None``.
[design: the engine counts the first as a failure and simply has nothing to
store for the second]
"""

from __future__ import annotations

import calendar
import math
import re

from gnsswire import nmea

from ssrcore.model.types import RoverPosition, UtcTime

# [spec: u-blox F9 HPG 1.51 2.7.5.1, 2.7.21.1; u-blox X20 HPG 2.11 2.7.5.1, 2.7.20.1;
#  "hhmmss.ss"] hours, minutes, seconds with an optional fraction.
_TIME = re.compile(r"(\d{2})(\d{2})(\d{2}(?:\.\d+)?)")
# [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.5.3, "(d)ddmm.mmmmm"] degrees,
# then minutes with a fraction.
_LATITUDE = re.compile(r"(\d{2})(\d{2}(?:\.\d*)?)")
_LONGITUDE = re.compile(r"(\d{3})(\d{2}(?:\.\d*)?)")
# [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.7.5.1, alt and sep "numeric", m]
# a signed decimal.
_NUMBER = re.compile(r"[+-]?\d+(?:\.\d*)?")

# [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.5.5, GGA quality: 1 autonomous
#  GNSS fix, 2 differential GNSS fix, 4 RTK fixed, 5 RTK float] The other values the
# section defines are 0 (no fix, also "GNSS fix, but user limits exceeded")
# and 6 (dead reckoning, also "dead reckoning fix, but user limits
# exceeded"); neither is a position measured within the receiver's limits.
GGA_FIX_QUALITIES = frozenset({1, 2, 4, 5})

# [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.7.5.1, altUnit and sepUnit
#  "M (meters, fixed field)"]
_METRES = "M"


def time_of_day(text: str) -> float:
    """Seconds of the UTC day of an ``hhmmss.ss`` field.

    [spec: u-blox F9 HPG 1.51 2.7.5.1 and 2.7.21.1, u-blox X20 HPG 2.11 2.7.5.1 and
    2.7.20.1, "hhmmss.ss", UTC; IS-GPS-200N with IRN-IS-200N-003 and -004
    20.3.3.5.2.4 b, a positive leap second is the second 23:59:60, so in the last
    minute of the day seconds run to below 61]
    """
    hours, minutes, seconds = _hms(text)
    return hours * 3600.0 + minutes * 60.0 + seconds


def _hms(text: str) -> tuple[int, int, float]:
    m = _TIME.fullmatch(text)
    if m is None:
        raise ValueError(f"UTC time {text!r} is not hhmmss.ss")
    hours, minutes, seconds = int(m[1]), int(m[2]), float(m[3])
    leap = hours == 23 and minutes == 59
    if hours > 23 or minutes > 59 or seconds >= (61.0 if leap else 60.0):
        raise ValueError(f"UTC time {text!r} out of range")
    return hours, minutes, seconds


def degrees(text: str, hemisphere: str, longitude: bool) -> float:
    """Signed decimal degrees of an NMEA latitude or longitude and its hemisphere letter.

    [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.5.3, degrees and minutes,
    the minutes divided by 60 (its example: 4722.80340 N is 47.38005667 degrees);
    2.7.5.1, NS and EW indicators]
    """
    m = (_LONGITUDE if longitude else _LATITUDE).fullmatch(text)
    if m is None:
        raise ValueError(f"{'longitude' if longitude else 'latitude'} {text!r} is not "
                         f"{'dddmm.mmmmm' if longitude else 'ddmm.mmmmm'}")
    whole, minutes = int(m[1]), float(m[2])
    value = whole + minutes / 60.0
    limit = 180.0 if longitude else 90.0
    if minutes >= 60.0 or value > limit:
        raise ValueError(f"{text!r} out of range")
    signs = {"E": 1.0, "W": -1.0} if longitude else {"N": 1.0, "S": -1.0}
    if hemisphere not in signs:
        raise ValueError(f"hemisphere {hemisphere!r} is not one of {sorted(signs)}")
    return signs[hemisphere] * value


def _metres(text: str, unit: str, name: str) -> float:
    if _NUMBER.fullmatch(text) is None or unit != _METRES:
        raise ValueError(f"GGA {name} {text!r} {unit!r} is not a number of metres")
    return float(text)


def _text(value: str | None) -> str:
    # [design: gnsswire decodes a sentence with fewer fields than its table, the
    # missing ones None; a field the sentence does not have is read as an empty one]
    return "" if value is None else value


def rover_position(gga: nmea.Gga) -> RoverPosition | None:
    """The rover position of a GGA sentence; None when the sentence has no usable fix.

    Latitude and longitude in radians, height on the ellipsoid as the
    altitude above mean sea level plus the geoid separation.
    [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.7.5.1 (alt "above mean sea
    level", sep "difference between ellipsoid and mean sea level"), 2.5.3, 2.5.5;
    design: radians, as every angle in the engine (DESIGN 4)]
    """
    t = _text
    if not t(gga.quality).isdigit():
        raise ValueError(f"GGA quality {gga.quality!r} is not a digit")
    if int(gga.quality) not in GGA_FIX_QUALITIES:
        return None
    latitude = degrees(t(gga.lat), t(gga.ns), longitude=False)
    longitude = degrees(t(gga.lon), t(gga.ew), longitude=True)
    height = _metres(t(gga.alt), t(gga.alt_unit), "alt") + _metres(t(gga.sep), t(gga.sep_unit), "sep")
    return RoverPosition(math.radians(latitude), math.radians(longitude), height,
                         time_of_day(t(gga.time)))


def utc_time(zda: nmea.Zda) -> UtcTime | None:
    """The UTC date and time of a ZDA sentence; None when the receiver has no time yet.

    [spec: u-blox F9 HPG 1.51 2.7.21.1 and u-blox X20 HPG 2.11 2.7.20.1: time
    "hhmmss.ss" UTC, day 1 to 31, month 1 to 12, year "yyyy"; design: a receiver without a valid time sends
    the fields null, which is no time rather than an error; the local zone
    fields are ignored, since the time and date fields are UTC whatever they
    hold; a field the sentence does not have (None) is no time, as an empty one]
    """
    fields = (zda.time, zda.day, zda.month, zda.year)
    if not all(fields):
        return None
    if not (len(zda.day) == 2 and len(zda.month) == 2 and len(zda.year) == 4
            and (zda.day + zda.month + zda.year).isdigit()):
        raise ValueError(f"ZDA date {zda.day!r} {zda.month!r} {zda.year!r} is not dd mm yyyy")
    day, month, year = int(zda.day), int(zda.month), int(zda.year)
    # [design: a date the calendar does not have, 31 April say, is malformed]
    if not 1 <= month <= 12 or not 1 <= day <= calendar.monthrange(year, month)[1]:
        raise ValueError(f"ZDA date {year:04d}-{month:02d}-{day:02d} does not exist")
    hour, minute, second = _hms(zda.time)
    return UtcTime(year, month, day, hour, minute, second)
