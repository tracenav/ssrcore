"""Which messages are produced at an epoch, and in which order.

DESIGN.md section 2.6.2, with the station messages of 2.6 and 2.6.4: 1005
(the ARP), 1033 (antenna and receiver descriptors) and, when GLONASS is
sent, 1230 (GLONASS code-phase biases), all as ``gnsswire.rtcm3``
dataclasses; then the MSMs of the epoch, GPS, GLONASS, Galileo, BeiDou,
QZSS, with the multiple message bit set on all but the last.

The station messages go out when due rather than every epoch, the
indicators of 1005 follow the systems sent, and the multiple message bit is
cleared on the MSM actually sent last.
"""

from __future__ import annotations

import dataclasses
import math
from fractions import Fraction
from typing import Iterable, Sequence

from gnsswire import rtcm3

from ssrcore import __version__
from ssrcore.model.types import GnssTime, System, VrsConfig, WEEK_SECONDS
from ssrcore.vrs.placement import Station

# [spec: RTCM 10403.3 Table 3.5-71, the first three digits of an MSM's number name
# its GNSS; design: DESIGN 2.6.2, GPS first (RTCM 3.1.3 item 3, 3.1.7), then
# GLONASS, Galileo, BeiDou, QZSS] MSM system prefix to (send rank, system).
MSM_SYSTEMS = {107: (0, System.GPS), 108: (1, System.GLONASS), 109: (2, System.GALILEO),
               112: (3, System.BEIDOU), 111: (4, System.QZSS)}

# [spec: RTCM 10403.3 3.1.6, network observations referred to the
# ADVNULLANTENNA, which the antenna descriptor then states; Table 3.4-1 DF031,
# 0 "use standard IGS model"; 3.5.12.5, an unknown serial number has a zero
# counter; design: DESIGN 2.6, the receiver is this program, named with its
# version in the firmware field (3.5.12.5: firmware strings follow the
# maker's convention)]
ANTENNA = b"ADVNULLANTENNA"
ANTENNA_SETUP = 0
ANTENNA_SERIAL = b""
RECEIVER = b"ssrcore"
FIRMWARE = __version__.encode("ascii")
RECEIVER_SERIAL = b""

# [spec: RTCM 10403.3 Table 3.4-1 DF021, "reserved ... all bits should be set to
# zero for now"] ITRF realisation year.
ITRF_YEAR = 0
# [spec: RTCM 10403.3 Table 3.4-1 DF141, 1: non-physical or computed reference station]
NON_PHYSICAL = 1
# [spec: RTCM 10403.3 3.1.4 and 3.5.16.2, DF142 applies to MSMs; design: DESIGN
# 2.6.2, every system's observations are synthesised at one instant of GPS
# time, one clock for all] single receiver oscillator.
SINGLE_OSCILLATOR = 1
# [spec: RTCM 10403.3 3.5.16.2, "The Quarter Cycle Indicator (DF364) ... is NOT
# applied to MSMs"; Table 3.4-1 DF364, 00 "correction status unspecified"]
QUARTER_CYCLE = 0

# [spec: SPARTN 8.9, "GLONASS biases are aligned to a theoretical zero bias
# receiver"; RTCM 10403.3 Table 3.4-1 DF421, 1: pseudorange and phase range
# aligned to the same measurement epoch; DF422, the mask from DF423 (MSb) to
# DF426; SPARTN Table 6.2 SF026, SF028, the GLONASS signals with biases are 1C
# and 2C, the C/A codes] 1230: aligned, with zero biases for L1 C/A and L2 C/A.
GLONASS_ALIGNED = 1
GLONASS_MASK = 0b1010

# [spec: RTCM 10403.3 3.1.4, epoch times aligned to GNSS system time within
# 1 µs] tolerance of an epoch time on the observation grid, s.
EPOCH_TOLERANCE = 1e-6


def msm_system(msm: rtcm3.Msm) -> System:
    """[spec: RTCM 10403.3 Table 3.5-71] the GNSS of an MSM from its message number.

    SBAS (110) is refused: SPARTN carries no SBAS corrections (SPARTN Table 6.1).
    """
    entry = MSM_SYSTEMS.get(msm.number // 10)
    if entry is None or not 1 <= msm.number % 10 <= 7:
        raise ValueError(f"message {msm.number} is not an MSM of GPS, GLONASS, Galileo, BeiDou or QZSS")
    return entry[1]


def message_1005(station: Station, systems: Iterable[System]) -> rtcm3.Msg1005:
    """The 1005 of a station: its ARP, ID and the systems sent.

    [spec: RTCM 10403.3 Table 3.5-6; Table 3.4-1 DF003, DF021, DF022 to DF024
    ("service supported"), DF025 to DF027, DF141, DF142, DF364; design:
    DESIGN 2.6.4, DF022 to DF024 the systems of the epoch's MSMs]
    """
    sent = set(systems)
    x, y, z = station.arp
    return rtcm3.Msg1005(number=1005, station=station.station_id, itrf_year=ITRF_YEAR,
                         gps=int(System.GPS in sent), glonass=int(System.GLONASS in sent),
                         galileo=int(System.GALILEO in sent), reference_station=NON_PHYSICAL,
                         x=x, single_oscillator=SINGLE_OSCILLATOR, reserved=0, y=y,
                         quarter_cycle=QUARTER_CYCLE, z=z)


def message_1033(station: Station) -> rtcm3.Msg1033:
    """The 1033 of a station: antenna ADVNULLANTENNA, setup 0, no serial; receiver ssrcore.

    [spec: RTCM 10403.3 Table 3.5-36, 3.5.12.5, 3.1.6; Table 3.1-2 note 1,
    1033 preferred over 1008]
    """
    return rtcm3.Msg1033(number=1033, station=station.station_id,
                         antenna_count=len(ANTENNA), antenna=ANTENNA, antenna_setup=ANTENNA_SETUP,
                         antenna_serial_count=len(ANTENNA_SERIAL), antenna_serial=ANTENNA_SERIAL,
                         receiver_count=len(RECEIVER), receiver=RECEIVER,
                         firmware_count=len(FIRMWARE), firmware=FIRMWARE,
                         receiver_serial_count=len(RECEIVER_SERIAL), receiver_serial=RECEIVER_SERIAL)


def message_1230(station: Station) -> rtcm3.Msg1230:
    """The 1230 of a station: aligned, L1 C/A and L2 C/A code-phase biases zero.

    [spec: RTCM 10403.3 Table 3.5-109; Table 3.4-1 DF421 to DF426; SPARTN 8.9]
    """
    return rtcm3.Msg1230(number=1230, station=station.station_id, aligned=GLONASS_ALIGNED, reserved=0,
                         mask=GLONASS_MASK, l1_ca=0, l1_p=None, l2_ca=0, l2_p=None)


def station_messages(station: Station, systems: Iterable[System]) -> list:
    """1005, 1033 and, when GLONASS is among ``systems``, 1230, in that order.

    [spec: RTCM 10403.3 Table 3.1-2, GNSS RTK: 1005 or 1006, 1033, and 1230
    "if GLONASS data is transmitted"; design: DESIGN 2.6.2 items 1 and 2]
    """
    sent = set(systems)
    messages = [message_1005(station, sent), message_1033(station)]
    if System.GLONASS in sent:
        messages.append(message_1230(station))
    return messages


def order_msms(msms: Sequence[rtcm3.Msm], station_id: int) -> list[rtcm3.Msm]:
    """The MSMs of one epoch in send order, each with the station's ID and the
    multiple message bit set on all but the last.

    MSMs of one system keep their given order (a system whose cells exceed
    the 64 of DF396 needs more than one).
    [spec: RTCM 10403.3 Table 3.4-1 DF393, 1 when more MSMs follow for the
    same time and station, 0 on the last; DF003, the ID links observations to
    the station messages; 3.1.3 item 1, one message per system "if
    possible"; design: DESIGN 2.6.2 item 3 for the order of systems]
    """
    for m in msms:
        msm_system(m)
    ordered = sorted(msms, key=lambda m: MSM_SYSTEMS[m.number // 10][0])   # stable: keeps a system's order
    last = len(ordered) - 1
    return [dataclasses.replace(m, station=station_id, multiple=0 if k == last else 1)
            for k, m in enumerate(ordered)]


class Schedule:
    """The epoch grid and the times the station messages are due.

    [design: DESIGN 2.6.2]

    - Observation epochs lie on whole multiples of 1/``rate`` seconds of GPS
      time (counted from the GPS epoch, so the grid is the same in every
      week), within 1 µs (RTCM 3.1.4).
    - The station messages are due at the first epoch, at every new
      placement, and otherwise at the epoch after which the next would come
      later than ``station_interval`` after the last sending: so no gap
      exceeds the interval, and with an interval that is a multiple of the
      observation period they go out exactly every interval.
    - Epochs strictly increase. [design: a stream of observations does not go
      back in time; a replay starts a new VirtualBase]
    """

    def __init__(self, config: VrsConfig) -> None:
        """[spec: RTCM 10403.3 3.1.7, 0.5 to 2 Hz (checked by VrsConfig)]"""
        self._rate = Fraction(config.rate).limit_denominator(1000)
        self._period = 1 / self._rate
        self._interval = config.station_interval
        self._last_epoch: GnssTime | None = None
        self._last_station: GnssTime | None = None

    @property
    def period(self) -> float:
        """[spec: RTCM 10403.3 3.1.7] seconds between observation epochs, 1/rate."""
        return float(self._period)

    def is_epoch(self, t: GnssTime) -> bool:
        """[design: DESIGN 2.6.2, epochs on whole multiples of the period in GPS time;
        RTCM 10403.3 3.1.4 for the 1 µs]"""
        whole = Fraction(t.week * WEEK_SECONDS + t.second)
        offset = float(whole % self._period) + t.fraction
        period = float(self._period)
        offset = math.fmod(offset, period)
        return min(offset, period - offset) <= EPOCH_TOLERANCE

    def next_epoch(self, t: GnssTime) -> GnssTime:
        """[design: DESIGN 2.6.2] the first epoch of the grid at or after ``t``."""
        whole = Fraction(t.week * WEEK_SECONDS + t.second)
        offset = float(whole % self._period) + t.fraction
        period = float(self._period)
        steps = math.ceil((offset - EPOCH_TOLERANCE) / period)
        epoch = t + (steps * period - offset)
        return epoch

    def advance(self, t: GnssTime) -> None:
        """Take ``t`` as the next epoch; ValueError when it is off the grid or not after the last.

        [design: DESIGN 2.6.2, see the class]
        """
        if not self.is_epoch(t):
            raise ValueError(f"{t!r} is not on the {self.period} s epoch grid")
        if self._last_epoch is not None and t - self._last_epoch <= 0.0:
            raise ValueError(f"epoch {t!r} does not follow {self._last_epoch!r}")
        self._last_epoch = t

    def station_due(self, t: GnssTime, new: bool) -> bool:
        """Whether epoch ``t`` (already advanced to) carries the station messages; records it if so.

        [design: DESIGN 2.6.1 (on every move) and 2.6.2 (every
        ``station_interval``), see the class; RTCM 10403.3 3.5.16.2, a station
        that has moved sends its position before its observations]
        """
        if self._last_epoch is None or t - self._last_epoch != 0.0:
            raise ValueError(f"advance to {t!r} before asking whether the station messages are due")
        due = (new or self._last_station is None
               or (t - self._last_station) + float(self._period) > self._interval + EPOCH_TOLERANCE)
        if due:
            self._last_station = t
        return due
