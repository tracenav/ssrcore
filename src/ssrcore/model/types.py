"""Value types shared by every layer above the kernels: satellite, signal, time,
configuration, the typed inputs of the engine and the result of feeding one.

DESIGN.md sections 2.2 and 2.5. Every type is a frozen dataclass; those that
serve as keys (``Sat``, ``Signal``, ``GnssTime``) and every input (so that
rule S8 can compare messages) are hashable. Values are in SI units except
where a field says otherwise and why. Nothing here computes a correction:
scaling of wire integers is the adapters' (``ssrcore.adapters``), arithmetic
the toolkit's.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import IntEnum, StrEnum
from typing import Any

# [spec: SPARTN Table 8.0] speed of light, m/s: one definition, the
# correction kernel's (DESIGN 4), re-exported here for the model.
from ssrcore.toolkit.ssr import SPEED_OF_LIGHT

# [spec: SPARTN Table 6.1] the constellations, numbered as SPARTN numbers its
# subtypes: one enumeration, the toolkit's (the kernels' SSR_SYSTEM_ codes),
# re-exported here for the model.
# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] one week, s.
from ssrcore.toolkit.timescale import WEEK_SECONDS, System


# ---- satellites -----------------------------------------------------------------------


# [spec: SPARTN Table 6.2 SF011, SF012, SF093, SF094, SF095] the widest satellite
# mask of each system; position + 1 is the PRN (GLONASS: the orbit slot), so
# these are the highest satellite numbers SPARTN can name.
_MAX_SATELLITE = {System.GPS: 64, System.GLONASS: 63, System.GALILEO: 64,
                  System.BEIDOU: 64, System.QZSS: 64}


@dataclass(frozen=True, order=True, slots=True)
class Sat:
    """A satellite: its system and its number, PRN or (GLONASS) orbit slot.

    Numbers are those SPARTN uses: satellite mask position + 1 (Table 6.2
    SF011, SF012, SF093 to SF095), so QZSS satellites are 1, 2, ... and not
    the PRN 193, 194, ... of their ranging codes.
    """

    system: System
    number: int

    def __post_init__(self) -> None:
        """[spec: SPARTN Table 6.2 SF011, SF012, SF093 to SF095] number from 1 to the widest mask."""
        object.__setattr__(self, "system", System(self.system))
        if not isinstance(self.number, int) or isinstance(self.number, bool):
            raise TypeError(f"satellite number must be an int, not {type(self.number).__name__}")
        if not 1 <= self.number <= _MAX_SATELLITE[self.system]:
            raise ValueError(f"{self.system.name} satellite number {self.number} outside "
                             f"1 to {_MAX_SATELLITE[self.system]}")

    @classmethod
    def gps(cls, prn: int) -> Sat:
        """[spec: SPARTN Table 6.2 SF011] a GPS satellite by PRN."""
        return cls(System.GPS, prn)

    @classmethod
    def glonass(cls, slot: int) -> Sat:
        """[spec: SPARTN Table 6.2 SF012] a GLONASS satellite by orbit slot."""
        return cls(System.GLONASS, slot)

    @classmethod
    def galileo(cls, prn: int) -> Sat:
        """[spec: SPARTN Table 6.2 SF093] a Galileo satellite by PRN."""
        return cls(System.GALILEO, prn)

    @classmethod
    def beidou(cls, prn: int) -> Sat:
        """[spec: SPARTN Table 6.2 SF094] a BeiDou satellite by PRN."""
        return cls(System.BEIDOU, prn)

    @classmethod
    def qzss(cls, number: int) -> Sat:
        """[spec: SPARTN Table 6.2 SF095] a QZSS satellite by its SPARTN mask number."""
        return cls(System.QZSS, number)


# ---- time -----------------------------------------------------------------------------


@dataclass(frozen=True, order=True, slots=True, init=False)
class GnssTime:
    """A time on the GPS time scale: full week number since 1980-01-06 (not
    rolled over), whole second of week in [0, 604800) and the fraction of
    that second in [0, 1).

    Times on other scales are converted to GPS time on input (rule T1), so
    every time inside the engine is one of these. Differences are taken in
    seconds before any further arithmetic (DESIGN 4): ``t1 - t2`` is a
    float, ``t + dt`` a ``GnssTime``. ``seconds`` is the second of week as
    one double, for the kernels. [design: DESIGN 4, time; a double second of
    week near 604800 has a spacing of 2^-33 s, 35 mm of light travel, so the
    whole seconds are kept apart and the fraction keeps about 10^-16 s]
    """

    week: int
    second: int
    fraction: float

    def __init__(self, week: int, seconds: float, fraction: float = 0.0) -> None:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] week from the GPS epoch, seconds
        within the week; ``seconds`` may carry a fraction, and ``fraction`` adds to it.
        [design: DESIGN 4, the whole seconds split off exactly]"""
        if not isinstance(week, int) or isinstance(week, bool):
            raise TypeError(f"week must be an int, not {type(week).__name__}")
        whole, frac = _split(seconds)
        extra, frac = _split(frac + float(fraction))
        whole += extra
        if week < 0:
            raise ValueError(f"week {week} before the GPS epoch")
        if not 0 <= whole < _WEEK:
            raise ValueError(f"seconds of week {whole + frac} outside [0, 604800)")
        object.__setattr__(self, "week", week)
        object.__setattr__(self, "second", whole)
        object.__setattr__(self, "fraction", frac)

    @property
    def seconds(self) -> float:
        """[design: DESIGN 4] the second of week as one double, for the kernels,
        whose arithmetic multiplies time differences only by rates."""
        return self.second + self.fraction

    def __repr__(self) -> str:
        """[design: readable times in test reports]"""
        return f"GnssTime(week={self.week}, second={self.second}, fraction={self.fraction!r})"

    @classmethod
    def normalised(cls, week: int, seconds: float, fraction: float = 0.0) -> GnssTime:
        """The same instant with seconds brought into [0, 604800) by whole weeks.

        [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, beginning or end of week
        crossovers]
        """
        whole, frac = _split(seconds)
        extra, frac = _split(frac + float(fraction))
        weeks, whole = divmod(whole + extra, _WEEK)
        return cls(week + weeks, whole, frac)

    def __add__(self, dt: float) -> GnssTime:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] this time moved by dt
        seconds, across weeks.

        [design: DESIGN 4, the whole seconds of ``dt`` are added exactly and its
        fraction to the fraction, so no resolution is lost]
        """
        if isinstance(dt, GnssTime) or not isinstance(dt, (int, float)):
            return NotImplemented
        whole, frac = _split(dt)
        return GnssTime.normalised(self.week, self.second + whole, self.fraction + frac)

    def __sub__(self, other: GnssTime | float):
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] seconds from ``other`` to
        this time, or this time moved back.

        Weeks and whole seconds are differenced as integers and the
        fractions apart, so that the result is exact up to its own rounding.
        """
        if isinstance(other, GnssTime):
            whole = (self.week - other.week) * _WEEK + (self.second - other.second)
            return whole + (self.fraction - other.fraction)
        if isinstance(other, (int, float)):
            return self + (-other)
        return NotImplemented


# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] one week in whole seconds.
_WEEK = 604800


def _split(value: float) -> tuple[int, float]:
    """[design: DESIGN 4] (whole, fraction) with fraction in [0, 1), both exact:
    the fraction of a double is exact by Sterbenz's lemma."""
    if isinstance(value, int):
        return value, 0.0
    v = float(value)
    if v != v or v in (float("inf"), float("-inf")):
        raise ValueError(f"time {value} is not finite")
    whole = math.floor(v)
    frac = v - whole
    if frac >= 1.0:   # cannot happen for a finite double, kept as a guard
        whole, frac = whole + 1, 0.0
    return int(whole), frac


# ---- signals --------------------------------------------------------------------------

# Carrier frequency per system and RINEX 3 band digit, Hz (exact integers).
# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.1.1 (L1 = 154 x 10.23 MHz,
#  L2 = 120 x 10.23 MHz);
#  RTCM Table 3.5-73 (GPS L5; QZSS L1, L2, L5; BeiDou B1, B2, B3);
#  Galileo OS SIS ICD 2.2 Table 2 (E1, E6, E5, E5a, E5b); BDS ICD 4.2.1 (B1I);
#  TM-23 2.2.4 Table 2.12 (BeiDou B1C at 1575.42 MHz, B2a+b at 1191.795 MHz);
#  TM-23 2.2.3 (BeiDou B2a shares the E5a band, B2b the E5b band)]
# Band digits are those of RINEX 3, which SPARTN Table 6.2 SF102 to SF107
# names for its signals; RTCM Table 3.5-73 lists the same digits as its
# reference signals (L1C, L2P, L5I; L2I, L7I, L6I for BeiDou; L8I for E5).
_CARRIER_HZ = {
    (System.GPS, "1"): 1575420000, (System.GPS, "2"): 1227600000, (System.GPS, "5"): 1176450000,
    (System.GALILEO, "1"): 1575420000, (System.GALILEO, "5"): 1176450000,
    (System.GALILEO, "7"): 1207140000, (System.GALILEO, "8"): 1191795000,
    (System.GALILEO, "6"): 1278750000,
    (System.BEIDOU, "2"): 1561098000, (System.BEIDOU, "1"): 1575420000,
    (System.BEIDOU, "5"): 1176450000, (System.BEIDOU, "7"): 1207140000,
    (System.BEIDOU, "8"): 1191795000, (System.BEIDOU, "6"): 1268520000,
    (System.QZSS, "1"): 1575420000, (System.QZSS, "2"): 1227600000, (System.QZSS, "5"): 1176450000,
}
# [spec: GLONASS ICD 3.3.1.1] f = f0 + K * delta f per sub-band, Hz.
_GLONASS_HZ = {"1": (1602000000, 562500), "2": (1246000000, 437500)}
# [spec: GLONASS ICD 3.3.1.1 and Table 3.1] frequency channel numbers K.
GLONASS_CHANNELS = range(-7, 7)

# [spec: SPARTN Table 6.2 SF025 to SF028, SF102 to SF107] the signals of each
# system's bias masks, in mask order (position 0 first). Phase bias masks name
# the L code and code bias masks the C code of the same signal at the same
# position; positions after these are spare. Galileo 8Q (E5) and 6C (E6) at
# positions 3 and 4 are those SPARTN v2.0.3 defines and v2.0.2 left spare.
SPARTN_SIGNALS = {
    System.GPS: ("1C", "2W", "2L", "5Q"),
    System.GLONASS: ("1C", "2C"),
    System.GALILEO: ("1C", "5Q", "7Q", "8Q", "6C"),
    System.BEIDOU: ("2I", "5P", "7I", "6I", "1P", "7P", "8P"),
    System.QZSS: ("1C", "2L", "5Q"),
}


@dataclass(frozen=True, order=True, slots=True)
class Signal:
    """A signal a receiver tracks: system and RINEX 3 code without the
    observation type (band digit and attribute letter, ``"1C"``).

    The signals SPARTN carries biases for are class attributes named as in
    DESIGN 2.1, ``Signal.GPS_L1C``, ``Signal.BEIDOU_L5P`` and so on. Any other
    code of a known band may be named (``Signal(System.GPS, "2X")``) so that
    the caller can ask for it; it has a frequency but no SPARTN bias, and is
    not substituted for another code of the same band (rule C9).
    """

    system: System
    code: str

    def __post_init__(self) -> None:
        """[spec: SPARTN Table 6.2 SF102 to SF107, RINEX 3 observation codes; design: any
        attribute letter A to Z is accepted, since the attribute list of RINEX is not
        staged and an unknown attribute only means that no bias is found for it]"""
        object.__setattr__(self, "system", System(self.system))
        code = self.code
        if not (isinstance(code, str) and len(code) == 2 and "A" <= code[1] <= "Z"):
            raise ValueError(f"RINEX signal code must be a band digit and an attribute letter, not {code!r}")
        if not self._known_band(self.system, code[0]):
            raise ValueError(f"{self.system.name} has no band {code[0]!r}")

    @staticmethod
    def _known_band(system: System, band: str) -> bool:
        if system == System.GLONASS:
            return band in _GLONASS_HZ
        return (system, band) in _CARRIER_HZ

    @property
    def band(self) -> str:
        """[spec: RTCM Table 3.5-73, RINEX band digit of the reference signals] the band digit."""
        return self.code[0]

    @property
    def code_observation(self) -> str:
        """[spec: SPARTN Table 6.2 SF027, SF028, SF105 to SF107] the pseudorange code, ``"C1C"``."""
        return "C" + self.code

    @property
    def phase_observation(self) -> str:
        """[spec: SPARTN Table 6.2 SF025, SF026, SF102 to SF104] the carrier phase code, ``"L1C"``."""
        return "L" + self.code

    @property
    def spartn_position(self) -> int | None:
        """[spec: SPARTN Table 6.2 SF025 to SF028, SF102 to SF107] position in the bias masks, or None."""
        codes = SPARTN_SIGNALS[self.system]
        return codes.index(self.code) if self.code in codes else None

    @classmethod
    def from_spartn(cls, system: System, position: int) -> Signal | None:
        """[spec: SPARTN Table 6.2 SF025 to SF028, SF102 to SF107] the signal at a bias mask
        position; None for a spare position."""
        codes = SPARTN_SIGNALS[System(system)]
        return cls(system, codes[position]) if 0 <= position < len(codes) else None

    def frequency(self, channel: int | None = None) -> float:
        """Carrier frequency in Hz; GLONASS needs the satellite's frequency channel K.

        [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.1.1; GLONASS ICD 3.3.1.1; Galileo OS
        SIS ICD 2.2 Table 2; BDS ICD 4.2.1; RTCM Table 3.5-73; TM-23 2.2.3, 2.2.4]
        """
        if self.system == System.GLONASS:
            if channel is None:
                raise ValueError("a GLONASS frequency needs the frequency channel")
            if channel not in GLONASS_CHANNELS:
                raise ValueError(f"GLONASS frequency channel {channel} outside -7 to 6")
            f0, df = _GLONASS_HZ[self.band]
            return float(f0 + channel * df)
        return float(_CARRIER_HZ[(self.system, self.band)])

    def wavelength(self, channel: int | None = None) -> float:
        """[spec: SPARTN Table 8.0 (c)] carrier wavelength in metres, c / f."""
        return SPEED_OF_LIGHT / self.frequency(channel)


for _system, _codes in SPARTN_SIGNALS.items():
    for _code in _codes:
        setattr(Signal, f"{_system.name}_L{_code}", Signal(_system, _code))
del _system, _codes, _code


# ---- configuration --------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MaxAge:
    """Maximum age per kind of correction, seconds (DESIGN 2.5).

    [design: SPARTN gives no age limit for these; the defaults are bounded by
    the 510 s of SF005 and the 450 s of SF068]
    """

    orbit: float = 120.0
    clock: float = 30.0
    bias: float = 300.0
    atmosphere: float = 120.0

    def __post_init__(self) -> None:
        """[spec: SPARTN Table 6.2 SF005, SF068] each limit positive and below the
        issue-of-update lifetime (510 s; 450 s for atmosphere)."""
        for name, limit in (("orbit", 510.0), ("clock", 510.0), ("bias", 510.0), ("atmosphere", 450.0)):
            value = getattr(self, name)
            if not 0.0 < value < limit:
                raise ValueError(f"max_age.{name} = {value} outside (0, {limit}) s")


class MsmType(IntEnum):
    """The MSM the virtual base station sends. [spec: RTCM Table 3.1-2, MSM7 the full
    high-precision service; design: DESIGN 2.6.2, MSM4 kept for receivers and links
    that take the standard-resolution full message without rates]"""

    MSM4 = 4
    MSM7 = 7


@dataclass(frozen=True, slots=True)
class VrsConfig:
    """Parameters of the virtual base station (DESIGN 2.6), each with its default."""

    # [design: DESIGN 2.6.1, the remaining baseline error at 5 km is well below a wide-lane cycle]
    move_distance: float = 5000.0
    # [design: DESIGN 2.6.1, no move on a single bad position; a time, so the same at every rate]
    # seconds the rover must have been beyond move_distance, counted in epochs with a position.
    move_confirm: float = 10.0
    # [design: DESIGN 2.6.1] ECEF metres; None lets the station follow the rover.
    pinned_position: tuple[float, float, float] | None = None
    # [spec: RTCM 3.1.7, 0.5 to 2 Hz for network observations] epochs per second.
    rate: float = 1.0
    # [design: DESIGN 2.6.2, within the 15 s of RTCM 3.1.7] seconds between 1005, 1033, 1230.
    station_interval: float = 10.0
    # [spec: RTCM Table 3.1-2, MSM7 the full high-precision service; design: DESIGN 2.6.2]
    msm: MsmType = MsmType.MSM7
    # [design: DESIGN 2.6.3] radians.
    elevation_mask: float = math.radians(5.0)
    # [spec: RTCM DF003, 0 to 4095]
    station_id: int = 0
    # [design: DESIGN Q21, the receiver decides; the record carries the fix flag]
    float_phase: bool = True

    def __post_init__(self) -> None:
        """[spec: RTCM 3.1.7 (rate), DF003 (station ID); design: DESIGN 2.6 for the others]"""
        object.__setattr__(self, "msm", MsmType(self.msm))
        if not 0.5 <= self.rate <= 2.0:
            raise ValueError(f"rate {self.rate} Hz outside 0.5 to 2 Hz (RTCM 3.1.7)")
        if not 0 <= self.station_id <= 4095:
            raise ValueError(f"station_id {self.station_id} outside 0 to 4095 (DF003)")
        if not 0.0 < self.station_interval <= 15.0:
            raise ValueError(f"station_interval {self.station_interval} s outside (0, 15] (RTCM 3.1.7)")
        if self.move_distance <= 0.0 or not self.move_confirm >= 0.0:
            raise ValueError("move_distance must be positive and move_confirm not negative")
        if not 0.0 <= self.elevation_mask < math.pi / 2:
            raise ValueError(f"elevation_mask {self.elevation_mask} rad outside [0, pi/2)")
        if self.pinned_position is not None:
            object.__setattr__(self, "pinned_position", tuple(float(v) for v in self.pinned_position))
            if len(self.pinned_position) != 3:
                raise ValueError("pinned_position is an ECEF (x, y, z)")


@dataclass(frozen=True, slots=True)
class SpartnSource:
    """An atmosphere source of ``Config.atmosphere``: the SPARTN store. ``stream`` None
    is the satellite stream, bound by its SIOU (rule A1); a named stream is that
    stream's atmosphere under its own SIOU (rule A10).
    [design: DESIGN 2.5, rule A9; a value, so that ``Config`` stays frozen and replayable]"""

    stream: StreamId | None = None

    def __post_init__(self) -> None:
        """[design: DESIGN 2.5, a stream is a ``StreamId`` (rule S1)]"""
        if self.stream is not None and not isinstance(self.stream, StreamId):
            raise TypeError(f"SpartnSource.stream must be a StreamId, not {type(self.stream).__name__}")


@dataclass(frozen=True, slots=True)
class ModelSource:
    """An atmosphere source of ``Config.atmosphere``: an a-priori model; ``"unb3m"``, the
    only one, gives the troposphere of UNB3m. [design: DESIGN 2.5, rule A9]"""

    model: str = "unb3m"

    def __post_init__(self) -> None:
        """[design: DESIGN 2.7, UNB3m is the one model the toolkit has]"""
        if self.model != "unb3m":
            raise ValueError(f"no atmosphere model {self.model!r}; the toolkit has 'unb3m'")


@dataclass(frozen=True, slots=True)
class CallerSource:
    """An atmosphere source of ``Config.atmosphere``: the provider object the caller
    passes as ``Engine(config, providers={name: provider})``, whose ``name`` is this
    one. [design: DESIGN 2.1, the object is not part of the configuration, its name is]"""

    name: str

    def __post_init__(self) -> None:
        """[design: DESIGN 2.5, the name keys the caller's object and the record's origin]"""
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("CallerSource.name must be a non-empty string")


@dataclass(frozen=True, slots=True)
class Config:
    """Everything that configures an engine; every default is here so that a
    replay can reproduce a run (DESIGN 2.5).

    Angles in radians, times in seconds.
    """

    max_age: MaxAge = field(default_factory=MaxAge)
    # [design: DESIGN 2.5, geometry only below this elevation; 0 keeps every satellite above the horizon]
    elevation_cutoff: float = 0.0
    # [design: DESIGN Q7, SF043 is at 0 m on the ellipsoid, so a receiver above it needs the reduction]
    troposphere_height_reduction: bool = True
    # [design: DESIGN 3.3 rule S2, switch only on a clear lead]
    stream_switch: float = 5.0
    # [design: DESIGN 3.3 rule T5, messages held until the time is known]
    hold_limit: int = 256
    vrs: VrsConfig = field(default_factory=VrsConfig)
    # [design: DESIGN 3.3 rule S10] the stream of orbit, clock and biases; None: the active stream (S2).
    satellite_stream: StreamId | None = None
    # [design: DESIGN 3.3 rule A9] atmosphere sources in walk order; the default is the engine before the seam.
    atmosphere: tuple[SpartnSource | ModelSource | CallerSource, ...] = (SpartnSource(),)

    def __post_init__(self) -> None:
        """[design: DESIGN 2.5, values outside these ranges have no meaning]"""
        if not -math.pi / 2 <= self.elevation_cutoff < math.pi / 2:
            raise ValueError(f"elevation_cutoff {self.elevation_cutoff} rad outside [-pi/2, pi/2)")
        if self.stream_switch < 0.0:
            raise ValueError("stream_switch must not be negative")
        if self.hold_limit < 1:
            raise ValueError("hold_limit must be at least 1")
        if self.satellite_stream is not None and not isinstance(self.satellite_stream, StreamId):
            raise TypeError("satellite_stream must be a StreamId or None")
        object.__setattr__(self, "atmosphere", tuple(self.atmosphere))
        names = set()
        for source in self.atmosphere:
            if not isinstance(source, (SpartnSource, ModelSource, CallerSource)):
                raise TypeError(f"not an atmosphere source: {type(source).__name__}")
            if isinstance(source, CallerSource):
                if source.name in names:
                    raise ValueError(f"CallerSource {source.name!r} named twice")
                names.add(source.name)


# ---- inputs: SPARTN -------------------------------------------------------------------
#
# Each SPARTN input carries the transport fields that identify its stream
# (TF010, TF011) and its time tag as sent (TF008, TF009): the engine resolves
# a 16-bit tag and converts the tag to GPS time (rules T1 to T3), because that
# needs the receipt time and the time base, which an adapter does not have.


@dataclass(frozen=True, slots=True)
class SpartnTag:
    """A SPARTN time tag as sent, in the system time of the message's constellation
    (GPS time for GAD).

    [spec: SPARTN Table 7.1 TF008, TF009] ``full`` is TF008: True for the
    32-bit tag (seconds since 2010-01-01 00:00:00), False for the 16-bit tag
    (seconds of the half day, 0 to 43199).
    """

    value: int
    full: bool

    def __post_init__(self) -> None:
        """[spec: SPARTN Table 7.1 TF009] range of each tag type."""
        limit = 1 << 32 if self.full else 43200
        if not 0 <= self.value < limit:
            raise ValueError(f"time tag {self.value} outside 0 to {limit - 1}")


@dataclass(frozen=True, slots=True)
class StreamId:
    """The stream a SPARTN message belongs to (rule S1).

    [spec: SPARTN Table 7.1 TF010 solution ID, TF011 solution processor ID]
    """

    solution: int
    processor: int


@dataclass(frozen=True, slots=True)
class OrbitCorrection:
    """An OCB orbit block (SPARTN Table 6.5) in physical units.

    ``radial``, ``along``, ``cross`` are δO in metres (SF020). ``yaw`` is SF021
    in degrees, as DESIGN 2.4 records it; None when the header's SF008 is 0 or
    the value is the invalid 0x3F (rule C8). ``ephemeris_type`` is the OCB
    header's (SF016, SF017, SF096 to SF098), which with ``iode`` names the
    ephemeris the block corrects (rules E1, E2).

    For a format with rates (rule C11), ``rate`` is δȮ (radial, along, cross)
    in m/s, signed as δO, and ``reference`` its reference time t0 in GPS time;
    δO is then the value at t0. Both None for SPARTN, which sends no rates.
    """

    iode: int
    radial: float
    along: float
    cross: float
    yaw: float | None = None
    ephemeris_type: int = 0
    rate: tuple[float, float, float] | None = None
    reference: GnssTime | None = None

    def __post_init__(self) -> None:
        """[spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6), DF365 (t0); IGS SSR v1.00 4.1, IDF013;
        design: rates need their reference time, and a reference time alone means nothing]"""
        if self.rate is not None:
            object.__setattr__(self, "rate", tuple(float(v) for v in self.rate))
            if len(self.rate) != 3:
                raise ValueError("rate is (radial, along, cross) in m/s")
        if (self.rate is None) != (self.reference is None):
            raise ValueError("rate and reference are given together")
        if self.reference is not None and not isinstance(self.reference, GnssTime):
            raise TypeError("reference must be a GnssTime")


@dataclass(frozen=True, slots=True)
class ClockCorrection:
    """An OCB clock block (SPARTN Table 6.6) in physical units.

    ``correction`` is δC in metres (SF020). ``iode_continuity`` is SF022 as
    seconds (its special values 0, 1, 5, ... 320 s). ``ure`` is the SF024
    index as sent (0 unknown, 7 above 1 m), since its values are classes, not
    a scale.

    For a format with a clock polynomial (rule C11), ``correction`` is C0 at
    the reference time ``reference`` (t0, GPS time), ``rate`` C1 in m/s and
    ``acceleration`` C2 in m/s², signed as ``correction``; a term not sent is
    None and counts as zero. All None for SPARTN, which sends no rates.
    """

    iode_continuity: float
    correction: float
    ure: int
    rate: float | None = None
    acceleration: float | None = None
    reference: GnssTime | None = None

    def __post_init__(self) -> None:
        """[spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8), DF376 (t0); IGS SSR v1.00 4.2, IDF019;
        design: a polynomial needs its reference time, and a reference time alone means nothing]"""
        polynomial = self.rate is not None or self.acceleration is not None
        if polynomial != (self.reference is not None):
            raise ValueError("rate or acceleration and reference are given together")
        if self.reference is not None and not isinstance(self.reference, GnssTime):
            raise TypeError("reference must be a GnssTime")

    @property
    def polynomial(self) -> bool:
        """[design: rule C11, a clock with C1 or C2 is evaluated at a time; SPARTN's is not]"""
        return self.reference is not None


@dataclass(frozen=True, slots=True)
class PhaseBias:
    """One phase bias (SPARTN Table 6.12): metres (SF020), fix flag SF023,
    continuity SF015 as seconds."""

    signal: Signal
    value: float
    fixed: bool
    continuity: float


@dataclass(frozen=True, slots=True)
class CodeBias:
    """One code bias (SPARTN Tables 6.7 to 6.11): metres (SF029)."""

    signal: Signal
    value: float


@dataclass(frozen=True, slots=True)
class OcbSatellite:
    """One satellite block of an OCB message (SPARTN Table 6.4).

    With ``dnu`` set (SF013) the block carries nothing else (rule C6).
    ``continuity`` is SF015 as seconds. ``phase_biases`` and ``code_biases``
    are None when the block has no bias block (SF014), and empty when it has
    one with no signal of a known position.
    """

    sat: Sat
    dnu: bool
    continuity: float | None = None
    orbit: OrbitCorrection | None = None
    clock: ClockCorrection | None = None
    phase_biases: tuple[PhaseBias, ...] | None = None
    code_biases: tuple[CodeBias, ...] | None = None

    def __post_init__(self) -> None:
        """[spec: SPARTN Table 6.2 SF013, no other field of the block with DNU set]"""
        if self.dnu and any(v is not None for v in (self.continuity, self.orbit, self.clock,
                                                     self.phase_biases, self.code_biases)):
            raise ValueError("a satellite block with DNU set carries no other field (SF013)")


@dataclass(frozen=True, slots=True)
class OcbInput:
    """An OCB message, SM 0-0 to 0-4 (SPARTN 6.5.1, Table 6.3), in physical units."""

    system: System
    stream: StreamId
    tag: SpartnTag
    siou: int                 # SF005
    end_of_set: bool          # SF010
    yaw_present: bool         # SF008
    datum: int                # SF009: 0 ITRF, 1 provider defined
    ephemeris_type: int       # SF016, SF017, SF096, SF097 or SF098 as sent
    satellites: tuple[OcbSatellite, ...]


@dataclass(frozen=True, slots=True)
class OcbHeaderInput:
    """The header of an OCB message whose subtype has no table (SPARTN Table 6.1:
    SM 0-5 to 0-15 are TBD), read for its end-of-set flag.

    [spec: SPARTN 8.4 item 5, decoders "should be designed so as to decode the
    SF010 field of all current and future SM 0 message possibilities"; Table
    6.3 SF005, SF010, SF009; Table 7.1 TF007 to TF011]
    """

    subtype: int
    stream: StreamId
    tag: SpartnTag
    siou: int                 # SF005
    end_of_set: bool          # SF010
    datum: int                # SF009


@dataclass(frozen=True, slots=True)
class TroposphereBlock:
    """An HPAC troposphere block (SPARTN Table 6.16) in physical units.

    ``hydrostatic`` is the area's average vertical hydrostatic delay at 0 m
    on the ellipsoid in metres (SF043 with its 2.3 m offset). ``coefficients``
    are those of the equation type SF041, in the order T00, T01, T10, T11 and
    the units SPARTN 8.12 uses (m, m/deg, m/deg²). ``residuals`` are the grid
    residuals in metres, None per node holding the invalid value and None
    altogether without a grid block. ``quality`` is the SF042 index.
    """

    equation_type: int
    quality: int
    hydrostatic: float
    coefficients: tuple[float, ...]
    residuals: tuple[float | None, ...] | None = None


@dataclass(frozen=True, slots=True)
class IonosphereSatellite:
    """An HPAC ionosphere satellite block (SPARTN Table 6.20) in physical units.

    ``coefficients`` C00, C01, C10, C11 as the equation type has them, in
    TECU, TECU/deg, TECU/deg² (SPARTN 8.11); ``residuals`` in TECU as for the
    troposphere; ``quality`` the SF055 index.
    """

    sat: Sat
    quality: int
    coefficients: tuple[float, ...]
    residuals: tuple[float | None, ...] | None = None


@dataclass(frozen=True, slots=True)
class IonosphereBlock:
    """An HPAC ionosphere block (SPARTN Table 6.19): equation type SF054 and satellites."""

    equation_type: int
    satellites: tuple[IonosphereSatellite, ...]


@dataclass(frozen=True, slots=True)
class HpacArea:
    """One area of an HPAC message (SPARTN Tables 6.14, 6.15): area ID SF031,
    grid point count SF039, troposphere and ionosphere when present (SF040)."""

    area_id: int
    grid_points: int
    troposphere: TroposphereBlock | None = None
    ionosphere: IonosphereBlock | None = None


@dataclass(frozen=True, slots=True)
class HpacInput:
    """An HPAC message, SM 1-0 to 1-4 (SPARTN 6.5.2, Table 6.13), in physical units."""

    system: System
    stream: StreamId
    tag: SpartnTag
    siou: int                 # SF005
    aiou: int                 # SF068
    areas: tuple[HpacArea, ...]


@dataclass(frozen=True, slots=True)
class AreaDefinition:
    """A GAD area (SPARTN Table 6.24) in physical units.

    ``latitude`` is the northern-most latitude (SF032) and ``longitude`` the
    western-most longitude (SF033) of the grid, ``*_nodes`` the node counts
    (SF034, SF035) and ``*_spacing`` the node spacing (SF036, SF037), all in
    degrees. [design: degrees, because SPARTN 8.10 and 8.11 compute area
    offsets in degrees (DESIGN 4)]
    """

    area_id: int
    latitude: float
    longitude: float
    latitude_nodes: int
    longitude_nodes: int
    latitude_spacing: float
    longitude_spacing: float


@dataclass(frozen=True, slots=True)
class GadInput:
    """A GAD message, SM 2-0 (SPARTN 6.5.3, Table 6.23), in physical units; its tag is in GPS time."""

    stream: StreamId
    tag: SpartnTag
    siou: int                 # SF005
    aiou: int                 # SF068
    areas: tuple[AreaDefinition, ...]


# ---- terms and the reasons they are absent (DESIGN 2.4) --------------------------------


class Term(StrEnum):
    """The terms of the record that can be absent, by field name.
    [design: DESIGN 2.4, the record states for each absent term why]"""

    SATELLITE = "satellite"
    GEOMETRY = "geometry"
    ORBIT = "orbit"
    CLOCK_CORRECTION = "clock_correction"
    CODE_BIAS = "code_bias"
    PHASE_BIAS = "phase_bias"
    TROPOSPHERE = "troposphere"
    TROPOSPHERE_GRID = "troposphere.grid"
    IONOSPHERE = "ionosphere"
    IONOSPHERE_GRID = "ionosphere.grid"
    EARTH_ROTATION = "modelled.earth_rotation"
    SHAPIRO = "modelled.shapiro"
    RELATIVITY_CLOCK = "modelled.relativity_clock"
    SOLID_TIDE = "modelled.solid_tide"
    WINDUP = "modelled.windup"
    FREQUENCY = "frequency"


class Missing(StrEnum):
    """Why a term is absent, each by the rule of DESIGN that removes it.
    [design: DESIGN 2.4 names the first six; the others are the remaining
    rules of DESIGN 2.5 and 3.3 that can remove a term. One enumeration, here,
    so that the stores of the model give the record's own reasons (DESIGN 3.1:
    the model may not import the record)]"""

    NO_EPHEMERIS = "no_ephemeris"                 # 3.4 step 3, nothing stored
    OUTSIDE_VALIDITY = "outside_validity"         # 2.5 validity of the ephemeris
    IODE_MISMATCH = "iode_mismatch"               # E1
    EPHEMERIS_TYPE = "ephemeris_type"             # E2
    UNHEALTHY = "unhealthy"                       # E5
    STALE = "stale"                               # 2.5 ages
    OUTSIDE_AREA = "outside_area"                 # A4
    LATITUDE_LIMIT = "latitude_limit"             # A4, beyond ±85° (SPARTN 8.10)
    DNU = "dnu"                                   # C6
    SIOU_MISMATCH = "siou_mismatch"               # S3, C5, A1
    AIOU_MISMATCH = "aiou_mismatch"               # A1, A2
    NO_CORRECTIONS = "no_corrections"             # 3.4 step 1
    IODE_CHANGE = "iode_change"                   # C4, clock held
    NO_BIAS = "no_bias"                           # C5, C9
    NO_ATMOSPHERE = "no_atmosphere"               # A1, A6
    INVALID_GRID_NODE = "invalid_grid_node"       # A7
    GRID_MISMATCH = "grid_mismatch"               # A7, residual count not the GAD's node count
    NO_LEAP_SECONDS = "no_leap_seconds"           # T4
    BELOW_ELEVATION = "below_elevation"           # 2.5, Config.elevation_cutoff
    NOT_CONVERGED = "not_converged"               # 2.3, light-time iteration; 3.4, orbit kernel
    NO_FREQUENCY = "no_frequency"                 # GLONASS without a frequency channel



# ---- inputs: broadcast ephemerides, time, position -----------------------------------


@dataclass(frozen=True, slots=True)
class KeplerEphemeris:
    """A GPS, Galileo, BeiDou or QZSS broadcast ephemeris in physical units.

    ``message`` is the RTCM message it came from (1019, 1042, 1044, 1045,
    1046), which says the ephemeris type (rule E2). ``iod`` is the value a
    SPARTN orbit block names (rule E1): GPS and QZSS IODE, Galileo IODnav,
    BeiDou IODC = mod(toc/720, 240) from the BDT toc. ``toe`` and ``toc`` are
    on the GPS time scale, the week rollover resolved.

    Angles are in semicircles and angular rates in semicircles/s, as
    broadcast. [design: the kernel converts with the value of π its ICD
    defines (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2 gives 3.1415926535898), so the
    input holds no π] ``sqrt_a`` in m^½, harmonic corrections in m and rad, clock
    terms in s, s/s, s/s². ``group_delays`` are the broadcast group delays in
    seconds in message order (GPS TGD; Galileo BGD E5a/E1, BGD E5b/E1 as the
    message has them; BeiDou TGD1, TGD2), kept although rule E4 does not apply
    them. ``fit_interval`` is the GPS and QZSS fit interval flag; 0 elsewhere.
    ``iodc`` is the GPS and QZSS 10-bit IODC, which with a fit interval flag of
    1 gives the fit's length (IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-XII); None
    elsewhere.
    """

    sat: Sat
    message: int
    iod: int
    toe: GnssTime
    toc: GnssTime
    sqrt_a: float
    e: float
    m0: float
    delta_n: float
    omega: float
    omega0: float
    omega_dot: float
    i0: float
    idot: float
    cuc: float
    cus: float
    crc: float
    crs: float
    cic: float
    cis: float
    af0: float
    af1: float
    af2: float
    group_delays: tuple[float, ...]
    health: int
    accuracy: int
    fit_interval: int = 0
    iodc: int | None = None

    def __post_init__(self) -> None:
        """[spec: RTCM Tables 3.5-21, 3.5-110, 3.5-111, 3.5-113; GLONASS has its own type]"""
        if self.sat.system == System.GLONASS:
            raise ValueError("a GLONASS ephemeris is a GlonassEphemeris")


@dataclass(frozen=True, slots=True)
class GlonassEphemeris:
    """A GLONASS broadcast ephemeris (RTCM 1020, GLONASS ICD 4.4) in physical units.

    ``tb`` is the reference time on the GPS time scale; ``iod`` the value a
    SPARTN orbit block names (SF019, DESIGN Q3). Position, velocity and
    lunisolar acceleration are PZ-90 in m, m/s, m/s². ``tau`` is τn (s),
    ``gamma`` γn, ``delta_tau`` Δτn (s), ``p1`` the interval P1 in seconds.
    """

    sat: Sat
    frequency_channel: int
    iod: int
    tb: GnssTime
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    acceleration: tuple[float, float, float]
    tau: float
    gamma: float
    delta_tau: float
    health: int
    p1: float

    def __post_init__(self) -> None:
        """[spec: GLONASS ICD 3.3.1.1 (channels -7 to 6), 4.4]"""
        if self.sat.system != System.GLONASS:
            raise ValueError("a GlonassEphemeris is for a GLONASS satellite")
        if self.frequency_channel not in GLONASS_CHANNELS:
            raise ValueError(f"frequency channel {self.frequency_channel} outside -7 to 6")


@dataclass(frozen=True, slots=True)
class LeapSeconds:
    """GPS − UTC in whole seconds, and the next change when one is announced.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 (ΔtLS, ΔtLSF and
    the time of the change); u-blox F9 HPG 1.51 3.15.26 and u-blox X20 HPG 2.11
    3.14.22, NAV-TIMELS] ``effective`` is the GPS time from which ``next`` applies.
    """

    current: int
    next: int | None = None
    effective: GnssTime | None = None

    def __post_init__(self) -> None:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4, a change has a value and a time]"""
        if (self.next is None) != (self.effective is None):
            raise ValueError("an announced change needs both its value and its time")

    def at(self, t: GnssTime) -> int:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] GPS − UTC at GPS time t,
        the change applied from ``effective``."""
        if self.next is not None and self.effective is not None and t >= self.effective:
            return self.next
        return self.current


@dataclass(frozen=True, slots=True)
class UtcTime:
    """A UTC date and time of day (NMEA 0183 ZDA as gnsswire decodes it)."""

    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float


@dataclass(frozen=True, slots=True)
class RoverPosition:
    """A rover position (NMEA 0183 GGA): WGS-84 latitude and longitude in
    radians, ellipsoidal height in metres (GGA altitude plus geoid
    separation), and the UTC time of day in seconds the fix is for."""

    latitude: float
    longitude: float
    height: float
    utc_seconds: float


class FailureKind(StrEnum):
    """Why a frame did not become a message (DESIGN 3.5). [design: the cases of DESIGN 3.5]"""

    ENCRYPTED = "encrypted"
    DECODE = "decode"
    UNSUPPORTED = "unsupported"


class Protocol(StrEnum):
    """The wire protocol a frame was read in (DESIGN 2.2, framing helper).

    [spec: SPARTN Table 7.1 TF001; RTCM 4 Table 4-1; u-blox F9 HPG 1.51 and
    u-blox X20 HPG 2.11 3.2 (UBX) and 2.1 (NMEA)]
    """

    SPARTN = "spartn"
    RTCM3 = "rtcm3"
    UBX = "ubx"
    NMEA = "nmea"


@dataclass(frozen=True, slots=True)
class FrameFailure:
    """A frame that failed, reported so that it is counted (rule S7).

    ``protocol`` is the frame's; ``message_type`` and ``subtype`` as far as
    the frame could be read: SPARTN TF002 and TF007; an RTCM message number
    (DF002) with ``subtype`` None, or for message 4076 the IGS SSR sub-type
    (IDF002) when the frame is long enough; a UBX class and ID as class × 256 + ID
    with ``subtype`` None; nothing for NMEA. ``received`` is the receipt
    time when the caller has it. [design: DESIGN 3.5; the protocol keeps
    numbers of different protocols apart]
    """

    kind: FailureKind
    protocol: Protocol
    message_type: int | None = None
    subtype: int | None = None
    received: GnssTime | None = None

    def __post_init__(self) -> None:
        """[design: DESIGN 3.5, the kinds and protocols are closed sets]"""
        object.__setattr__(self, "kind", FailureKind(self.kind))
        object.__setattr__(self, "protocol", Protocol(self.protocol))


@dataclass(frozen=True, slots=True)
class Carried:
    """A message together with the transport it arrived on.

    The framing helper wraps the receiver messages whose meaning depends on
    what came before them on the same transport: the navigation words of
    UBX-RXM-SFRBX, which are joined into ephemerides per transport
    (``adapters.navdata``). [design: the navigation words of two receivers
    must not meet, and only the framer knows the transport (DESIGN 2.2)]
    """

    transport: str
    message: Any


# ---- the result of feeding one input --------------------------------------------------


class FeedOutcome(StrEnum):
    """What ``Engine.feed`` did with an input. [design: DESIGN 2.1, every rejection visible]"""

    STORED = "stored"            # stored; may replace an older item
    SUPERSEDED = "superseded"    # not stored: the store holds a newer item
    HELD = "held"                # kept aside until the time or leap seconds are known (T4, T5)
    DUPLICATE = "duplicate"      # equal to a stored message (S8)
    COUNTED = "counted"          # a frame failure, counted (S7)
    REJECTED = "rejected"        # not usable; ``reason`` says why


class FeedReason(StrEnum):
    """Why an input was not simply stored, by the rule that decided.
    [design: DESIGN 2.1, one reason per rule that can refuse or hold]"""

    UNSUPPORTED_SYSTEM = "unsupported_system"   # DESIGN 2.2
    NO_TIME = "no_time"                         # T5
    NO_LEAP_SECONDS = "no_leap_seconds"         # T4
    HOLD_OVERFLOW = "hold_overflow"             # T5, an older held message dropped
    STALE = "stale"                             # older than the issue-of-update lifetime (2.5)
    OLDER = "older"                             # S5, C4: the store holds a newer one
    CONFLICT = "conflict"                       # S8, same key and time tag, different body
    UNHEALTHY = "unhealthy"                     # E5, stored but not used
    DNU = "dnu"                                 # C6
    FRAME_FAILURE = "frame_failure"             # S7
    INCONSISTENT = "inconsistent"               # 2.2, fields that cannot form an input (a cutover)
    NOT_AN_INPUT = "not_an_input"               # 2.2, a message with a table that is not an input


@dataclass(frozen=True, slots=True)
class FeedResult:
    """The result of feeding one input: what happened, why, and to which satellites or areas.

    [design: DESIGN 2.1, ``feed`` never raises for bad data, so the result
    is where every effect is visible]
    """

    outcome: FeedOutcome
    reason: FeedReason | None = None
    sats: tuple[Sat, ...] = ()
    areas: tuple[int, ...] = ()
    detail: str = ""

    def __post_init__(self) -> None:
        """[design: a refusal or a hold always says why]"""
        object.__setattr__(self, "outcome", FeedOutcome(self.outcome))
        if self.reason is not None:
            object.__setattr__(self, "reason", FeedReason(self.reason))
        if self.outcome in (FeedOutcome.REJECTED, FeedOutcome.HELD, FeedOutcome.SUPERSEDED,
                            FeedOutcome.COUNTED) and self.reason is None:
            raise ValueError(f"outcome {self.outcome.value} needs a reason")
