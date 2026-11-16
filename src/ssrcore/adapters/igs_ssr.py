"""IGS SSR (RTCM proprietary message 4076) to the engine's correction inputs.

IGS State Space Representation Format v1.00 (October 2020), "IGS SSR"
below; RTCM 10403.3 3.5.13 where it says more. gnsswire decodes the frames
(``rtcm3.IgsOrbit`` and the others); this module turns each decoded message
into the input the engine already takes: an ``OcbInput`` of orbit, clock
and bias blocks for the corrections store, or a ``VtecInput`` for the
ionosphere provider of ``model.igs_vtec``. The rules, IG1 to IG12 of DESIGN
3.3, are stated where they are applied:

- IG1 stream identity (``IgsStreamId``), IG2 issue of update and end of set
  (``IgsSsr.input``), IG3 epoch and reference time (``epoch_time``,
  ``reference_time``), IG4 orbit and its pairing with a clock (``_orbit``,
  ``_pair``), IG5 clock (``_clock_polynomial``, ``_high_rate``,
  ``_clock_correction``), IG6 biases (``_biases``), IG7 signals
  (``signal``), IG8 issue of data (``ephemeris_iod``), IG9 URA
  (``ure_index``), IG10 yaw (``_yaw_at``), IG11 reference point
  (``IgsReference``), IG12 what is left out.

The adapter keeps what one message needs of earlier ones (the orbit update
interval, the latest clock polynomial, the other half of a bias block, the
yaw, the URA, the discontinuity counters), so it is an object the engine
holds, as it holds ``NavData``; it reads no clock and opens nothing.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from gnsswire import rtcm3

from ssrcore.model.ephemeris import EphemerisStore, own_scale
from ssrcore.model.igs_vtec import VtecInput, VtecLayer, quality
from ssrcore.model.types import (
    ClockCorrection,
    CodeBias,
    GlonassEphemeris,
    GnssTime,
    MaxAge,
    OcbInput,
    OcbSatellite,
    OrbitCorrection,
    PhaseBias,
    Sat,
    Signal,
    SpartnTag,
    StreamId,
    System,
    WEEK_SECONDS,
)
from ssrcore.toolkit import orbits, ssr, timescale

# [spec: IGS SSR v1.00 Table 5, the messages of sub-types 21 to 127 and 201]
IGS_SATELLITE_MESSAGES = (rtcm3.IgsOrbit, rtcm3.IgsClock, rtcm3.IgsCombined, rtcm3.IgsHighRateClock,
                          rtcm3.IgsCodeBias, rtcm3.IgsPhaseBias, rtcm3.IgsUra)
IGS_MESSAGES = IGS_SATELLITE_MESSAGES + (rtcm3.IgsVtec,)

# [spec: IGS SSR v1.00 Table 5, the system of each group of sub-types; design: SBAS
# is not a system of the engine (DESIGN 2.2), so its messages are refused]
_SYSTEMS = {"GPS": System.GPS, "GLONASS": System.GLONASS, "Galileo": System.GALILEO,
            "QZSS": System.QZSS, "BDS": System.BEIDOU}

# [spec: IGS SSR v1.00 Table 7 IDF013 to IDF022, IDF025, IDF028; resolutions in m, m/s, m/s²]
_RADIAL, _ALONG = 1e-4, 4e-4
_DOT_RADIAL, _DOT_ALONG = 1e-6, 4e-6
_C0, _C1, _C2 = 1e-4, 1e-6, 2e-8
_HIGH_RATE = 1e-4
_CODE_BIAS = 0.01
_PHASE_BIAS = 1e-4
# [spec: IGS SSR v1.00 Table 7 IDF026, 1/256 semicircle; IDF027, 1/8192 semicircle/s]
_YAW_DEG = 180.0 / 256.0
_YAW_RATE_DEG = 180.0 / 8192.0
# [spec: IGS SSR v1.00 Table 7 IDF039, IDF040, 0.005 TECU; −163.84 TECU "not available";
# IDF036, 10 km]
_TEC = 0.005
_TEC_NOT_AVAILABLE = -32768
_LAYER_HEIGHT = 10e3
# [spec: IGS SSR v1.00 Table 7 IDF004, the update intervals start at 00:00:00 of each day]
_DAY = 86400
# [spec: RTCM 10403.3 3.5.13 DF380, DF382: identifiers IGS SSR v1.00 IDF024 leaves
# reserved and RTCM defines, with the same numbering for every identifier both
# define: GPS 9 L2 L2C (M+L), Galileo 11 E5 I and 12 E5 Q; design: the IGS table
# is RTCM's extended, so where only RTCM names an identifier its name is taken]
_RTCM_SIGNALS = {(System.GPS, 9): "2X", (System.GALILEO, 11): "8I", (System.GALILEO, 12): "8Q"}
# [spec: IGS SSR v1.00 Table 7 IDF034, the special values: 000 000 unknown, 111 111 above 5466.5 mm]
_URA_UNKNOWN, _URA_ABOVE = 0, 63
# [spec: IGS SSR v1.00 Table 7 IDF031, a 4-bit counter, 0 to 15]
_COUNTER = 16


class ReferencePoint(StrEnum):
    """The satellite reference point of a stream's orbits: the antenna phase centre of
    the system's reference frequency, or the centre of mass.
    [spec: IGS SSR v1.00 4.1, Table 2]"""

    APC = "apc"
    COM = "com"


@dataclass(frozen=True, slots=True)
class IgsStreamId(StreamId):
    """IG1. The stream of an IGS SSR message: ``solution`` the SSR Solution ID
    (IDF009), ``processor`` the SSR Provider ID (IDF008). Never equal to a SPARTN
    ``StreamId`` with the same numbers, so the two formats' streams are apart
    (rule S1) and every rule keyed by stream (S2 to S10, A10, A11) applies as is.
    [spec: IGS SSR v1.00 Table 7 IDF008, "to identify a SSR service", "globally
    unique"; IDF009, "different SSR services of one SSR provider"; design: the
    pair plays the part of SPARTN's TF010, TF011; a subclass keeps the engine's
    one stream type]"""

    @classmethod
    def of(cls, message: Any) -> IgsStreamId:
        """[spec: IGS SSR v1.00 Table 7 IDF008, IDF009] the stream of a decoded message."""
        return cls(message.solution, message.provider)


@dataclass(frozen=True, slots=True)
class IgsReference:
    """IG11. The caller's statement of a stream's satellite reference point, fed to
    the engine like any input. Until it is fed, or when it says centre of mass,
    the stream's satellite messages are refused (``not_an_input``); its VTEC
    (IM201) is taken either way, being free of the reference point.
    [spec: IGS SSR v1.00 4.1, "Services using the IGS SSR format can refer to
    either the satellite CoM or the average APC of a GNSS specific reference
    frequency"; no field of the messages says which ("the actual CRS will be
    identified by the configured stream of the service provider"), and a CoM
    service needs "the full satellite attitude ... as well as the CoM/APC
    offsets", messages "under development"; design: the record holds an
    antenna phase centre (DESIGN 2.4), so only a stream declared APC is used,
    and the declaration is an input so a replay carries it]"""

    provider: int
    solution: int
    point: ReferencePoint

    def __post_init__(self) -> None:
        """[spec: IGS SSR v1.00 Table 7 IDF008 0 to 65535, IDF009 0 to 15]"""
        object.__setattr__(self, "point", ReferencePoint(self.point))
        if not 0 <= self.provider <= 65535 or not 0 <= self.solution <= 15:
            raise ValueError("provider 0 to 65535 (IDF008), solution 0 to 15 (IDF009)")

    @property
    def stream(self) -> IgsStreamId:
        """[spec: IGS SSR v1.00 IDF008, IDF009] the stream the statement is for."""
        return IgsStreamId(self.solution, self.provider)


@dataclass(frozen=True, slots=True)
class Adapted:
    """What one message gives the engine: the input, and whether it completes a
    block the engine already has at the same epoch (a bias block's other half, a
    high-rate clock on its polynomial), so that its replacement there is no conflict.
    [design: DESIGN 2.1, every effect visible in the feed result]"""

    item: OcbInput | VtecInput
    completes: bool = False
    detail: str = ""


@dataclass(frozen=True, slots=True)
class Kept:
    """A message that changes only the adapter's own state (a URA, kept for the clock
    corrections that follow), or repeats the previous message of its sub-type and
    stream exactly. [design: DESIGN 2.1; S8, a repeat is a duplicate]"""

    detail: str
    duplicate: bool = False


@dataclass(frozen=True, slots=True)
class Refused:
    """A message the adapter cannot turn into an input, and why.
    [design: DESIGN 2.1, the engine turns it into a rejection with its reason]"""

    reason: str
    unsupported_system: bool = False


# ---- rules on one value ------------------------------------------------------------------


def system_of(message: Any) -> System | None:
    """[spec: IGS SSR v1.00 Table 5, sub-types 21 to 47 GPS and GLONASS, 61 to 107
    Galileo, QZSS and BDS, 121 to 127 SBAS] the engine's system of a satellite
    message, None for SBAS."""
    name, _ = rtcm3.IGS_SUBTYPES[message.subtype]
    return _SYSTEMS.get(name)


def satellite(system: System, satellite_id: int) -> Sat | None:
    """IG1 for satellites. [spec: IGS SSR v1.00 Table 7 IDF011: GPS, GLONASS and BDS
    1 to 63 and 0 = 64; Galileo 1 to 36; QZSS 1 to 10 = PRN 193 to 202, which are
    the engine's QZSS numbers 1 to 10 (``Sat``, SPARTN SF095)] the satellite, or None
    for an ID the document reserves or the engine cannot name (GLONASS 64)."""
    name = next(k for k, v in _SYSTEMS.items() if v == system)
    number = rtcm3.IGS_SATELLITES[name].get(satellite_id)
    if number is None:
        return None
    if system == System.QZSS:
        number -= 192
    try:
        return Sat(system, number)
    except ValueError:
        return None


def signal(system: System, signal_id: int) -> Signal | None:
    """IG7. [spec: IGS SSR v1.00 Table 7 IDF024, the GNSS Signal and Tracking Mode
    Identifier, RINEX 3 band and attribute per system; RTCM 10403.3 DF380, DF382 for
    three identifiers IGS reserves (``_RTCM_SIGNALS``)] the engine's signal, or None
    for an identifier neither document defines or a band the engine has no
    frequency for (GLONASS G3, G1a, G2a; QZSS L6)."""
    name = next(k for k, v in _SYSTEMS.items() if v == system)
    code = rtcm3.IGS_SIGNALS[name].get(signal_id) or _RTCM_SIGNALS.get((system, signal_id))
    if code is None:
        return None
    try:
        return Signal(system, code)
    except ValueError:
        return None


def _place(week: int, seconds: float, near: GnssTime) -> GnssTime:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV, a seconds-of-week
    time is within half a week of the time it is used at; design: as ``adapters.rtcm3`` places an
    ephemeris, by the time near receipt]"""
    t = GnssTime.normalised(week, seconds)
    dt = t - near
    if dt > WEEK_SECONDS / 2:
        t = t - WEEK_SECONDS
    elif dt < -WEEK_SECONDS / 2:
        t = t + WEEK_SECONDS
    return t


def epoch_time(message: Any, system: System | None, near: GnssTime) -> GnssTime:
    """IG3. The SSR epoch time IDF003 as a GPS time, its week the one that puts it
    within half a week of ``near``.

    IDF003 counts "full seconds since the beginning of the week of continuous
    time scale with no offset from GPS, Galileo, QZSS, SBAS, UTC leap seconds
    from GLONASS, -14 s offset from BDS": GPS time for GPS, Galileo, QZSS and
    GLONASS, and for BDS the scale 14 s behind it, BDT, which is converted.
    [spec: IGS SSR v1.00 Table 7 IDF003, IDF004 ("the update intervals are
    aligned to the GPS time scale for all GNSS"); BDS ICD 3.3, BDT = GPS time −
    14 s, through ``toolkit.timescale``; design: the field note can be read as
    GPS time for every system or as BDT for BDS; the recorded IGS stream sends
    BDS epochs 14 s below its GPS and Galileo epochs, on the 5 s grid of its
    update interval only as BDT, so BDT is the reading taken]"""
    if system == System.BEIDOU:
        bdt_week, _ = timescale.gps_to_bdt(near.week, near.seconds)
        week, sow = timescale.bdt_to_gps(bdt_week, float(message.epoch))
        return _place(week, sow, near)
    return _place(near.week, float(message.epoch), near)


def update_interval(message: Any) -> float:
    """[spec: IGS SSR v1.00 Table 7 IDF004] the update interval in seconds."""
    return float(rtcm3.update_interval(message.interval))


def reference_time(message: Any, epoch: GnssTime) -> GnssTime:
    """IG3. t0 = epoch + ½ update interval; the epoch itself for interval "0".
    [spec: IGS SSR v1.00 4, "The reference time for the correction terms is
    computed from the SSR Epoch Time 1 s data field (IDF003) plus half the SSR
    Update Interval. Exception is SSR Update Interval "0", which uses the SSR
    Epoch Time as reference time"; IDF013, IDF019]"""
    if message.interval == 0:
        return epoch
    return epoch + update_interval(message) / 2.0


def tag_of(t: GnssTime) -> SpartnTag:
    """[design: the engine's inputs carry a SPARTN tag (DESIGN 2.2); an IGS epoch,
    whole GPS seconds, is the 32-bit tag of the GPS scale, which the time base
    turns back into the same GPS time (T1)]"""
    return SpartnTag(timescale.spartn_tag_from_gps(t.week, t.seconds, System.GPS), True)


def ephemeris_type(sat: Sat) -> int:
    """IG8. The ephemeris an IGS orbit corrects, as the SPARTN ephemeris type the
    store looks up (rule E2): Galileo I/NAV (1, message 1046), BeiDou D2 for GEO
    (1) and D1 otherwise (0), the one message of every other system (0).
    [spec: IGS SSR v1.00 4.1 Table 1 (GPS NAV, GLONASS M, Galileo I/NAV, BDS
    D1/D2, QZSS NAV); BDS ICD 5.1.1, D2 from GEO satellites; SPARTN Table 6.2
    SF096, SF097 through ``ephemeris_message``]"""
    if sat.system == System.GALILEO:
        return 1
    if sat.system == System.BEIDOU and orbits.beidou_geo(sat.number):
        return 1
    return 0


def ephemeris_iod(sat: Sat, iod: int, ephemerides: EphemerisStore) -> int | None:
    """IG8. The issue of data the ephemeris store names (rule E1) for an IDF012 value,
    or None when the value names none.

    GPS and QZSS send the IODE and GLONASS, with bit 7 clear, tb: the store's
    values. Galileo sends the 8 LSBs of IODnav and BDS mod(toe/720, 240): the
    stored issue that matches is found (Galileo among the I/NAV issues), and
    without one the value is taken as sent, which is the IODnav itself below
    256. [spec: IGS SSR v1.00 Table 7 IDF012; DESIGN 2.2, the store keeps
    Galileo's IODnav and BeiDou's mod(toc/720, 240) (SPARTN SF100); design: an
    issue received after the correction is matched by the next correction, one
    update interval later; GLONASS with bit 7 set has no definition]"""
    if sat.system == System.GLONASS:
        return None if iod & 0x80 else iod & 0x7F
    if sat.system == System.GALILEO:
        for stored in ephemerides.issues(sat):
            e = stored.ephemeris
            if stored.message == 1046 and e.iod & 0xFF == iod:
                return e.iod
        return iod
    if sat.system == System.BEIDOU:
        for stored in ephemerides.issues(sat):
            e = stored.ephemeris
            if isinstance(e, GlonassEphemeris):
                continue
            toe = own_scale(System.BEIDOU, e.toe)[1]
            if int(toe) // 720 % 240 == iod:
                return e.iod
        return iod
    return iod


def ure_index(ura: int) -> int:
    """IG9. The SPARTN SF024 index the record carries for an IDF034 URA: the
    smallest whose bound is not below the URA's, 7 (above 1 m) beyond the last,
    0 for unknown. [spec: IGS SSR v1.00 Table 7 IDF034, URA [mm] ≤
    3^CLASS (1 + VALUE/4) − 1, 000 000 unknown, 111 111 above 5466.5 mm; SPARTN
    Table 6.2 SF024; design: the record's URE field is SPARTN's index (DESIGN
    2.5), and the coarser of the two scales keeps the bound true]"""
    if ura == _URA_UNKNOWN:
        return 0
    if ura == _URA_ABOVE:
        return 7
    bound = (3.0 ** (ura >> 3) * (1.0 + (ura & 7) / 4.0) - 1.0) / 1000.0
    for index in range(1, 7):
        step = ssr.field(24, index)
        if step.status == ssr.Status.OK and step.value >= bound - 1e-12:
            return index
    return 7


# ---- the adapter -------------------------------------------------------------------------


@dataclass(slots=True)
class _Part:
    """A message's satellite part kept for the messages that follow."""

    epoch: GnssTime
    interval: float
    iod_ssr: int
    value: Any


class IgsSsr:
    """The adapter of one engine: IGS SSR messages in, typed inputs out (``input``).

    ``max_age`` is the engine's (``Config.max_age``): a kept part (the other
    half of a bias block, a yaw, a URA, a clock polynomial under a high-rate
    clock) is combined with a later message only while its age is within the
    limit of its kind or its own update interval, whichever is longer.
    [spec: IGS SSR v1.00 4, "The validity interval of SSR parameters is at
    least the SSR Update Interval"; design: DESIGN 2.5 age limits]
    """

    def __init__(self, max_age: MaxAge | None = None) -> None:
        """[design: DESIGN 2.1, state only from the inputs; counts for the status]"""
        self.max_age = MaxAge() if max_age is None else max_age
        self.references: dict[IgsStreamId, ReferencePoint] = {}
        self.left_out: Counter[str] = Counter()
        self._orbit_interval: dict[tuple[IgsStreamId, System], float] = {}
        self._datum: dict[IgsStreamId, int] = {}
        self._clock: dict[tuple[IgsStreamId, Sat], _Part] = {}
        self._code: dict[tuple[IgsStreamId, Sat], _Part] = {}
        self._phase: dict[tuple[IgsStreamId, Sat], _Part] = {}
        self._yaw: dict[tuple[IgsStreamId, Sat], _Part] = {}
        self._ura: dict[tuple[IgsStreamId, Sat], _Part] = {}
        self._counters: dict[tuple[IgsStreamId, Sat, Signal], tuple[int, GnssTime]] = {}
        self._last: dict[tuple[int, IgsStreamId], tuple[Any, Adapted]] = {}
        self._issued: dict[tuple[IgsStreamId, Sat], int] = {}
        self._held: dict[tuple[IgsStreamId, Sat], _Part] = {}

    def declare(self, reference: IgsReference) -> None:
        """IG11. [design: see ``IgsReference``] the caller's statement, latest kept."""
        self.references[reference.stream] = reference.point

    def refusal(self, message: Any) -> Refused | None:
        """IG11, IG12. Why a satellite message cannot be used, None when it can.
        [spec: IGS SSR v1.00 4.1, Table 5; design: see ``IgsReference`` and ``system_of``]"""
        if system_of(message) is None:
            return Refused("SBAS is not a system of the engine (DESIGN 2.2)", unsupported_system=True)
        stream = IgsStreamId.of(message)
        point = self.references.get(stream)
        if point is None:
            return Refused(f"no reference point declared for {stream} (IgsReference, IGS SSR 4.1)")
        if point == ReferencePoint.COM:
            return Refused(f"{stream} gives centre-of-mass orbits; the record holds an antenna phase "
                           "centre and IGS SSR 4.1 has no message for the CoM-APC offset yet")
        return None

    def input(self, message: Any, near: GnssTime, ephemerides: EphemerisStore) -> Adapted | Kept | Refused:
        """IG2. One decoded 4076 message as an input of the engine, ``near`` a GPS time
        near its epoch (receipt, or the time base's current time).

        The OCB's issue of update is the IOD SSR, so a change of the generating
        configuration starts a new issue (rule S3); end of set is never claimed,
        since IDF005 = 0 ends one sub-type's sequence and not every message of
        the epoch, so a set ends by a newer epoch (rule S5). A message equal to
        the previous one of its sub-type and stream is a duplicate (S8) and
        changes nothing, not even a block another message has completed since.
        [spec: IGS SSR v1.00 Table 7 IDF005 "0 – last message of a sequence with
        the same IGS Message Number and SSR Epoch Time", IDF007; RTCM 3.5.13.11;
        design: DESIGN 3.3 S3, S5, S8]"""
        stream = IgsStreamId.of(message)
        key = (message.subtype, stream)
        last = self._last.get(key)
        if last is not None and last[0] == message:
            return Kept(f"repeat of the previous IGS sub-type {message.subtype} message", duplicate=True)
        if isinstance(message, rtcm3.IgsVtec):
            adapted = self._vtec(message, stream, near)
        else:
            refused = self.refusal(message)
            if refused is not None:
                return refused
            adapted = self._satellites(message, stream, near, ephemerides)
        if not isinstance(adapted, Refused):
            self._last[key] = (message, adapted)
        return adapted

    # ---- satellite messages ----------------------------------------------------------

    def _satellites(self, message: Any, stream: IgsStreamId, near: GnssTime,
                    ephemerides: EphemerisStore) -> Adapted | Kept | Refused:
        """[spec: IGS SSR v1.00 8.1 to 8.3] one OCB of the message's satellites."""
        system = system_of(message)
        epoch = epoch_time(message, system, near)
        interval = update_interval(message)
        if isinstance(message, (rtcm3.IgsOrbit, rtcm3.IgsCombined)):
            # [spec: IGS SSR v1.00 Table 7 IDF004, consistency of orbit and clock through
            # the update interval; IDF006, the CRS of the orbit messages]
            self._orbit_interval[(stream, system)] = interval
            self._datum[stream] = message.regional
        blocks: list[OcbSatellite] = []
        completes = False
        yaw_present = held = False
        for s in message.satellites:
            sat = satellite(system, s.satellite)
            if sat is None:
                self.left_out["satellite id"] += 1
                continue
            orbit = clock = phase = code = None
            if isinstance(message, (rtcm3.IgsOrbit, rtcm3.IgsCombined)):
                orbit = self._orbit(s, sat, message, stream, epoch, ephemerides)
            if isinstance(message, (rtcm3.IgsClock, rtcm3.IgsCombined)):
                clock = self._clock_polynomial(s, sat, message, stream, epoch, interval)
            if isinstance(message, (rtcm3.IgsOrbit, rtcm3.IgsClock)):
                orbit = self._pair(stream, sat, message.iod_ssr, epoch, interval, orbit, clock)
                held |= orbit is None and (stream, sat) in self._held
            elif orbit is not None:
                self._issued[(stream, sat)] = orbit.iode
                self._held.pop((stream, sat), None)
            yaw_present |= orbit is not None and orbit.yaw is not None
            if isinstance(message, rtcm3.IgsHighRateClock):
                clock = self._high_rate(s, sat, message, stream, epoch)
                completes |= clock is not None
            if isinstance(message, (rtcm3.IgsCodeBias, rtcm3.IgsPhaseBias)):
                phase, code = self._biases(s, sat, message, stream, epoch, interval)
                completes = True
            if isinstance(message, rtcm3.IgsUra):
                self._ura[(stream, sat)] = _Part(epoch, interval, message.iod_ssr, ure_index(s.ura))
                continue
            if orbit is None and clock is None and phase is None and code is None:
                continue
            blocks.append(OcbSatellite(sat, False, None, orbit, clock, phase, code))
        if isinstance(message, rtcm3.IgsUra):
            return Kept("URA kept for the clock corrections that follow (IG9)")
        if not blocks and held:
            return Kept("orbits of a new IOD held for the clock of their epoch (IG4)")
        if not blocks:
            return Refused("no satellite of the message can be used (IG12)")
        # [design: the OCB's system is the scale of its tag (TF009, T1), and an IGS
        # epoch is on the GPS scale once converted (IG3); the satellites name theirs]
        ocb = OcbInput(System.GPS, stream, tag_of(epoch), message.iod_ssr, False, yaw_present,
                       self._datum.get(stream, 0), ephemeris_type(blocks[0].sat), tuple(blocks))
        return Adapted(ocb, completes)

    def _fresh(self, part: _Part | None, epoch: GnssTime, iod_ssr: int, limit: float) -> bool:
        """[spec: IGS SSR v1.00 4, validity at least the update interval; IDF007, parts
        of one configuration only; design: DESIGN 2.5 age limits]"""
        return (part is not None and part.iod_ssr == iod_ssr
                and 0.0 <= epoch - part.epoch <= max(limit, part.interval))

    def _orbit(self, s: Any, sat: Sat, message: Any, stream: IgsStreamId, epoch: GnssTime,
               ephemerides: EphemerisStore) -> OrbitCorrection | None:
        """IG4. δO and δȮ in m and m/s at t0, signed as the engine's (X = X_b − δX), with
        the ephemeris the IOD names and the yaw of IG10. A BeiDou GEO orbit is left out.
        [spec: IGS SSR v1.00 4.1, X_orbit = X_broadcast − δX, δO(t) = δO + δȮ (t − t0),
        the frame of SPARTN 8.8 (C1); Table 7 IDF013 to IDF018. design: IGS SSR 4.1
        builds the frame of a GEO satellite from ṙ_b + Ω̇_e × X_b, which the engine's
        C1 (ṙ_b) does not, so such a correction would be applied in another frame]"""
        if sat.system == System.BEIDOU and orbits.beidou_geo(sat.number):
            self.left_out["BeiDou GEO orbit"] += 1
            return None
        iod = ephemeris_iod(sat, s.iod, ephemerides)
        if iod is None:
            self.left_out["GLONASS IOD with bit 7 set"] += 1
            return None
        t0 = reference_time(message, epoch)
        return OrbitCorrection(
            iod, s.radial * _RADIAL, s.along * _ALONG, s.cross * _ALONG,
            yaw=self._yaw_at(stream, sat, message.iod_ssr, t0),
            ephemeris_type=ephemeris_type(sat),
            rate=(s.dot_radial * _DOT_RADIAL, s.dot_along * _DOT_ALONG, s.dot_cross * _DOT_ALONG),
            reference=t0)

    def _pair(self, stream: IgsStreamId, sat: Sat, iod_ssr: int, epoch: GnssTime, interval: float,
              orbit: OrbitCorrection | None, clock: ClockCorrection | None) -> OrbitCorrection | None:
        """IG4. The orbit to send in this message's block. An orbit of the IOD last sent
        goes at once; an orbit of a new IOD (or the first) is held and sent in the
        block of the first clock of its epoch or later, the newest held replacing
        an older one, and dropped if no such clock comes within the orbit's age.
        [spec: IGS SSR v1.00 4.2, the clock correction is "w.r.t. GNSS broadcast
        ephemeris" of the IOD the orbit message names; RTCM 10403.3 3.5.13.2, the
        orbit and clock messages of one epoch together give the correction; 4, "It is
        not allowed to send corrections for more than one IOD for the same orbit
        Update Interval". design: the store pairs the newest orbit with the newest
        clock and rule C4 rejects only a clock later than its orbit's IOD, so an
        orbit of a new IOD sent alone would meet the previous epoch's clock, made
        against the old IOD's broadcast clock, and give a wrong clock (on the IGS
        recording, up to 0.22 m for Galileo, until the next clock); held, the old
        orbit and clock stay paired until the new ones are both in. The orbit then
        takes the clock's tag, its age counted from that epoch, its t0 unchanged]"""
        key = (stream, sat)
        if orbit is not None:
            if self._issued.get(key) == orbit.iode:
                self._held.pop(key, None)
                return orbit
            self._held[key] = _Part(epoch, interval, iod_ssr, orbit)
            return None
        held = self._held.get(key)
        if clock is None or held is None:
            return None
        if held.iod_ssr != iod_ssr or epoch < held.epoch:
            return None
        del self._held[key]
        if epoch - held.epoch > max(self.max_age.orbit, held.interval):
            self.left_out["orbit of a new IOD without a clock"] += 1
            return None
        self._issued[key] = held.value.iode
        return held.value

    def _yaw_at(self, stream: IgsStreamId, sat: Sat, iod_ssr: int, t: GnssTime) -> float | None:
        """IG10. The yaw of the latest phase bias message, in degrees in [0, 360), at
        the orbit's t0, or None without one. [spec: IGS SSR v1.00 4.3, "positive
        rotation around the nadir-pointing z-axis" from the yaw origin, the frame of
        SPARTN 8.14 (C8); Table 7 IDF026, IDF027; 4, a correction term's reference
        time is epoch + ½ interval. design: the orbit block carries one yaw (C8), so
        the yaw and its rate are evaluated once, at the orbit's t0, the middle of the
        orbit's interval; a yaw is used while the phase biases it came with are]"""
        part = self._yaw.get((stream, sat))
        if not self._fresh(part, t, iod_ssr, self.max_age.bias):
            return None
        yaw, rate, t_yaw = part.value
        return (yaw + rate * (t - t_yaw)) % 360.0

    def _clock_polynomial(self, s: Any, sat: Sat, message: Any, stream: IgsStreamId, epoch: GnssTime,
                          interval: float) -> ClockCorrection:
        """IG5. C0, C1, C2 in m, m/s, m/s² at t0, signed as the engine's: negated.

        IGS SSR 4.2 writes t_satellite = t_broadcast − δC/c of satellite *times*,
        the system time "computed ... from broadcast clock parameters", which is
        t = t_sv − Δt_sv of the ICDs; so the clock offset is Δt_sv = Δt_b + δC/c.
        The engine's clock is SPARTN 8.7's, the same formula written of the
        offset (t_s = t_b − δC/c, a satellite clock c·Δt_b − δC in metres),
        so the engine's δC is the IGS −δC. The IGS recording shows it: at every
        one of its 33 broadcast issue changes the step of C0 is the negative of
        the step of c·Δt_b (within 1 cm), so the corrected clock is continuous
        only so read. [spec: IGS SSR v1.00 4.2, Table 7 IDF019 to IDF021;
        IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, t = t_sv − Δt_sv; SPARTN v2.0.3 8.7;
        design: the reading the recorded stream confirms, stated in DESIGN (IG5)]"""
        t0 = reference_time(message, epoch)
        polynomial = (-s.c0 * _C0, -s.c1 * _C1, -s.c2 * _C2, t0)
        self._clock[(stream, sat)] = _Part(epoch, interval, message.iod_ssr, polynomial)
        return self._clock_correction(stream, sat, message.iod_ssr, epoch, polynomial, 0.0)

    def _high_rate(self, s: Any, sat: Sat, message: Any, stream: IgsStreamId,
                   epoch: GnssTime) -> ClockCorrection | None:
        """IG5. A high-rate clock added to the polynomial of its interval: the clock
        message of the same stream, satellite and IOD SSR whose update interval
        contains the epoch. Without one it is left out.
        [spec: IGS SSR v1.00 4.2, "The High Rate Clock correction is added to the
        corresponding clock correction", δC = C0 + C1 (t − t0) + C2 (t − t0)² +
        δC_HighRate; Table 7 IDF022; design: the corresponding clock correction is
        the one whose interval holds the epoch (IG2's consistency by interval);
        negated with it, see ``_clock_polynomial``]"""
        part = self._clock.get((stream, sat))
        if part is None or part.iod_ssr != message.iod_ssr \
                or not 0.0 <= epoch - part.epoch < max(part.interval, 1.0):
            self.left_out["high-rate clock without its polynomial"] += 1
            return None
        return self._clock_correction(stream, sat, message.iod_ssr, epoch, part.value,
                                      -s.high_rate_clock * _HIGH_RATE)

    def _clock_correction(self, stream: IgsStreamId, sat: Sat, iod_ssr: int, epoch: GnssTime,
                          polynomial: tuple, high_rate: float) -> ClockCorrection:
        """IG5. The clock block: its IODE continuity is the time since the start of the
        orbit update interval holding the epoch, so that rule C4 pairs the clock only
        with an orbit of that interval or later; the URE index of IG9.
        [spec: IGS SSR v1.00 4, "It is not allowed to send corrections for more than
        one IOD for the same orbit Update Interval. The consistency of SSR Orbit and
        Clock messages is maintained through the SSR Update Interval"; Table 7 IDF004,
        the intervals start at 00:00:00 of the day, GPS time for every system; design:
        SPARTN SF022 is the same statement, "the same IODE has been used over the last
        X seconds", so the rule of C4 needs no change]"""
        c0, c1, c2, t0 = polynomial
        u = self._orbit_interval.get((stream, sat.system))
        if u is None:
            part = self._clock.get((stream, sat))
            u = part.interval if part is not None else 1.0
        continuity = float((epoch.second % _DAY) % int(u)) + epoch.fraction
        ura = self._ura.get((stream, sat))
        ure = ura.value if self._fresh(ura, epoch, iod_ssr, self.max_age.clock) else 0
        return ClockCorrection(continuity, c0 + high_rate, ure, c1, c2, t0)

    def _biases(self, s: Any, sat: Sat, message: Any, stream: IgsStreamId, epoch: GnssTime,
                interval: float) -> tuple[tuple[PhaseBias, ...] | None, tuple[CodeBias, ...] | None]:
        """IG6. A satellite's bias block, code and phase from the two messages.

        Values are signed as the engine's: IGS biases "should be added to the PR
        or PHR measurements ... to get corrected ranges", so the value a
        measurement contains (b_s, B_s, DESIGN 2.4) is the negative. The fix flag
        is IDF029. A phase bias's continuity is the time since its discontinuity
        counter IDF031 last changed (or since first seen), the meaning of SPARTN's
        SF015, so a change ends the signal's arc (C7). The block is both halves
        (C10 replaces a block wholly): a message carries the other kind's latest
        biases of the same IOD SSR within their age; a code message leaves out the
        phase biases when a newer phase message is due at its epoch, so that the
        phase block of that epoch is the first at its time and a counter change in
        it can end the arc (C7 counts a discontinuity only at a later time).
        [spec: IGS SSR v1.00 4.3, the sign, "The coherent use of code and phase
        signals relies on the availability of information for both phase and code
        biases"; Table 7 IDF025, IDF028, IDF029, IDF031; design: as stated]"""
        key = (stream, sat)
        if isinstance(message, rtcm3.IgsCodeBias):
            code = []
            for b in s.biases:
                sg = signal(sat.system, b.signal)
                if sg is None:
                    self.left_out["code bias signal"] += 1
                    continue
                code.append(CodeBias(sg, -b.bias * _CODE_BIAS))
            code = tuple(code)
            self._code[key] = _Part(epoch, interval, message.iod_ssr, code)
            phase_part = self._phase.get(key)
            phase = None
            if self._fresh(phase_part, epoch, message.iod_ssr, self.max_age.bias) \
                    and epoch < phase_part.epoch + phase_part.interval:
                phase = tuple(PhaseBias(sg, value, fixed, epoch - since)
                              for sg, value, fixed, since in phase_part.value)
            return phase, code
        phase = []
        for b in s.biases:
            sg = signal(sat.system, b.signal)
            if sg is None:
                self.left_out["phase bias signal"] += 1
                continue
            counter = self._counters.get((stream, sat, sg))
            if counter is None or counter[0] != b.discontinuity % _COUNTER:
                counter = (b.discontinuity % _COUNTER, epoch)
                self._counters[(stream, sat, sg)] = counter
            phase.append((sg, -b.bias * _PHASE_BIAS, bool(b.integer), counter[1]))
        self._phase[key] = _Part(epoch, interval, message.iod_ssr, tuple(phase))
        phase = tuple(PhaseBias(sg, value, fixed, epoch - since) for sg, value, fixed, since in phase)
        # [spec: IGS SSR v1.00 Table 7 IDF026, IDF027; 4, t0 = epoch + ½ interval]
        self._yaw[key] = _Part(epoch, interval, message.iod_ssr,
                               (s.yaw * _YAW_DEG, s.yaw_rate * _YAW_RATE_DEG, reference_time(message, epoch)))
        code_part = self._code.get(key)
        code = code_part.value if self._fresh(code_part, epoch, message.iod_ssr, self.max_age.bias) else None
        return phase, code

    # ---- IM201 -------------------------------------------------------------------------

    def _vtec(self, message: rtcm3.IgsVtec, stream: IgsStreamId, near: GnssTime) -> Adapted:
        """[spec: IGS SSR v1.00 8.4.1 Tables 24 to 27, IDF035 to IDF041; IDF003 in GPS
        time, IM201 belonging to no other system] the message for ``model.igs_vtec``."""
        epoch = epoch_time(message, None, near)
        layers = []
        for layer in message.layers:
            cosine = tuple(None if c == _TEC_NOT_AVAILABLE else c * _TEC for c in layer.cosine)
            sine = tuple(None if c == _TEC_NOT_AVAILABLE else c * _TEC for c in layer.sine)
            layers.append(VtecLayer(layer.height * _LAYER_HEIGHT, layer.degree_minus_1 + 1,
                                    layer.order_minus_1 + 1, cosine, sine))
        return Adapted(VtecInput(stream, epoch, update_interval(message), message.iod_ssr,
                                 quality(message.quality), tuple(layers)))


def is_igs(message: Any) -> bool:
    """[spec: IGS SSR v1.00 2, every IGS SSR message is RTCM message 4076] whether a
    decoded message is one of IGS SSR's."""
    return isinstance(message, IGS_MESSAGES)

