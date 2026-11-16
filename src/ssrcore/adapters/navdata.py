"""Navigation words from a receiver (UBX-RXM-SFRBX) to ephemeris inputs, through gnsswire.nav.

DESIGN.md section 2.2 and question Q15. A receiver that decodes the
broadcast navigation message hands it on unit by unit: a GPS subframe, a
Galileo I/NAV page, a BeiDou subframe, a GLONASS string. ``gnsswire.nav``
reads each unit and checks its code (parity, CRC, BCH or Hamming), and its
``Assembler`` yields an ephemeris only when every piece of one issue of data
is present. That ephemeris is named and scaled as the RTCM 3 ephemeris
message of its system, and is handed on *as that RTCM message*
(``rtcm3.Msg1019``, ``Msg1020``, ``Msg1042``, ``Msg1046``): the RTCM
adapter (``adapters.rtcm3``) then scales it, places its reference times and
the engine stores it with the same code as an ephemeris that arrived over
RTCM, so the health, validity and consistency rules E1 to E5, the holds of
T4 and T5 and every scaling are one code for both sources.

``NavData`` holds the assembly state, per named transport, so that the words
of two receivers never meet. Its state is a function of the units pushed,
in order, and of the times passed with them; nothing reads a clock. It is
bounded: per transport, one assembly per satellite of each system, and
the GLONASS almanac's channel per slot.

What it yields for each unit pushed, in order:

- an RTCM ephemeris message, when the unit completes a new issue;
- a ``LeapSeconds``, when the unit carries the GPS − UTC count and it
  differs from the last one this transport gave (GPS subframe 4 page 18,
  Galileo word 6, BeiDou subframe 5 page 10 or D2 page 102);
- a ``FrameFailure``, when the unit fails its code or does not fit its table;
- a rejected ``FeedResult``, when the unit cannot be used (a signal whose
  words no staged document lays out; an ephemeris whose satellite number
  contradicts the carrier's; a GLONASS ephemeris without a frequency channel).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from gnsswire import nav, rtcm3, ubx
from gnsswire.nav import beidou, carriers, galileo, glonass, gps

from ssrcore.model.types import (
    FailureKind,
    FeedOutcome,
    FeedReason,
    FeedResult,
    FrameFailure,
    GnssTime,
    LeapSeconds,
    Protocol,
)
from ssrcore.toolkit import timescale

__all__ = ["Carrier", "GLONASS_ROW", "NavData", "SFRBX", "STALE_PARTS"]

_DAY = timescale.DAY_SECONDS

# [design: every part of an issue is broadcast again within 30 s (a GPS LNAV
#  frame, IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.2; a Galileo I/NAV nominal sub-frame,
#  Galileo OS SIS ICD 2.2 4.3.3 Table 40; a BeiDou D1 frame and the ten D2 pages of subframe
#  1, BDS ICD 5.2.2, 5.3.1; a GLONASS frame, ICD 4.3.2), so an assembly that holds a
#  part older than ten repetitions is started over. That costs at most one
#  repetition; keeping the part could
#  let a part of an old issue meet a new issue that reuses its issue of data,
#  which IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.4.4 allows after six hours and the Galileo
#  OS SIS ICD 2.2 5.1.9.2 does not bound]
STALE_PARTS = 300.0

# [spec: GLONASS ICD 4.3.2, a string lasts 2 s and strings 1 to 5 open each 30 s
#  frame, so the five are received within 10 s. design: a longer gap between two
#  strings of one stream means strings were lost unseen (a receiver passes on
#  only the strings it decoded), and strings before and after the gap may
#  belong to different frames, so the row is broken. The gap is measured in
#  the times given with the units, and is found as well as they resolve it]
GLONASS_ROW = 10.0


@dataclass(frozen=True, slots=True)
class Carrier:
    """How one receiver message carries navigation units.

    ``message`` is the gnsswire class; ``protocol`` and ``number`` name it in a
    ``FrameFailure`` (rule S7); ``unit`` gives its ``gnsswire.nav.Received``;
    ``glonass`` the (slot, signal) of a GLONASS string, else None (read from the
    carrier's own fields, so that a string that fails its code is known);
    ``channel`` the GLONASS frequency channel + 7 (RTCM DF040) when the carrier
    gives it, else None; ``glonass_frame`` the GLONASS frame number within the
    superframe when the carrier gives it, else None.
    [design: one adapter for every carrier; the SkyTraq messages are another
    ``Carrier`` over ``gnsswire.nav.from_skytraq``]
    """

    message: type
    protocol: Protocol
    number: int
    unit: Callable[[Any], nav.Received]
    glonass: Callable[[Any], tuple[int, int] | None]
    channel: Callable[[Any], int | None]
    glonass_frame: Callable[[Any], int | None]


# [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 1.5.2 (Table 1), gnssId 6 is GLONASS]
_UBX_GLONASS = 6


def _sfrbx_glonass(m: ubx.RxmSfrbx) -> tuple[int, int] | None:
    """[spec: u-blox F9 HPG 1.51 1.5.4 and 3.17.9, u-blox X20 HPG 2.11 1.5.4 and 3.16.8,
    gnssId, svId (the GLONASS slot), sigId]"""
    return (m.sv_id, m.sig_id) if m.gnss_id == _UBX_GLONASS else None


def _sfrbx_channel(m: ubx.RxmSfrbx) -> int | None:
    """[spec: u-blox F9 HPG 1.51 3.17.9 and u-blox X20 HPG 2.11 3.16.8 (RXM-SFRBX
    freqId, "Only used for GLONASS: This is the frequency slot + 7 (range from 0 to
    13)"), which is RTCM DF040]"""
    return m.freq_id if m.gnss_id == _UBX_GLONASS else None


def _sfrbx_frame(m: ubx.RxmSfrbx) -> int | None:
    """[spec: UBX-13003221 R28 10.3, GLONASS word 4 holds the frame number in its 8 LSBs,
    0 when unknown]"""
    if m.gnss_id != _UBX_GLONASS or len(m.dwrd) != 4:
        return None
    return (m.dwrd[3] & 0xFF) or None


#: UBX-RXM-SFRBX. [spec: u-blox F9 HPG 1.51 3.17.9 and u-blox X20 HPG 2.11 3.16.8 (the
#: same message, of up to 16 words in both); UBX-13003221 R28 chapter 10 through
#: ``gnsswire.nav.from_sfrbx``; DESIGN 3.5, a UBX type is class × 256 + ID]
SFRBX = Carrier(ubx.RxmSfrbx, Protocol.UBX, (ubx.RxmSfrbx.MSG_CLASS << 8) | ubx.RxmSfrbx.MSG_ID,
                nav.from_sfrbx, _sfrbx_glonass, _sfrbx_channel, _sfrbx_frame)


@dataclass(slots=True)
class _Transport:
    """What one transport's units have built up. [design: per transport, so that two
    receivers never mix words]"""

    assemblies: dict[tuple, carriers.Assembler] = field(default_factory=dict)
    seen: dict[tuple, dict[tuple, GnssTime]] = field(default_factory=dict)
    almanac_even: dict[tuple, Any] = field(default_factory=dict)
    channel: dict[int, int] = field(default_factory=dict)
    almanac_health: dict[int, int] = field(default_factory=dict)
    leap: LeapSeconds | None = None


class NavData:
    """Navigation units of named transports to RTCM ephemeris messages, leap seconds,
    failures and refusals (module documentation).

    [design: DESIGN 2.2, Q15; the RTCM message is the hand-over so that one code
    converts and checks both sources]
    """

    def __init__(self, carriers_: tuple[Carrier, ...] = (SFRBX,)) -> None:
        """[design: the carriers taken, UBX-RXM-SFRBX by default]"""
        self.carriers = tuple(carriers_)
        self._transports: dict[str, _Transport] = {}

    @property
    def transports(self) -> tuple[str, ...]:
        """[design: DESIGN 2.2] the transports with assembly state, in order of first unit."""
        return tuple(self._transports)

    def end(self, transport: str) -> None:
        """Forget ``transport``'s state, as ``Framer.flush`` ends its bytes.
        [design: the state of a transport that has ended is not kept]"""
        self._transports.pop(transport, None)

    def push(self, message: Any, transport: str = "", near: GnssTime | None = None,
             received: GnssTime | None = None) -> list[Any]:
        """What one carrier message gives, in order (module documentation).

        ``near`` is a GPS time near the unit's receipt (the receipt time or the
        time base's current time, rule T5), or None when none is known yet: it
        ages the parts kept (``STALE_PARTS``) and resolves the week of an
        announced leap second. ``received`` is copied into a ``FrameFailure``.
        Raises ``TypeError`` for a message no carrier takes.
        [spec: gnsswire.nav (UBX-13003221 chapter 10; each system's ICD); design:
        DESIGN 2.1, a type that is not an input is a programming error]
        """
        carrier = next((c for c in self.carriers if isinstance(message, c.message)), None)
        if carrier is None:
            raise TypeError(f"not a navigation data carrier: {type(message).__name__}")
        state = self._transports.setdefault(transport, _Transport())
        try:
            unit = carrier.unit(message)
        except nav.Unsupported as e:
            # [design: DESIGN 2.2, words no staged document lays out are, to the engine,
            #  a system it does not support; refused visibly, not counted as a failure]
            return [FeedResult(FeedOutcome.REJECTED, FeedReason.UNSUPPORTED_SYSTEM, detail=str(e))]
        except nav.Error as e:
            # [spec: each ICD's code (IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-XIV,
            #  Galileo OS SIS ICD 2.2 5.1.9.4, BDS ICD 5.1.3, GLONASS ICD 4.7); design: DESIGN
            #  3.5, the carrier frame passed its checksum and its content cannot be read, a
            #  decode failure (rule S7); the unit never reaches an assembly]
            self._break(state, carrier, message)
            return [FrameFailure(FailureKind.DECODE, carrier.protocol, carrier.number, None, received)]
        out: list[Any] = []
        if unit.system == carriers.GLONASS:
            self._almanac(state, unit, carrier.glonass_frame(message))
        leap = _leap_seconds(unit, near)
        if leap is not None and leap != state.leap:
            state.leap = leap
            out.append(leap)
        eph = self._assembly(state, unit, near).push(unit)
        if eph is not None:
            out.append(self._rtcm(state, unit, eph, carrier.channel(message)))
        return out

    # ---- assembly ----------------------------------------------------------------------

    @staticmethod
    def _key(unit: nav.Received) -> tuple:
        """[design: one assembly per satellite of a system: Galileo E1-B and E5b-I words
        share it (Galileo OS SIS ICD 2.2 4.3.1, both carry I/NAV), BeiDou D1 and D2 share it (a
        satellite sends one, BDS ICD 5.1.1)]"""
        system = "beidou" if unit.system in (carriers.BEIDOU_D1, carriers.BEIDOU_D2) else unit.system
        return (system, unit.satellite)

    def _assembly(self, state: _Transport, unit: nav.Received, near: GnssTime | None) -> carriers.Assembler:
        """The satellite's assembly, started over when a part it holds arrived more
        than ``STALE_PARTS`` before ``near``, or, for GLONASS, when the last string
        arrived more than ``GLONASS_ROW`` before it. The time of each kind of part
        is its latest arrival, as the assembly keeps the latest of each kind.
        [design: STALE_PARTS, GLONASS_ROW]"""
        key = self._key(unit)
        if near is not None:
            parts = state.seen.setdefault(key, {})
            if unit.system == carriers.GLONASS:
                last = max(parts.values(), default=None)
                stale = last is not None and near - last > GLONASS_ROW
            else:
                stale = any(near - t > STALE_PARTS for t in parts.values())
            if stale:
                state.assemblies.pop(key, None)
                parts.clear()
            parts[_kind(unit)] = near
        assembly = state.assemblies.get(key)
        if assembly is None:
            assembly = state.assemblies[key] = carriers.Assembler()
        return assembly

    def _break(self, state: _Transport, carrier: Carrier, message: Any) -> None:
        """A GLONASS string that failed starts its stream's count over.

        [spec: GLONASS ICD 4.3, strings 1 to 5 of one frame follow each other;
        gnsswire.nav ``glonass.Assembler`` takes only strings 1 to 5 in a row;
        design: a failed string is a gap in the row, so strings before and after
        it, which may belong to different frames, are not joined]
        """
        stream = carrier.glonass(message)
        if stream is None:
            return
        state.assemblies.pop((carriers.GLONASS, stream[0]), None)
        state.almanac_even.pop(stream, None)

    def _almanac(self, state: _Transport, unit: nav.Received, frame: int | None) -> None:
        """The frequency channel and health of each slot, from almanac string pairs.

        An even almanac string (6 to 14) names the satellite nA and its Cn; the
        odd string right after it in the same stream gives its HnA, the carrier
        number, negative channels coded 25 to 31 for -7 to -1. Strings 14 and 15
        are taken only when the carrier gives the frame number, since in frame 5
        they are not almanac. [spec: GLONASS ICD 4.5, Tables 4.9, 4.10, 4.11;
        design: the channel and almanac health RTCM 1020 has (DF040, DF104) and
        the strings 1 to 5 lack, for a carrier that does not give the channel]
        """
        stream = (unit.satellite, unit.signal)
        message = unit.message
        m = message.m
        even = state.almanac_even.pop(stream, None)
        if m in (14, 15) and frame is None:
            return
        if isinstance(message, glonass.AlmanacEven):
            state.almanac_even[stream] = message
        elif isinstance(message, glonass.AlmanacOdd) and even is not None and even.m == m - 1:
            h = message.h
            channel = h - 32 if h >= 25 else h
            state.channel[even.na] = channel + 7
            state.almanac_health[even.na] = even.cn

    def _rtcm(self, state: _Transport, unit: nav.Received, eph: Any, channel: int | None) -> Any:
        """The RTCM message of an assembled ephemeris, or a refusal.

        [spec: gnsswire.nav ``rtcm_fields`` of each system (RTCM Tables 3.5-21,
        3.5-22, 3.5-111, 3.5-113); design: the satellite numbers in the data (Galileo
        word 4 SVID, Galileo OS SIS ICD 2.2 Table 45; GLONASS n, ICD 4.5 Table 4.5) must agree
        with the carrier's, or the words of two satellites were joined]
        """
        if isinstance(eph, gps.Ephemeris):
            return rtcm3.Msg1019(**gps.rtcm_fields(eph))
        if isinstance(eph, beidou.Ephemeris):
            return rtcm3.Msg1042(**beidou.rtcm_fields(eph))
        if isinstance(eph, galileo.InavEphemeris):
            if eph.satellite != unit.satellite:
                return _inconsistent(f"Galileo word 4 SVID {eph.satellite}, carrier SV {unit.satellite}")
            return rtcm3.Msg1046(**galileo.rtcm_fields(eph))
        # GLONASS. [spec: GLONASS ICD 4.5 Table 4.5 note (1), n is a GLONASS-M word;
        # design: 0, as an older satellite sends it, is replaced by the carrier's slot]
        n = eph.satellite or unit.satellite
        if n != unit.satellite:
            return _inconsistent(f"GLONASS string 4 n {eph.satellite}, carrier slot {unit.satellite}")
        if channel is None:
            channel = state.channel.get(n)
        if channel is None:
            # [design: the assembly is started over, so that the same issue is given
            #  again, at its next frame, once the almanac has given the channel]
            state.assemblies.pop(self._key(unit), None)
            return FeedResult(FeedOutcome.REJECTED, FeedReason.INCONSISTENT,
                              detail=f"GLONASS slot {n}: no frequency channel from the carrier or the "
                                     "almanac (RTCM DF040)")
        fields = glonass.rtcm_fields(eph, channel=channel, almanac_health=state.almanac_health.get(n))
        fields["satellite"] = n
        return rtcm3.Msg1020(**fields)


def _kind(unit: nav.Received) -> tuple:
    """The kind of part a unit is, of which an assembly keeps the latest.

    [spec: gnsswire.nav assemblers: GPS by subframe, Galileo by I/NAV word type,
    BeiDou by D1 subframe or D2 page, GLONASS by string number]
    """
    message = unit.message
    if unit.system == carriers.GALILEO_INAV:
        return ("word", message.word >> 122)       # Galileo OS SIS ICD 2.2 4.3.5: the first 6 bits are the type
    return (type(message).__name__, getattr(message, "m", None), getattr(message, "pnum", None))


def _inconsistent(detail: str) -> FeedResult:
    """[design: DESIGN 2.2, fields that cannot form an input are rejected, inconsistent]"""
    return FeedResult(FeedOutcome.REJECTED, FeedReason.INCONSISTENT, detail=detail)


# ---- leap seconds ------------------------------------------------------------------------

# [spec: BDS ICD 3.3, BDT is GPS time less 14 s, through toolkit.timescale (BDT week 0,
#  second 0 is GPS week 1356, second 14)]
_BDT_TO_GPS_SECONDS = round(timescale.bdt_to_gps(0, 0.0)[1])

# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.1.1 (Data ID number two,
#  binary 01, denotes the LNAV data structure of the appendix; IRN-IS-200N-003 removed
#  "and is the only valid value") and 6.4.6.4 (determine the Data ID before decoding
#  and using LNAV data); Table 20-V lays out the pages of that structure]
_GPS_LNAV_DATA_ID = 0b01


def _leap_seconds(unit: nav.Received, near: GnssTime | None) -> LeapSeconds | None:
    """GPS − UTC and an announced change, from a unit that carries them.

    GPS subframe 4 page 18 and Galileo word 6 give ΔtLS, ΔtLSF, WN_LSF (8 bits)
    and DN (1 to 7, "day one" the first of the week), the leap second at the end
    of day DN; their ΔtLS is GPS − UTC (GST and GPS time share whole seconds).
    BeiDou gives the same for BDT − UTC, DN from 0, so GPS − UTC is its count
    plus 14 s. When ΔtLS = ΔtLSF the effectivity time is not needed; otherwise
    WN_LSF is resolved nearest the week of ``near`` (none known: the current
    count alone, the page coming again within 12.5 min for GPS). A GPS page
    whose Data ID is not 2 is not read: its structure is not that of Table 20-V.
    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 (and its NOTE: when ΔtLS = ΔtLSF
    the effectivity time is not needed), 20.3.3.5.1.1 and 6.4.6.4 (``_GPS_LNAV_DATA_ID``);
    Galileo OS SIS ICD 2.2 5.1.7 Table 75 (DN 1 =
    Sunday to 7), 5.1.2 (GST was 13 s ahead of UTC at its start, as GPS time was, so
    GST − UTC is GPS − UTC); BDS ICD 5.2.4.18 Table 5-18 (DN 0 to 6), 3.3;
    ``adapters.ubx`` for the effectivity time, the end of day DN plus the new
    count]
    """
    message = unit.message
    if isinstance(message, gps.Subframe4Page18):
        if message.data_id != _GPS_LNAV_DATA_ID:
            return None
        offset, first_day, scale = 0, 1, "gps"
    elif unit.system == carriers.GALILEO_INAV and message.page_type == 0 and message.page_type_odd == 0:
        message = galileo.inav_word(message)
        if not isinstance(message, galileo.InavWord6):
            return None
        offset, first_day, scale = 0, 1, "gst"
    elif isinstance(message, beidou.Utc):
        offset, first_day, scale = _BDT_TO_GPS_SECONDS, 0, "bdt"
    else:
        return None
    current = message.dt_ls + offset
    if message.dt_lsf == message.dt_ls:
        return LeapSeconds(current)
    if near is None or not first_day <= message.dn <= first_day + 6:
        return LeapSeconds(current)
    # the start of week WN_LSF of the message's scale, on the GPS time scale
    if scale == "gps":
        start = GnssTime(timescale.resolve_week(message.wn_lsf, 8, near.week), 0.0)
    elif scale == "gst":
        wn = timescale.resolve_week(message.wn_lsf, 8, timescale.gps_to_gst(near.week, near.seconds)[0])
        start = GnssTime(*timescale.gst_to_gps(wn, 0.0))
    else:
        wn = timescale.resolve_week(message.wn_lsf, 8, timescale.gps_to_bdt(near.week, near.seconds)[0])
        start = GnssTime(*timescale.bdt_to_gps(wn, 0.0))
    end_of_day = (message.dn - first_day + 1) * _DAY
    return LeapSeconds(current, message.dt_lsf + offset, start + (end_of_day + message.dt_lsf))
