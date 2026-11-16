"""The engine: feed typed messages with explicit times, query line-of-sight records.

DESIGN.md sections 2.1 to 2.5 and 2.8. ``Engine.feed`` takes any input of
DESIGN 2.2, a gnsswire message or its adapter's typed input, with the GPS time
it was received when the caller has it, and puts it into the stores of the
model (DESIGN 3.2) through the time base (rules T1 to T5). ``Engine.query``
reads the stores and assembles the ``LineOfSight`` record of DESIGN 2.4: the
selection and the transmit time are the model's (``model.selection``), every
term is a toolkit function of the selected data, and every absent term has
its reason.

The engine opens nothing and reads no clock; its whole state is its inputs
and its configuration (DESIGN 2.1, I/O at the edge), so a replay of the same
inputs gives the same records.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from gnsswire import nmea, rtcm3, spartn, ubx

from ssrcore.adapters import igs_ssr
from ssrcore.adapters import nmea as nmea_adapter
from ssrcore.adapters.navdata import NavData
from ssrcore.adapters import rtcm3 as rtcm3_adapter
from ssrcore.adapters import spartn_atmosphere, spartn_ocb
from ssrcore.adapters import ubx as ubx_adapter
from ssrcore.model.atmosphere import Absent, AtmosphereStore
from ssrcore.model.corrections import Corrections
from ssrcore.model.ephemeris import EphemerisStore, message_of
from ssrcore.model.igs_vtec import IgsVtecIonosphere, VtecInput
from ssrcore.model.providers import (Binding, Look, SlantTec, Where, Zenith, build, walk_ionosphere,
                                     walk_troposphere)
from ssrcore.model.selection import Choice, Transmit, choose, transmit
from ssrcore.model.stream import Streams, StreamStatus
from ssrcore.model.timebase import TimeBase
from ssrcore.model.types import (
    Carried,
    Config,
    FailureKind,
    FeedOutcome,
    FeedReason,
    FeedResult,
    FrameFailure,
    GadInput,
    GlonassEphemeris,
    GnssTime,
    HpacInput,
    KeplerEphemeris,
    LeapSeconds,
    Missing,
    OcbHeaderInput,
    OcbInput,
    Protocol,
    RoverPosition,
    Sat,
    Signal,
    StreamId,
    System,
    Term,
    UtcTime,
)
from ssrcore.record import (
    Broadcast,
    ClockTerm,
    CodeBiasTerm,
    Geometry,
    IonosphereTerm,
    LineOfSight,
    Modelled,
    OrbitTerm,
    PhaseBiasTerm,
    Quality,
    SatelliteQuality,
    SatelliteState,
    StreamQuality,
    TroposphereTerm,
)
from ssrcore.toolkit import KernelError, astro, geodesy, ionosphere, relativity, ssr, tides, timescale, troposphere
from ssrcore.toolkit import windup as _windup

Vector = tuple[float, float, float]

# [spec: SPARTN Table 6.2 SF024] the URE index whose value is a class, not a number.
_SF_URE = 24

# [design: DESIGN 2.4] the terms of a satellite, absent together when there is no state.
_SATELLITE_TERMS = (Term.SATELLITE, Term.GEOMETRY, Term.ORBIT, Term.CLOCK_CORRECTION, Term.EARTH_ROTATION,
                    Term.SHAPIRO, Term.RELATIVITY_CLOCK, Term.SOLID_TIDE, Term.WINDUP, Term.TROPOSPHERE,
                    Term.IONOSPHERE)
# [design: DESIGN 2.5, Config.elevation_cutoff] the terms left out below the cutoff:
# all but the geometry and the satellite state.
_BELOW_TERMS = (Term.EARTH_ROTATION, Term.SHAPIRO, Term.RELATIVITY_CLOCK, Term.SOLID_TIDE, Term.WINDUP,
                Term.TROPOSPHERE, Term.IONOSPHERE)

# [spec: RTCM Tables 3.5-21, 3.5-22, 3.5-110 to 3.5-113] the RTCM messages that are ephemerides.
_EPHEMERIS_MESSAGES = (rtcm3.Msg1019, rtcm3.Msg1020, rtcm3.Msg1042, rtcm3.Msg1044, rtcm3.Msg1045, rtcm3.Msg1046)


@dataclass(frozen=True, slots=True)
class EngineStatus:
    """What ``Engine.status`` reports: the active stream's indicators, every stream's
    state, the time base and the count of inputs fed; the configured satellite
    stream (S10, None for the active one) and the atmosphere providers' names in
    walk order (A9).

    [design: DESIGN 2.1, ``status()`` gives the stream state without a satellite]
    """

    active: StreamQuality | None
    streams: tuple[StreamStatus, ...]
    current: GnssTime | None
    leap_seconds: LeapSeconds | None
    held: int
    dropped: int
    fed: int
    satellite_stream: StreamId | None = None
    atmosphere: tuple[str, ...] = ()


@dataclass(slots=True)
class _Epoch:
    """What every query of one epoch (``t_rx``, receiver, store version) shares (DESIGN 2.8)."""

    key: tuple
    doy: float
    sky: astro.Sky | None = None
    tide: Vector | None = None
    atmosphere: dict = field(default_factory=dict)   # the providers' own per-epoch values (Where.cache)
    wheres: dict = field(default_factory=dict)       # binding or reason -> Where


@dataclass(frozen=True, slots=True)
class _Receiver:
    """A receiver position and what is derived from it alone (DESIGN 2.8)."""

    position: Vector
    geodetic: Vector
    east: Vector
    north: Vector
    latitude_deg: float
    longitude_deg: float


class Engine:
    """Feed inputs, query records. [design: DESIGN 2.1; the whole state is the
    inputs and the configuration, so a replay gives the same records]"""

    def __init__(self, config: Config | None = None, providers: dict | None = None) -> None:
        """[design: DESIGN 2.1, a frozen ``Config`` and, for its ``CallerSource``s, the
        caller's provider objects; nothing else configures the engine]"""
        self.config = Config() if config is None else config
        if not isinstance(self.config, Config):
            raise TypeError(f"expected Config, not {type(self.config).__name__}")
        c = self.config
        self.timebase = TimeBase(hold_limit=c.hold_limit)
        self.streams = Streams(c.stream_switch)
        self.corrections = Corrections(self.streams)
        self.atmosphere = AtmosphereStore(self.streams, c.max_age.atmosphere)
        self.providers = build(c, self.atmosphere, providers)
        self.ephemerides = EphemerisStore()
        self.navdata = NavData()
        self.igs = igs_ssr.IgsSsr(c.max_age)
        self.rover: RoverPosition | None = None
        self._fed = 0
        self._version = 0
        self._choices: dict[Sat, Choice] = {}
        self._choices_version = 0
        self._epoch: _Epoch | None = None
        self._receiver: _Receiver | None = None
        self.selections = 0          # selection work done, for DESIGN 2.8 and its test

    # ==== feed ======================================================================

    def feed(self, message: Any, received: GnssTime | None = None) -> FeedResult:
        """Put one input into the stores; the result says what it changed (DESIGN 2.1).

        ``message`` is a gnsswire message or a typed input of DESIGN 2.2;
        ``received`` the GPS time it was received, when the caller has it
        (rule T5). Never raises for bad data; raises ``TypeError`` for a type
        that is no input. Held inputs that the input makes placeable in time
        are released and stored in arrival order (T4, T5).
        [design: DESIGN 2.1, 2.2; rules T1 to T5]
        """
        if received is not None and not isinstance(received, GnssTime):
            raise TypeError(f"received must be a GnssTime, not {type(received).__name__}")
        self._fed += 1
        self._version += 1
        if received is not None:
            self.timebase.observe_received(received)
        result = self._feed(message, received)
        released = self._release()
        if released:
            detail = f"{released} held input(s) released"
            result = FeedResult(result.outcome, result.reason, result.sats, result.areas,
                                f"{result.detail}; {detail}" if result.detail else detail)
        return result

    def _release(self) -> int:
        """[design: DESIGN 3.3 rule T5, held inputs released in arrival order when
        the time or the leap seconds become known]"""
        count = 0
        while True:
            ready = self.timebase.release()
            if not ready:
                return count
            for item, received in ready:
                self._feed(item, received)
                count += 1

    def _hold(self, item: Any, received: GnssTime | None, reason: FeedReason, system: System | None,
              tag=None, needs_time: bool = False) -> FeedResult:
        # [design: DESIGN 3.3 rules T4, T5: held, at most hold_limit, the oldest dropped]
        dropped = self.timebase.hold((item, received), system, tag, received, needs_time)
        if dropped is not None:
            return FeedResult(FeedOutcome.HELD, FeedReason.HOLD_OVERFLOW,
                              detail=f"held for {reason.value}; the oldest held input was dropped")
        return FeedResult(FeedOutcome.HELD, reason)

    def _feed(self, message: Any, received: GnssTime | None) -> FeedResult:
        """[design: DESIGN 2.2, one branch per input; every gnsswire message is turned
        into its typed input by its adapter, so typed and wire inputs meet one path]"""
        if isinstance(message, spartn.Message):
            return self._spartn(message, received)
        if isinstance(message, (OcbInput, HpacInput, GadInput)):
            return self._spartn_input(message, received)
        if isinstance(message, OcbHeaderInput):
            return self.streams.observe_header(message)
        if isinstance(message, FrameFailure):
            return self.streams.count_failure(message)
        if igs_ssr.is_igs(message):
            return self._igs(message, received)
        if isinstance(message, igs_ssr.IgsReference):
            self.igs.declare(message)
            return FeedResult(FeedOutcome.STORED, detail=f"{message.stream}: {message.point.value}")
        if isinstance(message, VtecInput):
            return self._vtec(message)
        if isinstance(message, _EPHEMERIS_MESSAGES):
            self._notice(Protocol.RTCM3, message.extension, message.number, None)
            return self._rtcm_ephemeris(message, received)
        if isinstance(message, (KeplerEphemeris, GlonassEphemeris)):
            return self.ephemerides.add(message, received)
        if isinstance(message, LeapSeconds):
            self.timebase.set_leap_seconds(message)
            return FeedResult(FeedOutcome.STORED)
        if isinstance(message, ubx.NavTimels):
            return self._leap(message)
        if isinstance(message, UtcTime):
            return self._utc(message)
        if isinstance(message, nmea.Zda):
            return self._nmea(message, received, nmea_adapter.utc_time, self._utc, FeedReason.NO_TIME)
        if isinstance(message, RoverPosition):
            return self._rover(message)
        if isinstance(message, nmea.Gga):
            return self._nmea(message, received, nmea_adapter.rover_position, self._rover,
                              FeedReason.NOT_AN_INPUT)
        if isinstance(message, Carried):
            return self._navdata(message.message, message.transport, received)
        if isinstance(message, ubx.RxmSfrbx):
            return self._navdata(message, "", received)
        if isinstance(message, (rtcm3.Msg1005, rtcm3.Msg1006, rtcm3.Msg1033, rtcm3.Msg1230, rtcm3.Msm,
                                ubx.RxmRawx, ubx.RxmPmp, ubx.RxmPmpV0)):
            # [design: DESIGN 2.2, station and observation messages and receiver data
            # that are not inputs are refused visibly]
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT, detail=type(message).__name__)
        raise TypeError(f"not an input of the engine: {type(message).__name__}")

    # ---- SPARTN ---------------------------------------------------------------------

    def _notice(self, protocol: Protocol, extension: Any, message_type: int | None, subtype: int | None) -> None:
        """[spec: SPARTN 8.4 item 4; RTCM 3.5.17; design: DESIGN 3.3 S9, counted]"""
        if extension is not None:
            self.streams.count_notice(protocol, "extension", message_type, subtype)

    def _spartn(self, message: spartn.Message, received: GnssTime | None) -> FeedResult:
        """[spec: SPARTN 6.5, 8.4 items 4 and 5; design: DESIGN 3.5, an adapter that
        refuses a value counts the message as a failed frame (S7)]"""
        body, subtype = message.body, message.description.subtype
        if isinstance(body, spartn.Bpac):
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT,
                              detail="BPAC (SM 3-0) waits for a bpac kernel (DESIGN 1)")
        self._notice(Protocol.SPARTN, getattr(body, "extension", None), message.message_type, subtype)
        try:
            if isinstance(body, spartn.Ocb):
                if spartn_ocb.dnu_table(message) == "8.2":
                    self.streams.count_notice(Protocol.SPARTN, "dnu_table_8_2", 0, subtype)
                typed = spartn_ocb.ocb_input(message)
            elif isinstance(body, spartn.Hpac):
                typed = spartn_atmosphere.from_hpac(message)
            elif isinstance(body, spartn.Gad):
                typed = spartn_atmosphere.from_gad(message)
            else:
                raise TypeError(f"not a SPARTN input: {type(body).__name__}")
        except (ValueError, KeyError) as e:
            failure = FrameFailure(FailureKind.DECODE, Protocol.SPARTN, message.message_type, subtype, received)
            result = self.streams.count_failure(failure)
            return FeedResult(result.outcome, result.reason, detail=f"{result.detail}: {e}")
        return self._spartn_input(typed, received)

    def _spartn_input(self, item: OcbInput | HpacInput | GadInput, received: GnssTime | None) -> FeedResult:
        """[spec: SPARTN Table 7.1 TF009; rules T1 to T5 through the time base; S1 to
        S9, C3 to C10, A1 to A8 through the stores]"""
        system = System.GPS if isinstance(item, GadInput) else item.system
        self.timebase.observe_tag(item.tag, system)
        reason = self.timebase.missing(item.tag, system, received)
        if reason is not None:
            return self._hold(item, received, reason, system, item.tag)
        t = self.timebase.resolve(item.tag, system, received)
        if isinstance(item, OcbInput):
            return self.corrections.feed(item, t, received)
        if isinstance(item, HpacInput):
            return self.atmosphere.feed_hpac(item, t)
        return self.atmosphere.feed_gad(item, t)

    # ---- IGS SSR -----------------------------------------------------------------------

    def _igs(self, message: Any, received: GnssTime | None) -> FeedResult:
        """An IGS SSR message through ``adapters.igs_ssr``: satellite corrections into the
        corrections store as an OCB (rules S1 to S9, C3 to C11 unchanged), IM201 to the
        VTEC providers. [spec: IGS SSR v1.00 through the adapter (rules IG1 to IG12);
        design: DESIGN 2.2, the time near receipt resolves the week of IDF003, held
        until a time is known (T5); a block another message of its epoch completes is
        no conflict (S8, IG6)]"""
        self._notice(Protocol.RTCM3, message.extension, message.number, message.subtype)
        near = received if received is not None else self.timebase.current
        if near is None:
            return self._hold(message, received, FeedReason.NO_TIME, None, needs_time=True)
        adapted = self.igs.input(message, near, self.ephemerides)
        if isinstance(adapted, igs_ssr.Refused):
            reason = FeedReason.UNSUPPORTED_SYSTEM if adapted.unsupported_system else FeedReason.NOT_AN_INPUT
            return FeedResult(FeedOutcome.REJECTED, reason, detail=adapted.reason)
        if isinstance(adapted, igs_ssr.Kept):
            return FeedResult(FeedOutcome.DUPLICATE if adapted.duplicate else FeedOutcome.STORED,
                              detail=adapted.detail)
        if isinstance(adapted.item, VtecInput):
            return self._vtec(adapted.item)
        result = self._spartn_input(adapted.item, received)
        if adapted.completes and result.reason == FeedReason.CONFLICT:
            result = FeedResult(result.outcome, None, result.sats, result.areas,
                                f"{result.detail}; completes the block of its epoch (IG6)")
        return result

    def _vtec(self, vtec: VtecInput) -> FeedResult:
        """[spec: IGS SSR v1.00 4.5.1; design: DESIGN 3.6, the ionosphere of
        IM201 is a provider's (``model.igs_vtec``), which the caller configures; without
        one the message is no input]"""
        targets = [p for p in self.providers if isinstance(p, IgsVtecIonosphere)]
        if not targets:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT,
                              detail="IM201 needs an IgsVtecIonosphere provider (CallerSource)")
        results = [p.feed(vtec) for p in targets]
        return next((r for r in results if r.outcome == FeedOutcome.STORED), results[0])

    # ---- ephemerides, time, position ------------------------------------------------

    def _rtcm_ephemeris(self, message: Any, received: GnssTime | None) -> FeedResult:
        """[spec: RTCM 3.5.8, 3.5.9, 3.5.18 to 3.5.20 through ``adapters.rtcm3``; design:
        DESIGN 2.2, the time near receipt resolves the week; rule T4 for GLONASS;
        T5, held until a time is known]"""
        glonass = isinstance(message, rtcm3.Msg1020)
        system = System.GLONASS if glonass else None
        near = received if received is not None else self.timebase.current
        if near is None:
            return self._hold(message, received, FeedReason.NO_TIME, system, needs_time=True)
        leap = self.timebase.leap_at(near)
        if glonass and leap is None:
            return self._hold(message, received, FeedReason.NO_LEAP_SECONDS, system, needs_time=True)
        try:
            eph = rtcm3_adapter.ephemeris(message, near, leap)
        except rtcm3_adapter.EphemerisError as e:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.INCONSISTENT, detail=str(e))
        return self.ephemerides.add(eph, received)

    def _navdata(self, message: Any, transport: str, received: GnssTime | None) -> FeedResult:
        """Navigation words of a receiver, through ``adapters.navdata``: an assembled
        ephemeris is fed as the RTCM message of its system, leap seconds to the
        time base, a word that fails its code is counted (S7).

        The result is that of the last item the unit gave (an ephemeris, when
        there is one, is the last), with the details of every item.
        [spec: gnsswire.nav (UBX-13003221 R28 chapter 10; IS-GPS-200N with IRN-IS-200N-003 and -004
        20.3.3, Galileo OS SIS ICD 2.2 4.3.5, BDS ICD 5.2.4, GLONASS ICD 4.4, 4.5); design:
        DESIGN 2.2 and Q15, the RTCM message is the hand-over so that the
        conversion, the holds of T4 and T5 and rules E1 to E5 are one code for
        both sources; the words are assembled per transport]
        """
        if not isinstance(message, ubx.RxmSfrbx):
            raise TypeError(f"not a navigation data carrier: {type(message).__name__}")
        near = received if received is not None else self.timebase.current
        results = []
        for item in self.navdata.push(message, transport, near, received):
            if isinstance(item, _EPHEMERIS_MESSAGES):
                results.append(self._rtcm_ephemeris(item, received))
            elif isinstance(item, LeapSeconds):
                self.timebase.set_leap_seconds(item)
                results.append(FeedResult(FeedOutcome.STORED, detail="leap seconds"))
            elif isinstance(item, FrameFailure):
                results.append(self.streams.count_failure(item))
            else:                       # a rejected FeedResult
                results.append(item)
        if not results:
            return FeedResult(FeedOutcome.STORED, detail="navigation unit taken into its assembly")
        last = results[-1]
        if len(results) == 1:
            return last
        return FeedResult(last.outcome, last.reason, last.sats, last.areas,
                          "; ".join(r.detail for r in results if r.detail))

    def _leap(self, message: ubx.NavTimels) -> FeedResult:
        """[spec: u-blox F9 HPG 1.51 3.15.26 and u-blox X20 HPG 2.11 3.14.22 through
        ``adapters.ubx``; rule T4]"""
        try:
            leap = ubx_adapter.leap_seconds(message)
        except ValueError as e:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.INCONSISTENT, detail=str(e))
        if leap is None:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NO_LEAP_SECONDS,
                              detail="no valid leap second count from a trusted source (rule T4)")
        self.timebase.set_leap_seconds(leap)
        return FeedResult(FeedOutcome.STORED)

    def _utc(self, utc: UtcTime) -> FeedResult:
        """[design: DESIGN 3.3 rule T5, a ZDA date and time, only with leap seconds]"""
        if self.timebase.observe_utc(utc) is None:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NO_LEAP_SECONDS)
        return FeedResult(FeedOutcome.STORED)

    def _rover(self, rover: RoverPosition) -> FeedResult:
        """[design: DESIGN 2.2, the rover position is kept for the VRS and the
        caller; rule T5, its time of day placed in the known day]"""
        self.rover = rover
        self.timebase.observe_rover(rover)
        return FeedResult(FeedOutcome.STORED)

    def _nmea(self, sentence: Any, received: GnssTime | None, adapt, store, empty: FeedReason) -> FeedResult:
        """[spec: u-blox F9 HPG 1.51 2.7.5.1, 2.7.21.1 and u-blox X20 HPG 2.11 2.7.5.1,
        2.7.20.1 through ``adapters.nmea``; design: a malformed sentence is a failed
        frame (S7), one without a fix or a time carries nothing]"""
        try:
            value = adapt(sentence)
        except ValueError as e:
            result = self.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.NMEA, None, None, received))
            return FeedResult(result.outcome, result.reason, detail=f"{result.detail}: {e}")
        if value is None:
            return FeedResult(FeedOutcome.REJECTED, empty, detail="nothing usable in the sentence")
        return store(value)

    # ==== status =====================================================================

    def status(self) -> EngineStatus:
        """[design: DESIGN 2.1, 2.5] the stream state and the time base, without a satellite."""
        active = self.streams.active
        return EngineStatus(
            None if active is None else self._stream_quality(active),
            tuple(self.streams.status(s) for s in self.streams.streams()),
            self.timebase.current, self.timebase.leap_seconds, self.timebase.held,
            self.timebase.dropped, self._fed, self.config.satellite_stream,
            tuple(p.name for p in self.providers))

    def _stream_quality(self, stream: StreamId) -> StreamQuality:
        """[spec: SPARTN Table 7.1 TF010, TF011; Table 6.2 SF005, SF009, SF010, SF068;
        design: DESIGN 2.5 stream indicators, as received]"""
        st = self.streams.status(stream)
        last = st.last_set
        return StreamQuality(stream.solution, stream.processor, st.siou, st.aiou, st.datum,
                             bool(last and last.end_of_set), st.failures, bool(last and last.incomplete),
                             self.streams.switches, st.datum_changes, st.notices)

    # ==== query ======================================================================

    def query(self, sat: Sat, t_rx: GnssTime, receiver: Vector,
              signals: tuple[Signal, ...] = ()) -> LineOfSight:
        """The line-of-sight record of ``sat`` received at ``t_rx`` (GPS time) by a
        receiver at ``receiver`` (ECEF, m), with the per-signal terms of ``signals``
        (every signal with a bias when empty).

        Pure with respect to the stores (DESIGN 2.1): the caches it keeps
        (DESIGN 2.8) hold values that depend only on the stores and the
        arguments. [spec: DESIGN 3.4 steps 1 to 6 and the specification of each
        term, DESIGN 2.4; design: DESIGN 2.3, 2.5]
        """
        if not isinstance(sat, Sat) or not isinstance(t_rx, GnssTime):
            raise TypeError("query(sat: Sat, t_rx: GnssTime, receiver, signals)")
        rcv = self._receiver_of(receiver)
        signals = tuple(signals)
        missing: dict = {}
        base = dict(sat=sat, t_rx=t_rx, receiver=rcv.position)

        channel = self.ephemerides.glonass_channel(sat) if sat.system == System.GLONASS else None
        if sat.system == System.GLONASS and self.timebase.leap_seconds is None:
            # [spec: DESIGN 3.3 rule T4, GLONASS queries answer no_leap_seconds]
            return self._absent(base, signals, Missing.NO_LEAP_SECONDS, channel)

        choice = self._choice(sat, t_rx)
        if choice.broadcast is None:
            return self._absent(base, signals, choice.broadcast_missing, channel, choice)
        selection = choice.selection
        applied = choice.corrections_missing is None
        orbit_block = selection.orbit if applied else None
        tx = transmit(choice.broadcast, None if orbit_block is None else orbit_block.value, rcv.position, t_rx)
        if not isinstance(tx, Transmit):
            return self._absent(base, signals, tx, channel, choice)

        # ---- satellite state: SPARTN 8.7, 8.8 -----------------------------------------
        stored = choice.broadcast.stored
        eph = stored.ephemeris
        ref = eph.tb if isinstance(eph, GlonassEphemeris) else eph.toe
        clock_block = selection.clock if applied else None
        drift = choice.broadcast.drift(tx.t_tx)
        sat_clock = tx.clock
        if clock_block is not None:
            c = clock_block.value
            d_c = c.correction
            if c.polynomial:
                # [spec: RTCM 10403.3 3.5.13.6 eqs. (3.12-7), (3.12-8) and their derivative;
                # rule C11, at t_tx]
                dt = tx.t_tx - c.reference
                c1, c2 = c.rate or 0.0, c.acceleration or 0.0
                d_c = ssr.clock_at(c.correction, c1, c2, dt)
                drift = ssr.corrected_clock(drift, ssr.clock_rate_at(c1, c2, dt))
            sat_clock = ssr.corrected_clock(tx.clock, d_c)
        broadcast = Broadcast(tx.state.position, tx.clock, eph.iod, message_of(eph), ref, t_rx - ref)
        satellite = SatelliteState(tx.position, tx.state.velocity, sat_clock, drift, broadcast, channel)
        orbit = clock = None
        if applied:
            o = orbit_block.value
            orbit = OrbitTerm(*tx.orbit, tx.correction, o.iode, o.yaw, orbit_block.time,
                              t_rx - orbit_block.time, o.rate, o.reference)
            ure = ssr.field(_SF_URE, c.ure)
            clock = ClockTerm(d_c, c.iode_continuity, c.ure,
                              ure.value if ure.status == ssr.Status.OK else None,
                              clock_block.time, t_rx - clock_block.time, c.rate, c.acceleration, c.reference)
        else:
            missing[(Term.ORBIT, None)] = missing[(Term.CLOCK_CORRECTION, None)] = choice.corrections_missing

        # ---- geometry: TM-23 B.3 in the frame of reception (Ashby 2003 section 2) -----
        tau = t_rx - tx.t_tx
        arrival = geodesy.rotate_z(tx.position, geodesy.WGS84_OMEGA_E * tau)
        unit, elevation, azimuth = geodesy.line_of_sight(arrival, rcv.position, rcv.geodetic)
        geometry = Geometry(tx.range, tx.range_rate, unit, azimuth, elevation)

        record = dict(base, signals=signals, t_tx=tx.t_tx, geometry=geometry, satellite=satellite,
                      orbit=orbit, clock_correction=clock)
        signals = self._signals(signals, selection)
        record["signals"] = signals
        frequency = self._frequencies(sat, signals, channel, missing)
        record["frequency"] = frequency
        record["quality"] = Quality(None if choice.stream is None else self._stream_quality(choice.stream),
                                    self._satellite_quality(selection))
        if elevation < self.config.elevation_cutoff:
            for term in _BELOW_TERMS:
                missing[(term, None)] = Missing.BELOW_ELEVATION
            for s in signals:
                missing[(Term.CODE_BIAS, s)] = missing[(Term.PHASE_BIAS, s)] = Missing.BELOW_ELEVATION
            return LineOfSight(**record, missing=missing)

        epoch = self._epoch_of(t_rx, rcv)
        record["modelled"] = self._modelled(sat, tx, rcv, unit, selection, epoch, missing)
        record["code_bias"], record["phase_bias"] = self._biases(choice, signals, t_rx, missing)
        record["troposphere"], record["ionosphere"] = self._atmosphere(
            choice, sat, t_rx, rcv, elevation, azimuth, frequency, epoch, missing)
        return LineOfSight(**record, missing=missing)

    # ---- the parts of a query -------------------------------------------------------

    def _absent(self, base: dict, signals: tuple[Signal, ...], reason: Missing, channel: int | None,
                choice: Choice | None = None) -> LineOfSight:
        """[design: DESIGN 2.4, without a state every term is absent with one reason]"""
        missing: dict = {(term, None): reason for term in _SATELLITE_TERMS}
        signals = self._signals(signals, None if choice is None else choice.selection)
        for s in signals:
            missing[(Term.CODE_BIAS, s)] = missing[(Term.PHASE_BIAS, s)] = reason
        frequency = self._frequencies(base["sat"], signals, channel, missing)
        stream = None if choice is None or choice.stream is None else self._stream_quality(choice.stream)
        return LineOfSight(**base, signals=signals, frequency=frequency, missing=missing,
                           quality=Quality(stream, None if choice is None else
                                           self._satellite_quality(choice.selection)))

    @staticmethod
    def _signals(asked: tuple[Signal, ...], selection) -> tuple[Signal, ...]:
        """[design: DESIGN 2.3, the signals asked, or every signal with a bias]"""
        if asked or selection is None or selection.biases is None:
            return asked
        b = selection.biases.value
        return tuple(sorted({x.signal for x in b.code} | {x.signal for x in b.phase}))

    @staticmethod
    def _frequencies(sat: Sat, signals: tuple[Signal, ...], channel: int | None, missing: dict) -> dict:
        """[spec: Signal.frequency, IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.1.1, GLONASS ICD
        3.3.1.1 and the others of DESIGN 2.4 row λs]"""
        out = {}
        for s in signals:
            if s.system == System.GLONASS and channel is None:
                missing[(Term.FREQUENCY, s)] = Missing.NO_FREQUENCY
            else:
                out[s] = s.frequency(channel)
        return out

    def _receiver_of(self, receiver: Vector) -> _Receiver:
        """[spec: TM-23 B.1.2, B.2.1 through ``toolkit.geodesy``; design: DESIGN 2.8, the last receiver kept]"""
        position = (float(receiver[0]), float(receiver[1]), float(receiver[2]))
        r = self._receiver
        if r is not None and r.position == position:
            return r
        if len(receiver) != 3:
            raise ValueError("receiver is an ECEF (x, y, z)")
        try:
            geodetic = geodesy.ecef_to_geodetic(position)
        except KernelError as e:
            raise ValueError(f"receiver {position} has no geodetic position: {e}") from None
        east, north, _ = geodesy.enu_basis(geodetic[0], geodetic[1])
        self._receiver = _Receiver(position, geodetic, east, north,
                                   math.degrees(geodetic[0]), math.degrees(geodetic[1]))
        return self._receiver

    def _choice(self, sat: Sat, t_rx: GnssTime) -> Choice:
        """[design: DESIGN 2.8, a choice is kept while it holds and nothing is fed]"""
        if self._choices_version != self._version:
            self._choices.clear()
            self._choices_version = self._version
        kept = self._choices.get(sat)
        if kept is not None and kept.holds(t_rx):
            return kept
        self.selections += 1
        choice = choose(self.streams, self.corrections, self.ephemerides, sat, t_rx, self.config.max_age,
                        self.config.satellite_stream)
        self._choices[sat] = choice
        return choice

    def _epoch_of(self, t_rx: GnssTime, rcv: _Receiver) -> _Epoch:
        """[design: DESIGN 2.8, what the queries of one epoch share, the last epoch kept]"""
        key = (self._version, t_rx, rcv.position)
        e = self._epoch
        if e is not None and e.key == key:
            return e
        # [spec: TM-23 5.4.2.2.1 eq. (5.70), day of year from January 0.0; design: the
        # day of the GPS date, 18 s from the UTC date, nothing to Niell's seasons]
        day = t_rx.second // 86400
        doy = timescale.day_of_year(t_rx.week * 7 + day) + (t_rx.seconds - day * 86400.0) / 86400.0
        self._epoch = _Epoch(key, doy)
        return self._epoch

    def _tide(self, epoch: _Epoch, t_rx: GnssTime, rcv: _Receiver) -> Vector | None:
        """[spec: TN36 7.1.1 through ``toolkit.tides``, the Sun and Moon of
        ``toolkit.astro.sky``; design: DESIGN 2.4 row R_tide, UT1 taken as UTC]"""
        if epoch.tide is None:
            leap = self.timebase.leap_at(t_rx)
            if leap is None:
                return None
            epoch.sky = astro.sky(t_rx.week, t_rx.seconds, -float(leap))
            epoch.tide = tides.solid_earth_tide(rcv.position, epoch.sky)
        return epoch.tide

    def _modelled(self, sat: Sat, tx: Transmit, rcv: _Receiver, unit: Vector, selection, epoch: _Epoch,
                  missing: dict) -> Modelled:
        """E, S, dt_rel, R_tide and w of DESIGN 2.4."""
        receiver = rcv.position
        # [spec: Ashby 2003 section 10, Eq. (96), the logarithm; SPARTN 8.16]
        shapiro = relativity.shapiro(tx.position, receiver)
        # [spec: SPARTN 8.7 with r, ṙ of 8.8; zero for GLONASS, whose broadcast clock holds it]
        rel = 0.0 if sat.system == System.GLONASS else relativity.clock(tx.state.position, tx.state.velocity)
        displacement = self._tide(epoch, epoch.key[1], rcv)
        tide = None
        if displacement is None:
            missing[(Term.SOLID_TIDE, None)] = Missing.NO_LEAP_SECONDS
        else:
            # [spec: TN36 7.1.1; DESIGN 2.4, R_tide = −u·Δr]
            tide = -(unit[0] * displacement[0] + unit[1] * displacement[1] + unit[2] * displacement[2])
        # [spec: SPARTN 8.14 frame turned by the SF021 yaw (C8, zero when absent);
        # TM-23 5.5 with the receiver's east and north dipoles; design: the positions
        # of their own epochs, the frame of reception turning the line of sight by
        # 5 µrad, under 1e-6 cycle]
        yaw = 0.0 if selection is None else math.radians(selection.yaw)
        try:
            x, y, _ = _windup.satellite_frame(tx.position, tx.state.velocity, yaw)
            w = _windup.windup(tx.position, receiver, rcv.east, rcv.north, x, y)
        except KernelError:
            w = None
            missing[(Term.WINDUP, None)] = Missing.NOT_CONVERGED
        return Modelled(tx.sagnac, shapiro, rel, tide, displacement, w)

    def _biases(self, choice: Choice, signals: tuple[Signal, ...], t_rx: GnssTime,
                missing: dict) -> tuple[dict, dict]:
        """b_s and B_s [spec: SPARTN 8.9; rules C5, C7, C9, C10]."""
        code, phase = {}, {}
        selection = choice.selection
        if selection is None or selection.biases is None:
            reason = choice.corrections_missing if selection is None else selection.bias_missing
            for s in signals:
                missing[(Term.CODE_BIAS, s)] = missing[(Term.PHASE_BIAS, s)] = reason
            return code, phase
        block = selection.biases
        age = t_rx - block.time
        for s in signals:
            c = block.value.code_bias(s)
            if c is None:
                missing[(Term.CODE_BIAS, s)] = Missing.NO_BIAS
            else:
                code[s] = CodeBiasTerm(c.value, block.time, age)
            p = block.value.phase_bias(s)
            if p is None:
                missing[(Term.PHASE_BIAS, s)] = Missing.NO_BIAS
            else:
                phase[s] = PhaseBiasTerm(p.value, p.fixed, p.continuity,
                                         selection.phase_discontinuities.get(s, 0), block.time, age)
        return code, phase

    def _where(self, choice: Choice, t_rx: GnssTime, rcv: _Receiver, epoch: _Epoch) -> Where:
        """What the providers are asked about: the receiver, the epoch and the binding of
        the satellite corrections (rule A1), once per epoch and binding.
        [design: DESIGN 3.6; DESIGN 2.8, the epoch's dict is ``Where.cache``]"""
        selection = choice.selection
        key = choice.corrections_missing if selection is None else (selection.stream, selection.issue)
        where = epoch.wheres.get(key)
        if where is None:
            lat, lon, height = rcv.geodetic
            binding = None if selection is None else Binding(selection.stream, selection.issue)
            where = Where(t_rx, lat, lon, rcv.latitude_deg, rcv.longitude_deg, height, epoch.doy, binding,
                          choice.corrections_missing if selection is None else None, epoch.atmosphere)
            epoch.wheres[key] = where
        return where

    def _atmosphere(self, choice: Choice, sat: Sat, t_rx: GnssTime, rcv: _Receiver, elevation: float,
                    azimuth: float, frequency: dict, epoch: _Epoch, missing: dict):
        """T and I_s of DESIGN 2.4, each from the first provider that answers [spec: SPARTN
        8.10 to 8.12 through the default provider; design: rule A9]."""
        where = self._where(choice, t_rx, rcv, epoch)
        return (self._troposphere(walk_troposphere(self.providers, where), rcv, elevation, epoch, missing),
                self._ionosphere(walk_ionosphere(self.providers, where, sat, Look(elevation, azimuth)),
                                 choice.selection, frequency, missing))

    def _troposphere(self, zenith: Zenith | Absent, rcv: _Receiver, elevation: float, epoch: _Epoch,
                     missing: dict) -> TroposphereTerm | None:
        """[spec: SPARTN 8.12; TM-23 5.4.2.2.1 (Niell); UNB3m eq. (8) for a model's zenith
        delays; the provider brought the hydrostatic delay to the receiver height (DESIGN Q7)]"""
        if isinstance(zenith, Absent):
            missing[(Term.TROPOSPHERE, None)] = zenith.reason
            return None
        if not elevation > 0.0:
            # [spec: TM-23 5.4.2.2.1, the mapping functions are defined above the horizon]
            missing[(Term.TROPOSPHERE, None)] = Missing.BELOW_ELEVATION
            return None
        lat, _, height = rcv.geodetic
        mapping = troposphere.niell_mapping(lat, height, epoch.doy, elevation)
        s = zenith.spartn
        if s is None:
            slant = troposphere.slant(zenith.hydrostatic, zenith.wet, 0.0, mapping.hydrostatic, mapping.wet)
            return TroposphereTerm(slant, zenith.hydrostatic, mapping.hydrostatic, mapping.wet, None, None, None,
                                   None, zenith.age, zenith.wet, zenith.origin, zenith.sigma)
        if s.grid is None:
            missing[(Term.TROPOSPHERE_GRID, None)] = s.grid_missing
        slant = troposphere.slant(zenith.hydrostatic, s.polynomial, s.grid or 0.0, mapping.hydrostatic, mapping.wet)
        return TroposphereTerm(slant, zenith.hydrostatic, mapping.hydrostatic, mapping.wet, s.polynomial, s.grid,
                               s.quality, s.area_id, s.age, zenith.wet, zenith.origin, zenith.sigma)

    def _ionosphere(self, tec: SlantTec | Absent, selection, frequency: dict,
                    missing: dict) -> IonosphereTerm | None:
        """[spec: SPARTN 8.11; rules A6, A7, A11; TM-23 eq. (5.40) for a slant TEC of
        another provider, as for SPARTN's]"""
        if isinstance(tec, Absent):
            missing[(Term.IONOSPHERE, None)] = tec.reason
            return None
        s = tec.spartn
        if s is None:
            delays = {sg: ionosphere.delay(tec.tec, f, 0.0)[0] for sg, f in frequency.items()}
        else:
            if s.grid is None:
                missing[(Term.IONOSPHERE_GRID, None)] = s.grid_missing
            grid = s.grid or 0.0
            delays = {sg: ionosphere.delay(s.polynomial, f, grid)[0] for sg, f in frequency.items()}
        for (term, sg) in list(missing):
            if term == Term.FREQUENCY:
                missing[(Term.IONOSPHERE, sg)] = Missing.NO_FREQUENCY
        # [spec: SPARTN 8.6; design: rule A11, consistent with the code biases only when
        # from their own stream; None when the record has no code biases]
        with_biases = None
        if selection is not None and selection.biases is not None and selection.biases.value.code:
            with_biases = tec.origin.stream == selection.stream
        if s is None:
            return IonosphereTerm(tec.tec, None, None, delays, None, None, tec.age, tec.origin, tec.sigma, with_biases)
        return IonosphereTerm(s.slant, s.polynomial, s.grid, delays, s.quality, s.area_id, s.age, tec.origin,
                              tec.sigma, with_biases)

    @staticmethod
    def _satellite_quality(selection) -> SatelliteQuality | None:
        """[spec: SPARTN Table 6.2 SF015, SF022, SF023, SF024, SF008; design: DESIGN 2.5,
        the satellite indicators as received; ``fixed`` when every phase bias of the
        block says fixed, since SF023 is sent per phase bias]"""
        if selection is None:
            return None
        phase = () if selection.biases is None else selection.biases.value.phase
        clock = selection.clock.value
        return SatelliteQuality(clock.ure, bool(phase) and all(p.fixed for p in phase),
                                0.0 if selection.continuity is None else selection.continuity,
                                clock.iode_continuity, selection.yaw_present, selection.discontinuities)

