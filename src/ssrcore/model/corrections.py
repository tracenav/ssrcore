"""SPARTN orbit, clock and bias store (OCB, SM 0).

DESIGN.md section 3.3, rules C3 to C9 (C1 and C2 are the toolkit's ``ssr``
arithmetic, applied by the engine). Keyed by stream and satellite (DESIGN
3.2): data of different streams are never combined (S1).

Per satellite of a stream the store keeps, for each issue of update (S3,
S4), the orbit, clock and bias blocks received, each with its time tag (GPS
time), issue and receipt time; the satellite's DNU state (C6); and its arcs:
how many discontinuities the satellite and each phase bias have had (C7, and
S2, S6 for the whole stream). Blocks older than the SIOU lifetime of 510 s
are dropped, since nothing may use them (DESIGN 2.5).

``select`` answers, for a satellite and a time, the newest combination of
orbit and clock of one issue that may be used, with the biases of the same
issue (DESIGN 3.4 step 2), or why there is none.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Generic, TypeVar

from ssrcore.model.stream import SIOU_LIFETIME, Issue, Streams
from ssrcore.model.types import (
    ClockCorrection,
    CodeBias,
    FeedOutcome,
    FeedReason,
    FeedResult,
    GnssTime,
    MaxAge,
    Missing,
    OcbInput,
    OcbSatellite,
    OrbitCorrection,
    PhaseBias,
    Sat,
    Signal,
    StreamId,
)

T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class Block(Generic[T]):
    """One stored block: its value as sent, its time tag (GPS time), its issue and receipt time.

    [design: DESIGN 3.2, corrections store "each with SIOU, time tag, receipt time"]
    """

    value: T
    time: GnssTime
    issue: Issue
    received: GnssTime | None = None


@dataclass(frozen=True, slots=True)
class Biases:
    """A satellite's bias block (SPARTN Tables 6.7 to 6.12): phase and code biases by signal.

    The whole block replaces the previous one. [design: a bias block is the
    provider's set of biases at its time; a signal it no longer carries is
    not kept from an older block, which could be of another time and issue]
    """

    phase: tuple[PhaseBias, ...] = ()
    code: tuple[CodeBias, ...] = ()

    def phase_bias(self, signal: Signal) -> PhaseBias | None:
        """The phase bias of exactly ``signal``, None if the block has none (rule C9).

        [spec: SPARTN Table 6.2 SF025, SF026, SF102 to SF104; design: DESIGN 3.3
        C9, no substitution of another tracking mode of the band]
        """
        return next((b for b in self.phase if b.signal == signal), None)

    def code_bias(self, signal: Signal) -> CodeBias | None:
        """The code bias of exactly ``signal``, None if the block has none (rule C9).

        [spec: SPARTN Table 6.2 SF027, SF028, SF105 to SF107; design: DESIGN 3.3
        C9, no substitution of another tracking mode of the band]
        """
        return next((b for b in self.code if b.signal == signal), None)


@dataclass(frozen=True, slots=True)
class Selection:
    """The corrections of one satellite usable at a time (DESIGN 3.4 step 2).

    ``orbit`` and ``clock`` are of one issue (S3) and pass rule C4; ``biases``
    are of the same issue (C5), or None with ``bias_missing`` saying why.
    ``yaw`` is the yaw angle to use in degrees (C8); ``yaw_present`` whether
    the orbit block carried one. ``continuity`` is SF015 of the satellite's
    latest block in seconds; ``discontinuities`` the satellite's arc count
    and ``phase_discontinuities`` each phase bias signal's (C7).
    [spec: SPARTN 6.2, 8.9, 8.14; design: DESIGN 3.4]
    """

    stream: StreamId
    sat: Sat
    issue: Issue
    orbit: Block[OrbitCorrection]
    clock: Block[ClockCorrection]
    biases: Block[Biases] | None
    bias_missing: Missing | None
    yaw: float
    yaw_present: bool
    continuity: float | None
    discontinuities: int
    phase_discontinuities: dict[Signal, int] = field(default_factory=dict)


class _IssueBlocks:
    """The blocks of one satellite and one issue, each list in time-tag order."""

    __slots__ = ("orbit", "clock", "bias")

    def __init__(self) -> None:
        self.orbit: list[Block[OrbitCorrection]] = []
        self.clock: list[Block[ClockCorrection]] = []
        self.bias: list[Block[Biases]] = []

    def lists(self):
        return (self.orbit, self.clock, self.bias)


class _Satellite:
    """Everything stored for one satellite of one stream."""

    __slots__ = ("issues", "dnu", "last", "continuity", "arcs", "signal_arcs", "signal_last")

    def __init__(self) -> None:
        self.issues: dict[Issue, _IssueBlocks] = {}   # in the order first seen
        self.dnu: GnssTime | None = None
        self.last: GnssTime | None = None
        self.continuity: float | None = None
        self.arcs = 0
        self.signal_arcs: dict[Signal, int] = {}
        self.signal_last: dict[Signal, GnssTime] = {}

    def end_arc(self) -> None:
        # [design: DESIGN 3.3 C7; a discontinuity of the satellite ends the arc of
        # each of its phase biases too, so one signal's count changes at every
        # discontinuity that touches that signal]
        self.arcs += 1
        for signal in self.signal_arcs:
            self.signal_arcs[signal] += 1

    def has_data(self) -> bool:
        return any(any(lst) for b in self.issues.values() for lst in b.lists())


# Outcomes of putting one block, for the feed result.
_STORED, _DUPLICATE, _CONFLICT, _OLDER = "stored", "duplicate", "conflict", "older"


def _put(blocks: list[Block], new: Block) -> str:
    # [design: DESIGN 3.3 S8 at block level: an equal block of the same time tag
    # is a duplicate, a different one replaces it and is a conflict; C3, blocks
    # are kept in time order and an older block than the newest is refused]
    if blocks and new.time < blocks[-1].time:
        return _OLDER
    if blocks and new.time == blocks[-1].time:
        if blocks[-1].value == new.value:
            return _DUPLICATE
        blocks[-1] = new
        return _CONFLICT
    blocks.append(new)
    return _STORED


def _at(blocks: list[Block[T]], t: GnssTime) -> Block[T] | None:
    # [spec: SPARTN Table 7.1 TF009, the tag is the time of the data; design: DESIGN 2.5,
    # an item applies from its time tag until the next] the newest block not after t.
    for b in reversed(blocks):
        if b.time <= t:
            return b
    return None


def _fresh(block: Block, t: GnssTime, limit: float) -> bool:
    # [spec: SPARTN Table 6.2 SF005, age below 510 s; design: DESIGN 2.5 per-kind
    # limits and C3, a correction carries no rate so its age is bounded, not extrapolated]
    age = t - block.time
    return 0.0 <= age <= limit and age < SIOU_LIFETIME


class Corrections:
    """The corrections store: OCB orbit, clock and bias blocks by stream and satellite.

    It shares a ``Streams`` store, which keeps the issues, sets and the active
    stream. [design: DESIGN 3.2, corrections store keyed by stream and satellite]
    """

    def __init__(self, streams: Streams) -> None:
        """[design: DESIGN 3.2, one stream store serves the corrections and the atmosphere stores]"""
        self.streams = streams
        self._sats: dict[StreamId, dict[Sat, _Satellite]] = {}

    # ---- feeding ----------------------------------------------------------------

    def feed(self, ocb: OcbInput, time: GnssTime, received: GnssTime | None = None) -> FeedResult:
        """Store an OCB message whose time tag is ``time`` (GPS time, rules T1 to T3).

        Rules S1 to S6 through the stream store; per satellite C6 (DNU), C7
        (arcs) and the storing of each block under the message's issue (S3).
        [spec: SPARTN 6.2, 6.5.1, Table 6.2 SF013, SF014, SF015; design: DESIGN 3.3]
        """
        seen = self.streams.observe_ocb(ocb, time)
        if seen is None:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.STALE,
                              detail=f"time tag {SIOU_LIFETIME:.0f} s or more before the stream's newest")
        sats = self._sats.setdefault(ocb.stream, {})
        if seen.switched_from is not None:
            self._carry_arcs(seen.switched_from, ocb.stream)
        if seen.datum_changed:
            for s in sats.values():
                s.end_arc()
        counts = {_STORED: 0, _DUPLICATE: 0, _CONFLICT: 0, _OLDER: 0, "dnu": 0, "refused": 0, "continuity": 0}
        changed: list[Sat] = []
        for block in ocb.satellites:
            s = sats.get(block.sat)
            if s is None:
                s = sats[block.sat] = _Satellite()
            outcomes = self._satellite(s, block, time, seen.issue, received)
            for o in outcomes:
                counts[o] += 1
            if any(o in (_STORED, _CONFLICT, "dnu", "continuity") for o in outcomes):
                changed.append(block.sat)
            self._prune(s, self.streams.state(ocb.stream).latest)
        return self._result(counts, tuple(changed))

    def _satellite(self, s: _Satellite, block: OcbSatellite, time: GnssTime, issue: Issue,
                   received: GnssTime | None) -> list[str]:
        if block.dnu:
            # [spec: SPARTN Table 6.2 SF013, "Users must immediately stop using any
            # satellite that has been flagged"; design: DESIGN 3.3 C6, every stored
            # correction of the satellite is removed, and the arc ends]
            s.issues.clear()
            s.dnu = time if s.dnu is None or time > s.dnu else s.dnu
            s.end_arc()
            return ["dnu"]
        if s.dnu is not None:
            if time <= s.dnu:
                return ["refused"]  # [design: C6, a block from before the flag cannot lift it]
            s.dnu = None
        self._continuity(s, block, time)
        blocks = s.issues.setdefault(issue, _IssueBlocks())
        outcomes = []
        if block.orbit is not None:
            outcomes.append(_put(blocks.orbit, Block(block.orbit, time, issue, received)))
        if block.clock is not None:
            outcomes.append(_put(blocks.clock, Block(block.clock, time, issue, received)))
        if block.phase_biases is not None or block.code_biases is not None:
            biases = Biases(block.phase_biases or (), block.code_biases or ())
            outcome = _put(blocks.bias, Block(biases, time, issue, received))
            if outcome in (_STORED, _CONFLICT):
                self._phase_continuity(s, biases.phase, time)
            outcomes.append(outcome)
        if not any(blocks.lists()):
            del s.issues[issue]
        # [spec: SPARTN Table 6.2 SF014: a block may carry no orbit, clock or bias;
        # its continuity indicator SF015 still counts for C7]
        return outcomes or ["continuity"]

    @staticmethod
    def _continuity(s: _Satellite, block: OcbSatellite, time: GnssTime) -> None:
        # [spec: SPARTN 6.2 and Table 6.2 SF015: "the data has had no discontinuity
        # within the last X seconds"; the state is reset "when a data gap ... greater
        # than the continuity time has occurred". SF015 = 0 gives X = 0, so every
        # newer block ends the arc]
        if s.last is not None and time > s.last and block.continuity is not None \
                and time - s.last > block.continuity:
            s.end_arc()
        if s.last is None or time >= s.last:
            s.last = time
            s.continuity = block.continuity

    @staticmethod
    def _phase_continuity(s: _Satellite, phase: tuple[PhaseBias, ...], time: GnssTime) -> None:
        # [spec: SPARTN Table 6.12 and Table 6.2 SF015, the continuity of each phase bias]
        for b in phase:
            last = s.signal_last.get(b.signal)
            s.signal_arcs.setdefault(b.signal, 0)
            if last is not None and time > last and time - last > b.continuity:
                s.signal_arcs[b.signal] += 1
            if last is None or time >= last:
                s.signal_last[b.signal] = time

    def _carry_arcs(self, old: StreamId, new: StreamId) -> None:
        # [design: DESIGN 3.3 S2, a switch of the active stream is a discontinuity for
        # every satellite; each count of the new stream is raised above the old
        # stream's, so the counts the record gives only ever increase]
        before = self._sats.get(old, {})
        after = self._sats.setdefault(new, {})
        for sat in set(before) | set(after):
            o = before.get(sat)
            n = after.get(sat)
            if n is None:
                n = after[sat] = _Satellite()
            n.arcs = max(n.arcs, o.arcs if o else 0) + 1
            for signal in set(n.signal_arcs) | set(o.signal_arcs if o else ()):
                n.signal_arcs[signal] = max(n.signal_arcs.get(signal, 0),
                                            o.signal_arcs.get(signal, 0) if o else 0) + 1

    @staticmethod
    def _prune(s: _Satellite, latest: GnssTime | None) -> None:
        # [spec: SPARTN Table 6.2 SF005; DESIGN 2.5: no SPARTN item is used at an
        # age of 510 s or more, so blocks that old relative to the stream's newest
        # time are dropped]
        if latest is None:
            return
        for issue in list(s.issues):
            blocks = s.issues[issue]
            for lst in blocks.lists():
                lst[:] = [b for b in lst if latest - b.time < SIOU_LIFETIME]
            if not any(blocks.lists()):
                del s.issues[issue]

    @staticmethod
    def _result(counts: dict[str, int], sats: tuple[Sat, ...]) -> FeedResult:
        # [design: DESIGN 2.1, the feed result makes every effect visible]
        detail = ", ".join(f"{k}={v}" for k, v in counts.items() if v)
        if counts[_STORED] or counts[_CONFLICT] or counts["dnu"] or counts["continuity"]:
            reason = FeedReason.CONFLICT if counts[_CONFLICT] else FeedReason.DNU if counts["dnu"] else None
            return FeedResult(FeedOutcome.STORED, reason, sats, detail=detail)
        if counts["refused"]:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.DNU, detail=detail)
        if counts[_OLDER]:
            return FeedResult(FeedOutcome.SUPERSEDED, FeedReason.OLDER, detail=detail)
        return FeedResult(FeedOutcome.DUPLICATE, detail=detail)

    # ---- reading ----------------------------------------------------------------

    def satellites(self, stream: StreamId) -> tuple[Sat, ...]:
        """The satellites of ``stream`` with any stored block, in order. [spec: SPARTN SF011, SF012, SF093 to SF095]"""
        return tuple(sorted(sat for sat, s in self._sats.get(stream, {}).items() if s.has_data()))

    def dnu(self, stream: StreamId, sat: Sat) -> GnssTime | None:
        """The time tag of the DNU flag in force for ``sat``, None if none. [spec: SPARTN Table 6.2 SF013]"""
        s = self._sats.get(stream, {}).get(sat)
        return None if s is None else s.dnu

    def discontinuities(self, stream: StreamId, sat: Sat) -> int:
        """The satellite's arc count (rule C7, with S2 and S6). [spec: SPARTN 6.2, Table 6.2 SF015]"""
        s = self._sats.get(stream, {}).get(sat)
        return 0 if s is None else s.arcs

    def phase_discontinuities(self, stream: StreamId, sat: Sat, signal: Signal) -> int:
        """One phase bias signal's arc count (rule C7). [spec: SPARTN Table 6.12, Table 6.2 SF015]"""
        s = self._sats.get(stream, {}).get(sat)
        return 0 if s is None else s.signal_arcs.get(signal, 0)

    def blocks(self, stream: StreamId, sat: Sat) -> dict[Issue, tuple[tuple[Block, ...], ...]]:
        """Every stored (orbit, clock, bias) block list of ``sat`` by issue, oldest issue first.

        [design: DESIGN 3.2, what the store holds, for inspection and tests]
        """
        s = self._sats.get(stream, {}).get(sat)
        if s is None:
            return {}
        return {issue: tuple(tuple(lst) for lst in b.lists()) for issue, b in s.issues.items()}

    def next_change(self, stream: StreamId, sat: Sat, t: GnssTime) -> GnssTime | None:
        """The earliest time tag after ``t`` of any block stored for ``sat``, None if none.

        Until then, and while nothing is fed, ``select`` at a later time can
        change only by the age limits of the blocks it chose (C3), since every
        block it may take applies from its tag (TF009). [design: DESIGN 2.8,
        the interval over which a selection is kept]
        """
        s = self._sats.get(stream, {}).get(sat)
        if s is None:
            return None
        later = [b.time for blocks in s.issues.values() for lst in blocks.lists() for b in lst if b.time > t]
        return min(later) if later else None

    def select(self, stream: StreamId, sat: Sat, t: GnssTime, max_age: MaxAge) -> Selection | Missing:
        """The newest usable orbit and clock of one issue for ``sat`` at ``t``, or why there is none.

        Issues are tried newest first (DESIGN 3.4 step 2). In each, the orbit
        and clock blocks that apply at ``t`` (the newest not after it) must
        both be within their age limits (C3) and the clock must not be later
        than the orbit by more than its IODE continuity SF022 (C4). Biases
        are those of the same issue (C5), within their own age limit. A
        satellite flagged DNU has nothing (C6).
        [spec: SPARTN 6.2; Table 6.2 SF005, SF013, SF022; design: DESIGN 3.4
        step 2, newest complete combination rather than newest message]
        """
        s = self._sats.get(stream, {}).get(sat)
        if s is None:
            return Missing.NO_CORRECTIONS
        if s.dnu is not None:
            return Missing.DNU
        reasons: set[Missing] = set()
        any_orbit = any_clock = False
        for issue in reversed(list(s.issues)):
            blocks = s.issues[issue]
            orbit, clock = _at(blocks.orbit, t), _at(blocks.clock, t)
            any_orbit |= orbit is not None
            any_clock |= clock is not None
            if orbit is None or clock is None:
                continue
            if not (_fresh(orbit, t, max_age.orbit) and _fresh(clock, t, max_age.clock)):
                reasons.add(Missing.STALE)
                continue
            # [spec: SPARTN 6.2 and Table 6.2 SF022: the same IODE was used over the
            # last X seconds before the clock's time; an orbit older than that may
            # name an IODE that has since changed, so a new one must be received]
            if clock.time - orbit.time > clock.value.iode_continuity:
                reasons.add(Missing.IODE_CHANGE)
                continue
            biases, bias_missing = self._biases(s, issue, t, max_age.bias)
            yaw = orbit.value.yaw
            return Selection(stream, sat, issue, orbit, clock, biases, bias_missing,
                             # [spec: SPARTN 8.14, "If the yaw angle is not provided, the
                             # satellite yaw angle should be considered zero"; Table 6.2
                             # SF008, SF021; DESIGN 3.3 C8]
                             0.0 if yaw is None else yaw, yaw is not None,
                             s.continuity, s.arcs, dict(s.signal_arcs))
        if Missing.IODE_CHANGE in reasons:
            return Missing.IODE_CHANGE
        if Missing.STALE in reasons:
            return Missing.STALE
        if any_orbit and any_clock:
            return Missing.SIOU_MISMATCH
        return Missing.NO_CORRECTIONS

    @staticmethod
    def _biases(s: _Satellite, issue: Issue, t: GnssTime, limit: float) -> tuple[Block[Biases] | None, Missing | None]:
        # [spec: SPARTN 6.2, SF005: parts sent separately combine only under one SIOU;
        # design: DESIGN 3.3 C5 and 3.4 step 2, biases of another issue are missing,
        # not substituted]
        bias = _at(s.issues[issue].bias, t)
        if bias is not None:
            return (bias, None) if _fresh(bias, t, limit) else (None, Missing.STALE)
        other = any(_at(b.bias, t) is not None for i, b in s.issues.items() if i != issue)
        return None, Missing.SIOU_MISMATCH if other else Missing.NO_BIAS
