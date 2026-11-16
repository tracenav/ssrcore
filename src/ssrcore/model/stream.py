"""Stream identity and health: solution and processor, issues of update, ends of set, failures.

DESIGN.md section 3.3, rules S1 to S7. A *stream* is the data of one
(solution ID TF010, solution processor ID TF011) pair (S1). For each stream
this module keeps what SPARTN's continuity controls (SPARTN 6.2) need: the
issue of update of every message (SIOU SF005, and AIOU SF068 for the
atmosphere store), the gaps after which an equal issue counts as
new (S4, A2), the OCB set in progress and how the last one ended (S5), the
datum (S6), and which stream is active (S2). Frames that failed are counted
here (S7); they change no store.

Times are GPS times (``GnssTime``) of the messages' time tags, already
resolved and converted by the caller (rules T1 to T3, the time base's);
nothing here reads a clock.

The interface for other stores (DESIGN 3.3, S3): ``Streams.observe`` takes a
message's stream, time and issue numbers and returns ``Issues``; two blocks
may be combined exactly when their ``Issue`` values are equal.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ssrcore.model.types import (
    FeedOutcome,
    FeedReason,
    FeedResult,
    FrameFailure,
    FailureKind,
    GnssTime,
    OcbHeaderInput,
    OcbInput,
    Protocol,
    StreamId,
)

# [design: DESIGN 3.3 S7] failures are counted per (protocol, kind, message type, subtype).
FailureKey = tuple[str, str, int | None, int | None]

# [spec: SPARTN Table 6.2 SF005] "User must assume that SIOU has changed over
# outages of 510 secs or more even if the SIOU before/after an outage matches."
SIOU_LIFETIME = 510.0

# [spec: SPARTN Table 6.2 SF068] "User must assume that an AIOU has changed over
# outages of 450 seconds or more even if the AIOU before and after the outage matches."
AIOU_LIFETIME = 450.0

# [spec: SPARTN Table 6.1, message type 0 is OCB]
OCB_MESSAGE_TYPE = 0


@dataclass(frozen=True, slots=True, order=True)
class Issue:
    """An issue of update as received within one unbroken stretch of a stream.

    ``value`` is the SIOU (SF005) or AIOU (SF068) as sent; ``era`` counts
    the outages of the field's lifetime before it (rules S4, A2). Two blocks
    may be combined only when their issues are equal (rules S3, A1).
    [spec: SPARTN 6.2; Table 6.2 SF005, SF068]
    """

    value: int
    era: int


@dataclass(frozen=True, slots=True)
class Issues:
    """The issues of one message: its SIOU, and its AIOU when it carries one (HPAC, GAD).

    [spec: SPARTN Tables 6.3, 6.13, 6.23 (SIOU in OCB, HPAC, GAD; AIOU in HPAC, GAD)]
    """

    siou: Issue
    aiou: Issue | None = None


class _Lifetime:
    """One issue-of-update field of a stream: the era and the time of the last message carrying it."""

    __slots__ = ("lifetime", "era", "last", "value")

    def __init__(self, lifetime: float) -> None:
        self.lifetime = lifetime
        self.era = 0
        self.last: GnssTime | None = None
        self.value: int | None = None

    def stale(self, time: GnssTime) -> bool:
        # [spec: SPARTN Table 6.2 SF005, SF068] a message from before an outage of
        # the lifetime cannot be told apart from one of an earlier, different issue.
        return self.last is not None and self.last - time >= self.lifetime

    def issue(self, value: int, time: GnssTime) -> Issue:
        # [spec: SPARTN Table 6.2 SF005, SF068] a gap of the lifetime or more starts a new era.
        if self.last is not None and time - self.last >= self.lifetime:
            self.era += 1
        if self.last is None or time >= self.last:
            self.last = time
            self.value = value
        return Issue(value, self.era)


@dataclass(frozen=True, slots=True)
class OcbSet:
    """An OCB set that has ended (rule S5): its time, whether a message with
    end-of-set SF010 = 1 ended it, and whether a frame failed while it was open
    (DESIGN 3.5). [spec: SPARTN Table 6.2 SF010; 8.4 item 5]"""

    time: GnssTime
    end_of_set: bool
    incomplete: bool


@dataclass(frozen=True, slots=True)
class OcbObservation:
    """What one OCB message changed in its stream, for the corrections store.

    ``issue`` is the message's SIOU issue (S3, S4). ``datum_changed`` is
    rule S6's discontinuity. ``switched_from`` is the stream that was active
    before this message made its own stream active (S2), whose satellites
    therefore all have a discontinuity (C7).
    [design: DESIGN 3.3 S2, S3, S4, S6; the corrections store acts on them]
    """

    issue: Issue
    datum_changed: bool = False
    switched_from: StreamId | None = None


@dataclass(frozen=True, slots=True)
class StreamStatus:
    """The state of one stream for ``engine.status()`` and the record's stream quality.

    [spec: SPARTN Table 7.1 TF010, TF011; Table 6.2 SF005, SF009, SF010, SF068;
    design: DESIGN 2.5 quality indicators, with the end of set and failure counts]
    """

    stream: StreamId
    active: bool
    siou: int | None
    aiou: int | None
    datum: int | None
    latest: GnssTime | None
    latest_ocb: GnssTime | None
    open_set: GnssTime | None
    last_set: OcbSet | None
    datum_changes: int
    failures: dict[FailureKey, int] = field(default_factory=dict)
    notices: dict[FailureKey, int] = field(default_factory=dict)


class StreamState:
    """What is known of one stream (rule S1). Read through ``Streams``; changed only by it.

    [spec: SPARTN 6.2, Table 7.1 TF010, TF011]
    """

    def __init__(self, stream: StreamId) -> None:
        """[spec: SPARTN 6.2] an empty stream, before its first message."""
        self.stream = stream
        self.latest: GnssTime | None = None
        self.siou = _Lifetime(SIOU_LIFETIME)
        self.aiou = _Lifetime(AIOU_LIFETIME)
        self.latest_ocb: GnssTime | None = None
        self.datum: int | None = None
        self.datum_changes = 0
        self.open_set: GnssTime | None = None
        self.open_incomplete = False
        self.last_set: OcbSet | None = None

    def _note_time(self, time: GnssTime) -> None:
        if self.latest is None or time > self.latest:
            self.latest = time

    def _ocb_set(self, time: GnssTime, end_of_set: bool) -> None:
        # [spec: SPARTN Table 6.2 SF010, "all messages of the current type (including
        # all subtypes) have been sent for a single epoch"; design: DESIGN 3.3 S5, a
        # newer time tag also ends the set, without the flag]
        if self.open_set is not None and time > self.open_set:
            self.last_set = OcbSet(self.open_set, False, self.open_incomplete)
            self.open_set, self.open_incomplete = None, False
        if self.open_set is not None and time < self.open_set:
            return  # a late message of a set already passed
        if self.open_set is None and self.last_set is not None and time <= self.last_set.time:
            return  # a late message of a set already ended
        self.open_set = time
        if end_of_set:
            self.last_set = OcbSet(time, True, self.open_incomplete)
            self.open_set, self.open_incomplete = None, False


class Streams:
    """The stream store: every stream seen, the active one, and the frame failures.

    [spec: SPARTN 6.2; design: DESIGN 3.2, stream store keyed by (solution ID,
    processor ID)]
    """

    def __init__(self, stream_switch: float = 5.0) -> None:
        """[design: DESIGN 3.3 S2, ``Config.stream_switch`` seconds of lead before a switch]"""
        self.stream_switch = float(stream_switch)
        self._streams: dict[StreamId, StreamState] = {}
        self._active: StreamId | None = None
        self._switches = 0
        self._failures: dict[FailureKey, int] = {}
        self._notices: dict[FailureKey, int] = {}

    # ---- reading ----------------------------------------------------------------

    @property
    def active(self) -> StreamId | None:
        """The active stream, None before the first OCB message. [design: DESIGN 3.3 S2]"""
        return self._active

    @property
    def switches(self) -> int:
        """How many times the active stream has changed. [design: DESIGN 3.3 S2]"""
        return self._switches

    @property
    def failures(self) -> dict[FailureKey, int]:
        """Frames that failed, counted per (protocol, kind, message type, subtype).
        [design: DESIGN 3.3 S7]"""
        return dict(self._failures)

    @property
    def notices(self) -> dict[FailureKey, int]:
        """Messages used although they hold what a reader may meet in more than one way,
        counted per (protocol, notice, message type, subtype): ``extension``, data
        after the message's table; ``dnu_table_8_2``, a do-not-use satellite laid
        out as SPARTN Table 8.2 sizes it. [design: DESIGN 3.3 S9]"""
        return dict(self._notices)

    def streams(self) -> tuple[StreamId, ...]:
        """Every stream seen, in the order first seen. [spec: SPARTN 6.2, TF010, TF011]"""
        return tuple(self._streams)

    def state(self, stream: StreamId) -> StreamState | None:
        """The state of one stream, None if never seen. [spec: SPARTN 6.2, TF010, TF011]"""
        return self._streams.get(stream)

    def status(self, stream: StreamId) -> StreamStatus | None:
        """A snapshot of one stream's quality indicators. [design: DESIGN 2.5 stream quality]"""
        s = self._streams.get(stream)
        if s is None:
            return None
        return StreamStatus(stream, stream == self._active, s.siou.value, s.aiou.value, s.datum, s.latest,
                            s.latest_ocb, s.open_set, s.last_set, s.datum_changes, dict(self._failures),
                            dict(self._notices))

    # ---- feeding ----------------------------------------------------------------

    def _state(self, stream: StreamId) -> StreamState:
        s = self._streams.get(stream)
        if s is None:
            s = self._streams[stream] = StreamState(stream)
        return s

    def stale(self, stream: StreamId, time: GnssTime) -> bool:
        """Whether a message of ``stream`` at ``time`` is 510 s or more older than its newest.

        Such a message is from before an outage of the SIOU lifetime as far as
        this stream is concerned, so its SIOU cannot be compared with the
        current one. [spec: SPARTN Table 6.2 SF005; design: DESIGN 2.5, any
        SPARTN item is used only while its age is below 510 s]
        """
        s = self._streams.get(stream)
        return s is not None and s.latest is not None and s.latest - time >= SIOU_LIFETIME

    def observe(self, stream: StreamId, time: GnssTime, siou: int, aiou: int | None = None) -> Issues | None:
        """Note a SPARTN message of ``stream`` at ``time``; its issues of update, or None if stale.

        Rule S3: SIOU is solution-wide, so OCB, HPAC and GAD messages share one
        SIOU state per stream. Rule S4: after a gap of 510 s or more between a
        stream's messages, the SIOU starts a new era. Rule A2 likewise for the
        AIOU, over 450 s between the messages that carry one; a message with an
        AIOU that is 450 s or more older than the stream's latest such message
        is stale, since its era cannot be known.
        [spec: SPARTN 6.2; Table 6.2 SF005, SF068; design: gaps are measured
        between time tags, since receipt times may be absent (rule T5), and the
        AIOU gap over messages that carry an AIOU, since only they continue the
        area data]
        """
        if self.stale(stream, time):
            return None
        current = self._streams.get(stream)
        if aiou is not None and current is not None and current.aiou.stale(time):
            return None
        s = self._state(stream)
        s._note_time(time)
        issue = s.siou.issue(siou, time)
        area = None if aiou is None else s.aiou.issue(aiou, time)
        return Issues(issue, area)

    def observe_ocb(self, ocb: OcbInput, time: GnssTime) -> OcbObservation | None:
        """Note an OCB message at ``time`` (its tag in GPS time); None if stale.

        Rules S2 (active stream), S3 and S4 (issue), S5 (sets), S6 (datum).
        [spec: SPARTN 6.2; Table 6.2 SF005, SF009, SF010; design: DESIGN 3.3 S2]
        """
        issues = self.observe(ocb.stream, time, ocb.siou)
        if issues is None:
            return None
        s = self._streams[ocb.stream]
        newest = s.latest_ocb is None or time >= s.latest_ocb
        datum_changed = False
        if newest:
            # [spec: SPARTN Table 6.2 SF009, the datum is constant per solution and
            # processor ID; design: DESIGN 3.3 S6, a change is a discontinuity]
            datum_changed = s.datum is not None and ocb.datum != s.datum
            s.datum_changes += datum_changed
            s.datum = ocb.datum
            s.latest_ocb = time
        s._ocb_set(time, ocb.end_of_set)
        return OcbObservation(issues.siou, datum_changed, self._elect(ocb.stream))

    def _elect(self, stream: StreamId) -> StreamId | None:
        # [design: DESIGN 3.3 S2, the stream with the newest OCB time tag, kept
        # until another is newer by more than stream_switch]
        if self._active is None:
            self._active = stream
            return None
        if stream == self._active:
            return None
        lead = self._streams[stream].latest_ocb - self._streams[self._active].latest_ocb
        if lead > self.stream_switch:
            previous, self._active = self._active, stream
            self._switches += 1
            return previous
        return None

    def count_failure(self, failure: FrameFailure) -> FeedResult:
        """Count a frame that failed (rule S7); it changes no store.

        A failed OCB frame marks the OCB set in progress of every stream as
        incomplete, or, between sets, the next set to open, since the frame
        carries no readable stream key. [design: DESIGN 3.3 S7 and 3.5; a set's
        frames are sent together, so a failure between sets most likely
        belongs to the next]
        """
        key = (failure.protocol.value, failure.kind.value, failure.message_type, failure.subtype)
        self._failures[key] = self._failures.get(key, 0) + 1
        if failure.protocol == Protocol.SPARTN and failure.message_type == OCB_MESSAGE_TYPE:
            for s in self._streams.values():
                s.open_incomplete = True
        return FeedResult(FeedOutcome.COUNTED, FeedReason.FRAME_FAILURE, detail="/".join(map(str, key)))

    def count_notice(self, protocol: Protocol, notice: str, message_type: int | None,
                     subtype: int | None) -> None:
        """Count a message that is used but holds something to report (rule S9).

        [spec: SPARTN 8.4 item 4 and RTCM 3.5.17, data after a table is for later
        versions and is skipped; SPARTN Table 6.2 SF013 and Table 8.2, the two
        layouts of a DNU block; design: DESIGN 3.3 S9, used and made visible]
        """
        key = (Protocol(protocol).value, notice, message_type, subtype)
        self._notices[key] = self._notices.get(key, 0) + 1

    def observe_header(self, header: OcbHeaderInput) -> FeedResult:
        """The header of an OCB message whose subtype has no table: counted as an
        unsupported frame (S7) without marking its set incomplete, and its
        end-of-set flag ends the set open in its stream (S5).

        [spec: SPARTN 8.4 item 5, SF010 is read in every SM 0 message; Table 6.2
        SF010; design: DESIGN 3.3 S5, the subtype's time scale is not defined
        (Table 6.1, TBD), so its tag is not compared with the set's; a set's
        messages are sent together, so the flag ends the set that is open. Data
        the engine cannot use is not missing from the set it can use]
        """
        key = (Protocol.SPARTN.value, FailureKind.UNSUPPORTED.value, OCB_MESSAGE_TYPE, header.subtype)
        self._failures[key] = self._failures.get(key, 0) + 1
        s = self._streams.get(header.stream)
        ended = header.end_of_set and s is not None and s.open_set is not None
        if ended:
            s.last_set = OcbSet(s.open_set, True, s.open_incomplete)
            s.open_set, s.open_incomplete = None, False
        return FeedResult(FeedOutcome.COUNTED, FeedReason.UNSUPPORTED_SYSTEM,
                          detail="/".join(map(str, key)) + (", end of set" if ended else ""))
