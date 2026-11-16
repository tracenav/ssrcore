"""The stream store: rules S1 to S7 (DESIGN 3.3), on OCB sequences encoded with gnsswire.

Messages are built from chosen values, encoded with ``gnsswire.spartn.encode``,
parsed, decoded and converted by ``adapters.spartn_ocb.ocb_input``. Each
message's GPS time is given by the test: the 32-bit tag's epoch, from
``toolkit.timescale``, plus the tag (rules T1, T3 on the GPS scale).
"""

from __future__ import annotations

from dataclasses import replace

from gnsswire import spartn

from ssrcore.adapters.spartn_ocb import ocb_input
from ssrcore.model.corrections import Corrections, Selection
from ssrcore.model.stream import AIOU_LIFETIME, SIOU_LIFETIME, Issue, Issues, OcbSet, Streams
from ssrcore.model.types import (
    FailureKind,
    FeedOutcome,
    FeedReason,
    FrameFailure,
    GnssTime,
    MaxAge,
    Missing,
    OcbHeaderInput,
    Protocol,
    SpartnTag,
    Sat,
    StreamId,
    System,
)
from ssrcore.toolkit import timescale

TAG0 = 500_000_000
T0 = GnssTime(*timescale.spartn_tag_to_gps(TAG0, System.GPS))     # 2010-01-01 GPS time + tag

SAT_BITS = {0: 32, 1: 24, 2: 36, 3: 37, 4: 10}


def ocb(numbers=(5,), *, system=0, tag=0, solution=1, processor=0, siou=1, eos=False, datum=0, continuity=7):
    """An OCB message with orbit (IODE 40, zero corrections) and clock for each satellite, through the wire."""
    mask = sum(1 << (SAT_BITS[system] - n) for n in numbers)
    sats = tuple(spartn.OcbSatellite(0, 0b110, continuity, spartn.Orbit(40, 8191, 8191, 8191, None),
                                     spartn.Clock(7, 8191 + 50, 1), None) for _ in sorted(numbers))
    body = spartn.Ocb(siou, int(eos), 0, 0, datum, 0, 0, mask, sats)
    d = spartn.Description(system, 1, TAG0 + tag, solution, processor, None, None, None, None)
    frame = spartn.parse(spartn.encode(spartn.Message(0, 0, 0, d, body)), 0)
    return ocb_input(spartn.decode(frame))


class Feed:
    def __init__(self, stream_switch=5.0):
        self.streams = Streams(stream_switch)
        self.store = Corrections(self.streams)

    def __call__(self, numbers=(5,), tag=0, **kw):
        return self.store.feed(ocb(numbers, tag=tag, **kw), T0 + tag)


A, B = StreamId(1, 0), StreamId(1, 1)


def test_rule_S1_streams_are_never_combined():
    """S1, SPARTN 6.2 and TF010, TF011: (solution, processor) is the stream; no data cross between streams."""
    f = Feed()
    f((5,), solution=1, processor=0)
    f((6,), solution=1, processor=1)
    f((7,), solution=2, processor=0)
    assert f.streams.streams() == (A, B, StreamId(2, 0))
    assert f.store.satellites(A) == (Sat.gps(5),) and f.store.satellites(B) == (Sat.gps(6),)
    assert f.store.select(B, Sat.gps(5), T0, MaxAge()) == Missing.NO_CORRECTIONS
    # an orbit in one stream and a clock in another give nothing
    g = Feed()
    msg = ocb((5,), processor=0)
    orbit_only = msg.satellites[0].__class__(Sat.gps(5), False, 320.0, orbit=msg.satellites[0].orbit)
    clock_only = msg.satellites[0].__class__(Sat.gps(5), False, 320.0, clock=msg.satellites[0].clock)
    g.store.feed(replace(msg, satellites=(orbit_only,)), T0)
    g.store.feed(replace(msg, stream=B, satellites=(clock_only,)), T0)
    assert g.store.select(A, Sat.gps(5), T0, MaxAge()) == Missing.NO_CORRECTIONS
    assert g.store.select(B, Sat.gps(5), T0, MaxAge()) == Missing.NO_CORRECTIONS


def test_rule_S2_active_stream_switches_only_on_a_clear_lead():
    """S2: the newest OCB time tag wins, but only by more than stream_switch; a switch is a discontinuity."""
    f = Feed(stream_switch=5.0)
    f((5, 6), tag=0, processor=0)
    assert f.streams.active == A
    f((5,), tag=3, processor=1)                      # B ahead by 3 s
    f((5,), tag=5, processor=1)                      # by 5 s: not more than 5 s
    assert f.streams.active == A and f.streams.switches == 0
    f((5,), tag=10, processor=0)                     # A newer again
    f((5,), tag=14, processor=1)                     # B ahead by 4 s
    assert f.streams.active == A
    before = {n: f.store.discontinuities(A, Sat.gps(n)) for n in (5, 6)}
    f((5,), tag=16, processor=1)                     # by 6 s: switch
    assert f.streams.active == B and f.streams.switches == 1
    for n in (5, 6):                                 # every satellite, also one B has not sent
        assert f.store.discontinuities(B, Sat.gps(n)) > before[n]
    assert f.streams.status(B).active and not f.streams.status(A).active
    f((5,), tag=20, processor=0)                     # A back but not ahead: stays B
    assert f.streams.active == B


def test_rule_S3_siou_is_solution_wide():
    """S3, SPARTN SF005: one SIOU state per stream for OCB, HPAC and GAD; blocks combine on equal SIOU."""
    s = Streams()
    o = s.observe_ocb(ocb(siou=9), T0)
    assert o.issue == Issue(9, 0)
    assert s.observe(A, T0 + 1, 9, aiou=3) == Issues(Issue(9, 0), Issue(3, 0))   # an HPAC of the same issue
    assert s.observe(A, T0 + 2, 10).siou != o.issue                              # a GAD of another
    assert s.status(A).siou == 10 and s.status(A).aiou == 3


def test_rule_S4_a_gap_of_510_s_makes_an_equal_siou_new():
    """S4, SPARTN SF005: after an outage of 510 s or more the SIOU is new even if equal."""
    f = Feed()
    f((5,), tag=0, siou=4)
    assert f.streams.observe(A, T0 + 509, 4).siou == Issue(4, 0)
    f((5,), tag=509 + int(SIOU_LIFETIME), siou=4)                           # a gap of exactly 510 s
    sel = f.store.select(A, Sat.gps(5), T0 + 1019, MaxAge())
    assert isinstance(sel, Selection) and sel.issue == Issue(4, 1)
    r = f((5,), tag=509, siou=4)                                           # 510 s before the newest
    assert r.outcome == FeedOutcome.REJECTED and r.reason == FeedReason.STALE


def test_rule_S4_stale_parts_never_combine_across_an_outage():
    """S4 with S3: an orbit before an outage and a clock after it, same SIOU, do not combine."""
    f = Feed()
    f((5,), tag=0, siou=4)
    late = ocb((5,), siou=4)
    clock_only = late.satellites[0].__class__(Sat.gps(5), False, 320.0, clock=late.satellites[0].clock)
    f.store.feed(replace(late, satellites=(clock_only,)), T0 + 600)
    assert list(f.store.blocks(A, Sat.gps(5))) == [Issue(4, 1)]
    assert f.store.select(A, Sat.gps(5), T0 + 600, MaxAge()) == Missing.NO_CORRECTIONS


def test_rule_A2_interface_aiou_gap_of_450_s():
    """A2 through the stream store (for the atmosphere store): the AIOU gap is over messages that carry an AIOU."""
    s = Streams()
    assert s.observe(A, T0, 1, aiou=2).aiou == Issue(2, 0)
    s.observe(A, T0 + 300, 1)                                                # an OCB-like message, no AIOU
    assert s.observe(A, T0 + 449, 1, aiou=2).aiou == Issue(2, 0)
    assert s.observe(A, T0 + 449 + AIOU_LIFETIME, 1, aiou=2).aiou == Issue(2, 1)


def test_rule_S5_end_of_set_and_lost_end_of_set():
    """S5, SPARTN SF010 and 8.4 item 5: a set ends with the flag, over all subtypes, or by a newer tag."""
    f = Feed()
    f((5,), tag=0, system=0, eos=False)
    assert f.streams.status(A).open_set == T0 and f.streams.status(A).last_set is None
    assert isinstance(f.store.select(A, Sat.gps(5), T0, MaxAge()), Selection)   # usable before the end of set
    f((3,), tag=0, system=2, eos=True)                                           # Galileo ends the set
    assert f.streams.status(A).last_set == OcbSet(T0, True, False)
    assert f.streams.status(A).open_set is None
    f((5,), tag=5, system=0, eos=False)                                          # its end-of-set message is lost
    f((5,), tag=10, system=0, eos=False)                                         # the newer tag ends it
    assert f.streams.status(A).last_set == OcbSet(T0 + 5, False, False)
    f((7,), tag=5, system=0, eos=True)                                           # a late message of that set
    assert f.streams.status(A).last_set == OcbSet(T0 + 5, False, False)
    assert f.streams.status(A).open_set == T0 + 10


def test_rule_S5_failed_frame_marks_the_set_incomplete():
    """S5 with DESIGN 3.5: a failed OCB frame marks the open set, or the next one, incomplete."""
    f = Feed()
    f((5,), tag=0, eos=False)
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 0, 0))
    f((6,), tag=0, eos=True)
    assert f.streams.status(A).last_set == OcbSet(T0, True, True)
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 0, 1))  # between sets
    f((5,), tag=5, eos=True)
    assert f.streams.status(A).last_set == OcbSet(T0 + 5, True, True)
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 1, 0))  # not OCB: no effect
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.RTCM3, 0))     # RTCM, not SPARTN OCB
    f((5,), tag=10, eos=True)
    assert f.streams.status(A).last_set == OcbSet(T0 + 10, True, False)


def test_rule_S5_end_of_set_read_from_a_subtype_without_a_table():
    """S5 with SPARTN 8.4 item 5 (DESIGN Q19): the header of SM 0-9 ends the open set,
    is counted as unsupported (S7), and does not make the set incomplete."""
    f = Feed()
    f((5,), tag=0, eos=False)
    header = OcbHeaderInput(9, A, SpartnTag(TAG0, True), 1, True, 0)
    r = f.streams.observe_header(header)
    assert r.outcome == FeedOutcome.COUNTED and r.reason == FeedReason.UNSUPPORTED_SYSTEM
    assert f.streams.status(A).last_set == OcbSet(T0, True, False)
    assert f.streams.failures == {("spartn", "unsupported", 0, 9): 1}
    f.streams.observe_header(replace(header, end_of_set=False))     # no flag: the set stays as it is
    f((5,), tag=5, eos=False)
    f.streams.observe_header(replace(header, end_of_set=False))
    assert f.streams.status(A).open_set == T0 + 5
    f.streams.observe_header(replace(header, stream=StreamId(9, 9)))  # a stream never seen: nothing to end
    assert f.streams.status(A).open_set == T0 + 5


def test_rule_S9_notices_are_counted():
    """S9: a message used with data after its table, or with a DNU block of Table 8.2, is counted."""
    f = Feed()
    f.streams.count_notice(Protocol.SPARTN, "extension", 0, 0)
    f.streams.count_notice(Protocol.SPARTN, "extension", 0, 0)
    f.streams.count_notice(Protocol.RTCM3, "extension", 1019, None)
    f((5,), tag=0)
    assert f.streams.notices == {("spartn", "extension", 0, 0): 2, ("rtcm3", "extension", 1019, None): 1}
    assert f.streams.status(A).notices == f.streams.notices


def test_rule_S6_datum_change_is_a_discontinuity_for_every_satellite():
    """S6, SPARTN SF009: the datum is reported; a change within a stream ends every satellite's arc."""
    f = Feed()
    f((5, 6), tag=0, datum=0)
    assert f.streams.status(A).datum == 0
    f((5, 6), tag=5, datum=0)
    assert [f.store.discontinuities(A, Sat.gps(n)) for n in (5, 6)] == [0, 0]
    f((5,), tag=10, datum=1)
    assert f.streams.status(A).datum == 1 and f.streams.status(A).datum_changes == 1
    assert [f.store.discontinuities(A, Sat.gps(n)) for n in (5, 6)] == [1, 1]


def test_rule_S7_failures_are_counted_and_change_no_store():
    """S7: failed frames are counted per (protocol, kind, type, subtype); nothing stored changes."""
    f = Feed()
    f((5,), tag=0, eos=True)
    before = f.store.blocks(A, Sat.gps(5)), f.streams.status(A)
    r = f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 1, 2))
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 1, 2))
    f.streams.count_failure(FrameFailure(FailureKind.ENCRYPTED, Protocol.SPARTN))
    f.streams.count_failure(FrameFailure(FailureKind.DECODE, Protocol.UBX, 0x0126))
    assert r.outcome == FeedOutcome.COUNTED and r.reason == FeedReason.FRAME_FAILURE
    assert f.streams.failures == {("spartn", "decode", 1, 2): 2, ("spartn", "encrypted", None, None): 1,
                                  ("ubx", "decode", 0x0126, None): 1}
    assert f.store.blocks(A, Sat.gps(5)) == before[0]
    after = f.streams.status(A)
    assert (after.siou, after.latest, after.last_set) == (before[1].siou, before[1].latest, before[1].last_set)
    assert after.failures == f.streams.failures
