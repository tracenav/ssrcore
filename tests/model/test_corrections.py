"""The OCB adapter and the corrections store: rules C3 to C9 (DESIGN 3.3), and S8 at block level.

Every message is built from chosen physical values, encoded to a SPARTN frame
with ``gnsswire.spartn.encode``, parsed and decoded back, and converted by
``adapters.spartn_ocb.ocb_input``. Times are GPS times; the time base (item
3a) is not needed because each test gives the message's GPS time itself,
taken here as the 32-bit tag's epoch, from ``toolkit.timescale``, plus the
tag (rules T1, T3 on the GPS scale).
"""

from __future__ import annotations

import dataclasses

import pytest
from gnsswire import spartn

from ssrcore import record
from ssrcore.adapters.spartn_ocb import dnu_table, mask_positions, ocb_input
from ssrcore.model.corrections import Biases, Corrections, Selection
from ssrcore.model.stream import Issue, Streams
from ssrcore.model.types import (FeedOutcome, FeedReason, GnssTime, MaxAge, Missing, OcbSatellite, Sat, Signal,
                                 StreamId, System)
from ssrcore.toolkit import timescale

# 2010-01-01 00:00:00 GPS time, the SPARTN TF009 epoch (week 1564, 432000 s).
EPOCH = GnssTime(*timescale.spartn_tag_to_gps(0, System.GPS))
assert EPOCH == GnssTime(1564, 432000.0)
TAG0 = 500_000_000
T0 = EPOCH + TAG0

SAT_BITS = {0: 32, 1: 24, 2: 36, 3: 37, 4: 10}   # SF011 ... SF095, size indicator 0
BIAS_BITS = {0: 6, 1: 5, 2: 8, 3: 8, 4: 6}       # SF025 ... SF107, size indicator 0


def sf020(metres: float) -> int:
    """SPARTN 6.3 inverse for SF020 (±16.382 m, 0.002 m)."""
    return round((metres + 16.382) / 0.002)


def sf029(metres: float) -> int:
    """SPARTN 6.3 inverse for SF029 (±20.46 m, 0.02 m)."""
    return round((metres + 20.46) / 0.02)


def sat(number, *, dnu=False, continuity=7, orbit=None, clock=None, phase=None, code=None):
    """A satellite block description: orbit (iode, r, a, c, yaw code or None);
    clock (SF022 code, metres, URE); phase [(position, fix, SF015 code, metres)];
    code [(position, metres)]."""
    return dict(number=number, dnu=dnu, continuity=continuity, orbit=orbit, clock=clock, phase=phase, code=code)


def message(sats, *, system=0, tag=0, solution=1, processor=0, siou=1, eos=False, yaw=False, datum=0,
            full=True) -> spartn.Message:
    mask = 0
    blocks = []
    for s in sorted(sats, key=lambda s: s["number"]):
        mask |= 1 << (SAT_BITS[system] - s["number"])
        if s["dnu"]:
            blocks.append(spartn.OcbSatellite(1, None, None, None, None, None))
            continue
        orbit = clock = bias = None
        if s["orbit"] is not None:
            iode, r, a, c, y = s["orbit"]
            orbit = spartn.Orbit(iode, sf020(r), sf020(a), sf020(c), y if yaw else None)
        if s["clock"] is not None:
            k, m, ure = s["clock"]
            clock = spartn.Clock(k, sf020(m), ure)
        if s["phase"] is not None or s["code"] is not None:
            phase = sorted(s["phase"] or ())
            code = sorted(s["code"] or ())
            width = BIAS_BITS[system]
            bias = spartn.Bias(
                0, sum(1 << (width - 1 - p) for p, *_ in phase),
                tuple(spartn.PhaseBias(f, k, sf020(v)) for _, f, k, v in phase),
                0, sum(1 << (width - 1 - p) for p, _ in code),
                tuple(sf029(v) for _, v in code))
        present = (orbit is not None) << 2 | (clock is not None) << 1 | (bias is not None)
        blocks.append(spartn.OcbSatellite(0, present, s["continuity"], orbit, clock, bias))
    body = spartn.Ocb(siou, int(eos), 0, int(yaw), datum, 0, 0, mask, tuple(blocks))
    value = TAG0 + tag if full else (TAG0 + tag) % 43200
    d = spartn.Description(system, int(full), value, solution, processor, None, None, None, None)
    return spartn.Message(0, 0, 0, d, body)


def wire(msg: spartn.Message) -> spartn.Message:
    """Through the wire and back: encode, parse, decode."""
    frame = spartn.parse(spartn.encode(msg), 0)
    assert isinstance(frame, spartn.Frame)
    return spartn.decode(frame)


def ocb(sats, **kw):
    return ocb_input(wire(message(sats, **kw)))


class Feed:
    """A corrections store and the stream store it shares, fed by tag offsets in seconds."""

    def __init__(self, stream_switch=5.0):
        self.streams = Streams(stream_switch)
        self.store = Corrections(self.streams)

    def __call__(self, sats, tag=0, **kw):
        return self.store.feed(ocb(sats, tag=tag, **kw), T0 + tag)

    def select(self, number, at, system=System.GPS, stream=StreamId(1, 0), max_age=MaxAge()):
        return self.store.select(stream, Sat(system, number), T0 + at, max_age)


G = StreamId(1, 0)
ORBIT = (40, 0.5, -0.25, 0.124, None)
CLOCK = (7, 1.234, 2)


# ---- adapter -------------------------------------------------------------------------


def test_adapter_scales_every_field_through_the_wire():
    """SPARTN 6.5.1, Table 6.2: values come back in metres, seconds and signals, satellites by mask position."""
    inp = ocb([sat(5, continuity=3, orbit=(40, 0.5, -0.25, 0.124, 10), clock=(4, 1.234, 2),
                   phase=[(0, 1, 2, 0.01), (2, 0, 7, -0.03)], code=[(0, 1.2), (1, -0.4)]),
               sat(12, dnu=True)], tag=0, siou=17, eos=True, yaw=True, datum=1, solution=9, processor=3)
    assert inp.system == System.GPS and inp.stream == StreamId(9, 3)
    assert inp.tag.full and inp.tag.value == TAG0
    assert (inp.siou, inp.end_of_set, inp.yaw_present, inp.datum) == (17, True, True, 1)
    s5, s12 = inp.satellites
    assert s5.sat == Sat.gps(5) and s12.sat == Sat.gps(12) and s12.dnu
    assert s5.continuity == 10.0                                    # SF015 code 3
    assert s5.orbit.iode == 40
    assert (s5.orbit.radial, s5.orbit.along, s5.orbit.cross) == pytest.approx((0.5, -0.25, 0.124), abs=1e-12)
    assert s5.orbit.yaw == 60.0                                     # SF021 code 10, 6 degrees
    assert s5.clock.iode_continuity == 30.0 and s5.clock.ure == 2   # SF022 code 4
    assert s5.clock.correction == pytest.approx(1.234, abs=1e-12)
    assert [(b.signal, b.fixed, b.continuity) for b in s5.phase_biases] == [
        (Signal.GPS_L1C, True, 5.0), (Signal.GPS_L2L, False, 320.0)]
    assert [b.value for b in s5.phase_biases] == pytest.approx([0.01, -0.03], abs=1e-12)
    assert [b.signal for b in s5.code_biases] == [Signal.GPS_L1C, Signal.GPS_L2W]
    assert [b.value for b in s5.code_biases] == pytest.approx([1.2, -0.4], abs=1e-12)


@pytest.mark.parametrize("system", list(System))
def test_adapter_satellite_numbers_are_mask_positions_plus_one(system):
    """SPARTN Table 6.2 SF011, SF012, SF093 to SF095: position 0 is the leftmost bit; number = position + 1."""
    inp = ocb([sat(1, clock=CLOCK), sat(SAT_BITS[system], clock=CLOCK)], system=system)
    assert inp.system == system
    assert [s.sat for s in inp.satellites] == [Sat(system, 1), Sat(system, SAT_BITS[system])]
    assert mask_positions(0b1001, 4) == (0, 3)


def test_adapter_drops_a_block_whose_value_is_outside_its_range():
    """SPARTN 6.3 and Table 6.2 SF020: code 16383 is beyond ±16.382 m; that block is treated as not sent."""
    msg = message([sat(3, orbit=ORBIT, clock=CLOCK)])
    s = msg.body.satellites[0]
    bad = spartn.Ocb(*[getattr(msg.body, f) for f in ("siou", "end_of_set", "reserved", "yaw_present", "datum",
                                                       "ephemeris_type", "mask_size", "mask")],
                     (spartn.OcbSatellite(0, s.present, s.continuity,
                                          spartn.Orbit(s.orbit.iode, 16383, s.orbit.along, s.orbit.cross, None),
                                          s.clock, None),))
    inp = ocb_input(wire(spartn.Message(0, 0, 0, msg.description, bad)))
    assert inp.satellites[0].orbit is None and inp.satellites[0].clock is not None


def test_adapter_refuses_what_is_not_ocb():
    """DESIGN 2.1: a type the adapter does not know is a programming error."""
    with pytest.raises(TypeError):
        ocb_input(dataclasses.replace(message([sat(1, clock=CLOCK)]), message_type=1))
    with pytest.raises(TypeError):
        ocb_input(dataclasses.replace(message([sat(1, clock=CLOCK)]), body=None))


def test_adapter_keeps_the_half_day_tag_as_sent():
    """SPARTN Table 7.1 TF008, TF009: the 16-bit tag is kept for the time base to resolve (T2)."""
    inp = ocb([sat(1, clock=CLOCK)], full=False)
    assert not inp.tag.full and inp.tag.value == TAG0 % 43200


# ---- store ---------------------------------------------------------------------------


def test_select_combines_orbit_clock_and_biases_of_one_issue():
    """DESIGN 3.4 step 2 with C5: the newest orbit and clock of one SIOU, biases of the same SIOU."""
    f = Feed()
    r = f([sat(5, orbit=ORBIT, clock=CLOCK, phase=[(0, 1, 7, 0.01)], code=[(0, 1.2)])], siou=3)
    assert r.outcome == FeedOutcome.STORED and r.sats == (Sat.gps(5),)
    sel = f.select(5, 1)
    assert isinstance(sel, Selection)
    assert sel.issue == Issue(3, 0) and sel.stream == G
    assert sel.orbit.value.iode == 40 and sel.orbit.time == T0
    assert sel.clock.value.correction == pytest.approx(1.234, abs=1e-12)
    assert sel.biases.value.code_bias(Signal.GPS_L1C).value == pytest.approx(1.2, abs=1e-12)
    assert sel.bias_missing is None
    assert f.store.satellites(G) == (Sat.gps(5),)


def test_rule_C3_corrections_are_constants_within_their_ages():
    """C3: no rate, no extrapolation; ages limited per kind (DESIGN 2.5, MaxAge)."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=CLOCK, code=[(0, 1.2)])])
    first, later = f.select(5, 0), f.select(5, 30)
    assert later.orbit == first.orbit and later.clock == first.clock     # same values at any age
    assert f.select(5, 30.5) == Missing.STALE                            # clock 30 s
    ages = MaxAge(orbit=10.0, clock=60.0, bias=5.0)
    assert f.select(5, 11, max_age=ages) == Missing.STALE                # orbit 10 s
    sel = f.select(5, 6, max_age=ages)
    assert sel.biases is None and sel.bias_missing == Missing.STALE      # bias 5 s
    assert f.select(5, -1) == Missing.NO_CORRECTIONS                     # nothing applies before its tag


def test_rule_C3_a_block_applies_from_its_tag_until_the_next():
    """DESIGN 2.5: an item applies from its time tag; a query between two tags uses the earlier block."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=(7, 1.0, 2))], tag=0)
    f([sat(5, orbit=ORBIT, clock=(7, 2.0, 2))], tag=10)
    assert f.select(5, 5).clock.value.correction == pytest.approx(1.0, abs=1e-12)
    assert f.select(5, 10).clock.value.correction == pytest.approx(2.0, abs=1e-12)


def test_rule_C4_clock_held_until_an_orbit_after_the_iode_change():
    """C4, SPARTN 6.2 and SF022: a clock later than the orbit by more than its IODE continuity waits."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=(4, 1.0, 2))], tag=0)                   # SF022 code 4: 30 s
    f([sat(5, clock=(4, 1.0, 2))], tag=30)
    assert isinstance(f.select(5, 30), Selection)                        # 30 s: not more than 30 s
    f([sat(5, clock=(4, 1.0, 2))], tag=31)
    assert f.select(5, 31) == Missing.IODE_CHANGE
    f([sat(5, clock=(0, 1.0, 2))], tag=32)                               # SF022 code 0: a change now
    assert f.select(5, 32) == Missing.IODE_CHANGE
    f([sat(5, orbit=(41, 0.5, -0.25, 0.124, None))], tag=32)             # orbit of the new IODE
    sel = f.select(5, 32)
    assert sel.orbit.value.iode == 41 and sel.clock.time == T0 + 32


def test_rule_C5_biases_only_from_the_same_siou():
    """C5 and S3: biases of another SIOU are missing, not taken."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=CLOCK)], siou=1)
    assert f.select(5, 0).bias_missing == Missing.NO_BIAS
    f([sat(5, code=[(0, 1.2)])], siou=2, tag=1)
    sel = f.select(5, 1)
    assert sel.issue == Issue(1, 0) and sel.biases is None and sel.bias_missing == Missing.SIOU_MISMATCH
    f([sat(5, code=[(0, 1.4)])], siou=1, tag=2)
    assert f.select(5, 2).biases.value.code_bias(Signal.GPS_L1C).value == pytest.approx(1.4, abs=1e-12)


def test_rule_C5_siou_change_mid_set_never_mixes_parts():
    """S3, C5, DESIGN 3.4 step 2: during a SIOU change the last complete issue is used."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=(7, 1.0, 2)), sat(6, orbit=ORBIT, clock=(7, 1.0, 2))], siou=1, tag=0, eos=True)
    # next epoch: the first message has the new SIOU, the second message (the rest of the set) is still due
    f([sat(5, orbit=(41, 0.6, 0, 0, None))], siou=2, tag=5)
    sel = f.select(5, 5)
    assert sel.issue == Issue(1, 0) and sel.orbit.value.iode == 40       # new orbit has no clock of SIOU 2
    f([sat(5, clock=(7, 2.0, 2)), sat(6, orbit=ORBIT, clock=(7, 2.0, 2))], siou=2, tag=5, eos=True)
    sel = f.select(5, 5)
    assert sel.issue == Issue(2, 0) and sel.orbit.value.iode == 41
    assert sel.clock.value.correction == pytest.approx(2.0, abs=1e-12)
    assert f.select(6, 5).issue == Issue(2, 0)


def test_rule_C5_orbit_and_clock_never_of_different_siou():
    """S3: an orbit of one SIOU and a clock of another give nothing."""
    f = Feed()
    f([sat(5, orbit=ORBIT)], siou=1)
    f([sat(5, clock=CLOCK)], siou=2)
    assert f.select(5, 0) == Missing.SIOU_MISMATCH


def test_rule_C6_dnu_removes_everything_until_a_block_without_it():
    """C6, SPARTN SF013: stop using the satellite immediately; a later block lifts the flag."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=CLOCK), sat(6, orbit=ORBIT, clock=CLOCK)], tag=0)
    r = f([sat(5, dnu=True)], tag=5)
    assert r.outcome == FeedOutcome.STORED and r.reason == FeedReason.DNU
    assert f.select(5, 5) == Missing.DNU and f.select(5, 1) == Missing.DNU
    assert f.store.blocks(G, Sat.gps(5)) == {} and f.store.dnu(G, Sat.gps(5)) == T0 + 5
    assert isinstance(f.select(6, 5), Selection)                         # the other satellite is untouched
    arcs = f.store.discontinuities(G, Sat.gps(5))
    r = f([sat(5, orbit=ORBIT, clock=CLOCK)], tag=4)                     # from before the flag
    assert r.outcome == FeedOutcome.REJECTED and r.reason == FeedReason.DNU
    assert f.select(5, 5) == Missing.DNU
    f([sat(5, orbit=ORBIT, clock=CLOCK)], tag=6)
    assert isinstance(f.select(5, 6), Selection) and f.store.dnu(G, Sat.gps(5)) is None
    assert arcs == 1


def test_dnu_block_in_either_layout_and_which_one_was_read():
    """SPARTN Table 6.2 SF013 (DNU alone) or Table 8.2 (DNU with SF014, SF015):
    gnsswire reads whichever fits and says which; the adapter gives a DNU block
    either way and ``dnu_table`` reports the layout."""
    plain = wire(message([sat(5, orbit=ORBIT, clock=CLOCK), sat(6, dnu=True)]))
    assert dnu_table(plain) == "6.2"
    assert dnu_table(wire(message([sat(5, orbit=ORBIT, clock=CLOCK)]))) is None
    body = plain.body
    sats = (body.satellites[0], spartn.OcbSatellite(1, 0, 7, None, None, None))
    nominal = wire(dataclasses.replace(plain, body=dataclasses.replace(body, satellites=sats)))
    assert dnu_table(nominal) == "8.2"
    assert ocb_input(nominal).satellites == ocb_input(plain).satellites
    assert ocb_input(nominal).satellites[1] == OcbSatellite(Sat.gps(6), True)


def test_rule_C7_satellite_arc_ends_after_a_gap_beyond_its_continuity():
    """C7, SPARTN 6.2 and SF015: a gap greater than the continuity time ends the arc."""
    f = Feed()
    count = lambda: f.store.discontinuities(G, Sat.gps(5))  # noqa: E731
    f([sat(5, continuity=3, clock=CLOCK)], tag=0)                        # 10 s
    assert count() == 0
    f([sat(5, continuity=3, clock=CLOCK)], tag=10)                       # gap 10 s: continuous
    f([sat(5, continuity=3, orbit=ORBIT)], tag=10)                       # same epoch, another message
    assert count() == 0
    f([sat(5, continuity=3, clock=CLOCK)], tag=21)                       # gap 11 s > 10 s
    assert count() == 1
    f([sat(5, continuity=0, clock=CLOCK)], tag=22)                       # SF015 = 0: no continuity
    assert count() == 2
    f([sat(5, continuity=7)], tag=23)                                    # no block, the indicator still counts
    assert count() == 2


def test_rule_C7_phase_bias_arcs_per_signal():
    """C7, SPARTN Table 6.12 SF015: each phase bias has its own continuity; a satellite reset ends every arc."""
    f = Feed()
    signals = lambda: (f.store.phase_discontinuities(G, Sat.gps(5), Signal.GPS_L1C),  # noqa: E731
                       f.store.phase_discontinuities(G, Sat.gps(5), Signal.GPS_L2W))
    f([sat(5, phase=[(0, 1, 2, 0.1), (1, 1, 7, 0.2)])], tag=0)           # 5 s and 320 s
    f([sat(5, phase=[(0, 1, 2, 0.1), (1, 1, 7, 0.2)])], tag=5)
    assert signals() == (0, 0)
    f([sat(5, phase=[(0, 1, 2, 0.1), (1, 1, 7, 0.2)])], tag=11)          # 6 s > 5 s for L1C only
    assert signals() == (1, 0)
    f([sat(5, continuity=0, phase=[(0, 1, 7, 0.1), (1, 1, 7, 0.2)])], tag=12)   # satellite reset
    assert signals() == (2, 1)
    f([sat(5, orbit=ORBIT, clock=CLOCK)], tag=12)
    assert f.select(5, 12).phase_discontinuities == {Signal.GPS_L1C: 2, Signal.GPS_L2W: 1}


def test_rule_C7_lost_end_of_set_and_gap_in_the_stream():
    """S5 with C7: a lost end-of-set message changes no correction; a gap beyond SF015 ends arcs."""
    f = Feed()
    f([sat(5, continuity=4, orbit=ORBIT, clock=CLOCK)], tag=0, eos=False)   # its end-of-set message lost
    f([sat(5, continuity=4, orbit=ORBIT, clock=CLOCK)], tag=5, eos=True)
    assert f.streams.status(G).last_set.time == T0 + 5 and f.store.discontinuities(G, Sat.gps(5)) == 0
    f([sat(5, continuity=4, orbit=ORBIT, clock=CLOCK)], tag=40, eos=True)    # 35 s > 30 s
    assert f.store.discontinuities(G, Sat.gps(5)) == 1


def test_rule_C8_yaw_used_when_flagged_and_valid_otherwise_zero():
    """C8, SPARTN SF008, SF021 and 8.14: yaw only with the header flag and a valid value; else zero."""
    f = Feed()
    f([sat(1, orbit=(1, 0, 0, 0, 59), clock=CLOCK), sat(2, orbit=(1, 0, 0, 0, 0x3F), clock=CLOCK),
       sat(3, orbit=(1, 0, 0, 0, 61), clock=CLOCK)], yaw=True)
    f([sat(4, orbit=(1, 0, 0, 0, None), clock=CLOCK)], yaw=False)
    one, two, three, four = (f.select(n, 0) for n in (1, 2, 3, 4))
    assert (one.yaw, one.yaw_present) == (354.0, True)
    assert one.orbit.value.yaw == 354.0
    for s in (two, three, four):                                         # invalid 0x3F, outside 0-354, absent
        assert (s.yaw, s.yaw_present, s.orbit.value.yaw) == (0.0, False, None)


def test_rule_C9_biases_by_exact_rinex_code():
    """C9, SPARTN SF025 to SF028: signals by mask position; spare positions and other tracking modes are absent."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=CLOCK, phase=[(2, 1, 7, 0.1), (5, 1, 7, 0.3)], code=[(2, 0.5), (4, 0.7)])])
    biases = f.select(5, 0).biases.value
    assert isinstance(biases, Biases)
    assert [b.signal for b in biases.phase] == [Signal.GPS_L2L]          # position 5 is spare
    assert [b.signal for b in biases.code] == [Signal.GPS_L2L]           # position 4 is spare
    assert biases.code_bias(Signal(System.GPS, "2X")) is None            # no substitution of C2X for C2L
    assert biases.phase_bias(Signal(System.GPS, "2X")) is None
    assert biases.code_bias(Signal.GPS_L2L).value == pytest.approx(0.5, abs=1e-12)


def test_galileo_e5_and_e6_biases_of_spartn_v2_0_3():
    """SPARTN v2.0.3 Table 6.2 SF102, SF105: positions 3 and 4 are L8Q/C8Q (E5) and
    L6C/C6C (E6), which v2.0.2 left spare; 5 to 14 stay spare. A v2.0.2 sender's
    block (positions 0 to 2) decodes as before."""
    positions = range(6)
    inp = ocb([sat(1, code=[(p, 0.02 * p) for p in positions],
                   phase=[(p, 1, 7, 0.01 * p) for p in positions])], system=int(System.GALILEO))
    s = inp.satellites[0]
    assert [b.signal for b in s.code_biases] == [Signal.GALILEO_L1C, Signal.GALILEO_L5Q, Signal.GALILEO_L7Q,
                                                 Signal(System.GALILEO, "8Q"), Signal(System.GALILEO, "6C")]
    assert [b.signal.phase_observation for b in s.phase_biases] == ["L1C", "L5Q", "L7Q", "L8Q", "L6C"]
    assert [b.signal.code_observation for b in s.code_biases] == ["C1C", "C5Q", "C7Q", "C8Q", "C6C"]
    assert s.code_biases[3].value == pytest.approx(0.06, abs=1e-12)          # position 5 is spare
    assert s.phase_biases[4].value == pytest.approx(0.04, abs=1e-12)
    old = ocb([sat(1, code=[(p, 0.5) for p in range(3)], phase=[(p, 1, 7, 0.1) for p in range(3)])],
              system=int(System.GALILEO))
    assert [b.signal.code for b in old.satellites[0].code_biases] == ["1C", "5Q", "7Q"]
    assert [b.signal.code for b in old.satellites[0].phase_biases] == ["1C", "5Q", "7Q"]


@pytest.mark.parametrize("system,codes", [
    (System.GLONASS, ("1C", "2C")), (System.GALILEO, ("1C", "5Q", "7Q", "8Q", "6C")),
    (System.BEIDOU, ("2I", "5P", "7I", "6I", "1P", "7P", "8P")), (System.QZSS, ("1C", "2L", "5Q"))])
def test_rule_C9_signals_of_each_system(system, codes):
    """C9, SPARTN SF026, SF028, SF102 to SF107: each system's mask order."""
    inp = ocb([sat(1, code=[(p, 0.0) for p in range(len(codes))],
                   phase=[(p, 0, 7, 0.0) for p in range(len(codes))])], system=system)
    s = inp.satellites[0]
    assert tuple(b.signal.code for b in s.code_biases) == codes
    assert tuple(b.signal.code for b in s.phase_biases) == codes


def test_rule_S8_equal_block_is_a_duplicate_different_block_a_conflict():
    """S8 at block level: the same message twice changes nothing; a different body at the same tag replaces."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=(7, 1.0, 2))], tag=0)
    r = f([sat(5, orbit=ORBIT, clock=(7, 1.0, 2))], tag=0)
    assert r.outcome == FeedOutcome.DUPLICATE
    r = f([sat(5, orbit=ORBIT, clock=(7, 1.5, 2))], tag=0)
    assert r.outcome == FeedOutcome.STORED and r.reason == FeedReason.CONFLICT
    assert f.select(5, 0).clock.value.correction == pytest.approx(1.5, abs=1e-12)
    r = f([sat(5, clock=(7, 9.0, 2))], tag=-1)                           # older than the stored clock
    assert r.outcome == FeedOutcome.SUPERSEDED and r.reason == FeedReason.OLDER


def test_blocks_older_than_the_siou_lifetime_are_dropped():
    """DESIGN 2.5, SPARTN SF005: nothing is used at 510 s, so nothing that old is kept."""
    f = Feed()
    f([sat(5, orbit=ORBIT, clock=CLOCK)], tag=0)
    f([sat(5, clock=CLOCK)], tag=300)
    f([sat(5, clock=CLOCK)], tag=600)
    (orbit, clock, _), = f.store.blocks(G, Sat.gps(5)).values()
    assert orbit == () and [b.time for b in clock] == [T0 + 300, T0 + 600]


def test_store_reasons_are_record_reasons():
    """DESIGN 2.4: the reasons the stores give are the record's own enumeration."""
    assert Missing is record.Missing


# ---- live-service facts (DESIGN 5.1) ------------------------------------------------------


def test_live_all_zero_orbit_with_the_whole_correction_in_the_clock():
    """A service may send orbit blocks of zeros and put the whole line-of-sight
    correction, computed for the user's position, in the clock field. A zero orbit
    is a block like any other (C1 moves the broadcast position by nothing), and
    the clock value is kept whole, whatever its size within SF020's range."""
    from ssrcore.toolkit import ssr

    f = Feed()
    f([sat(5, orbit=(40, 0.0, 0.0, 0.0, None), clock=(7, -14.876, 2))], tag=0)
    sel = f.select(5, 1)
    assert isinstance(sel, Selection)
    o = sel.orbit.value
    assert o.iode == 40 and max(abs(o.radial), abs(o.along), abs(o.cross)) < 1e-12
    assert sel.clock.value.correction == pytest.approx(-14.876, abs=1e-12)
    r, v = (15_600e3, 7_540e3, 20_140e3), (-1_500.0, 2_900.0, 400.0)
    assert ssr.corrected_position(r, v, (o.radial, o.along, o.cross)) == pytest.approx(r, abs=1e-9)
    # Next epoch: the provider recomputes for the user; only the clock changes.
    f([sat(5, orbit=(40, 0.0, 0.0, 0.0, None), clock=(7, -14.902, 2))], tag=5)
    assert f.select(5, 5).clock.value.correction == pytest.approx(-14.902, abs=1e-12)
    assert f.store.discontinuities(G, Sat.gps(5)) == 0


def test_live_reconnection_replays_the_last_30_s():
    """A new connection replays the last 30 s. The replayed blocks are duplicates
    (S8) or older than what is held: nothing stored changes, no arc ends, and the
    stream continues as if the link had not dropped."""
    def epoch(tag):
        return [sat(5, orbit=ORBIT, clock=(7, 1.0 + tag / 1000, 2), phase=[(0, 1, 7, 0.01)]),
                sat(7, orbit=(12, 0.1, 0.2, 0.3, None), clock=(7, -2.0 + tag / 1000, 2))]

    live, replayed = Feed(), Feed()
    for tag in range(0, 65, 5):
        live(epoch(tag), tag=tag, eos=True)
    for tag in range(0, 35, 5):
        replayed(epoch(tag), tag=tag, eos=True)
    for tag in range(5, 35, 5):                     # reconnect at 30 s: 5 to 30 s again
        r = replayed(epoch(tag), tag=tag, eos=True)
        assert r.outcome in (FeedOutcome.DUPLICATE, FeedOutcome.SUPERSEDED), (tag, r)
    for tag in range(35, 65, 5):
        replayed(epoch(tag), tag=tag, eos=True)
    for n in (5, 7):
        assert replayed.store.blocks(G, Sat.gps(n)) == live.store.blocks(G, Sat.gps(n))
        assert replayed.store.discontinuities(G, Sat.gps(n)) == live.store.discontinuities(G, Sat.gps(n)) == 0
        assert replayed.select(n, 60) == live.select(n, 60)
    assert replayed.store.phase_discontinuities(G, Sat.gps(5), Signal.GPS_L1C) == 0


# ---- corrections with rates (rule C11) -------------------------------------------------------


def test_blocks_with_rates_are_stored_as_values_and_age_from_their_tag():
    """The store keeps the rates and their t0 as part of the block (S8 compares them);
    ages, the age limits and the next change are measured from the time tag, not t0."""
    store = Corrections(Streams())
    base = ocb_input(wire(message([sat(5, orbit=(33, 0.5, -0.25, 0.124, None), clock=(7, 1.234, 2))])))
    s = base.satellites[0]
    t0 = T0 + 2.5                                                  # epoch time + ½ of a 5 s interval
    with_rates = dataclasses.replace(base, satellites=(dataclasses.replace(
        s, orbit=dataclasses.replace(s.orbit, rate=(0.001, -0.004, 0.002), reference=t0),
        clock=dataclasses.replace(s.clock, rate=0.001, acceleration=2e-5, reference=t0)),))
    assert store.feed(with_rates, T0).outcome == FeedOutcome.STORED
    found = store.select(StreamId(1, 0), Sat.gps(5), T0 + 29.0, MaxAge())
    assert isinstance(found, Selection)
    assert found.orbit.value == with_rates.satellites[0].orbit and found.orbit.time == T0
    assert found.clock.value.rate == 0.001 and found.clock.value.reference == t0
    assert store.select(StreamId(1, 0), Sat.gps(5), T0 + 31.0, MaxAge()) == Missing.STALE  # clock 30 s from the tag
    assert store.next_change(StreamId(1, 0), Sat.gps(5), T0) is None
    # The same tag without rates is another block: a conflict, not a duplicate (S8).
    assert store.feed(base, T0).outcome != FeedOutcome.DUPLICATE
