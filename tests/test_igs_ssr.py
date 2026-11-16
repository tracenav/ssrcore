"""IGS SSR (RTCM 4076): the adapter's rules IG1 to IG12 and the engine end to end.

Every message is synthetic, encoded with gnsswire from chosen values and fed
through ``adapters.framing.Framer`` with receipt times, on the scene of
tests/test_engine.py (G05 with IODE 33, R01, the receiver at −51.4°,
−162.7°). One test per mapping rule; then the engine's record against the
toolkit with the rate and the polynomial evaluated at t_tx; a satellite
stream of this format with the atmosphere of a SPARTN stream (S10, A9, A10,
A11); the IGS VTEC provider; and a SPARTN replay that IGS messages of an
undeclared stream leave unchanged.
"""

from __future__ import annotations

import dataclasses
import math

import pytest
from gnsswire import rtcm3

from ssrcore.adapters import igs_ssr
from ssrcore.adapters.framing import Framer
from ssrcore.adapters.igs_ssr import IgsReference, IgsSsr, IgsStreamId, ReferencePoint
from ssrcore.engine import Engine
from ssrcore.model.ephemeris import EphemerisStore, broadcast_state
from ssrcore.model.igs_vtec import IgsVtecIonosphere
from ssrcore.model.types import (
    CallerSource,
    Config,
    FailureKind,
    FeedOutcome,
    FeedReason,
    FrameFailure,
    GnssTime,
    LeapSeconds,
    Missing,
    OcbInput,
    Protocol,
    Sat,
    Signal,
    SpartnSource,
    StreamId,
    System,
)
from ssrcore.toolkit import geodesy, ionosphere, ssr, timescale
from ssrcore.toolkit import windup as wu
from test_engine import (
    G05, L1C, L2W, LEAP, R01, RECEIVER, STREAM, T, feed, frame_1019, frame_1020, gad, hpac, ocb,
    GPS_SAT, stream,
)

PROVIDER, SOLUTION = 300, 2
IGS = IgsStreamId(SOLUTION, PROVIDER)
EPOCH = T - 5.0                 # 87 395 s, on the 5 s grid of interval code 2
INTERVAL = 2                    # IDF004 code 2: 5 s
T0 = EPOCH + 2.5                # IGS SSR 4: epoch + ½ interval
REFERENCE = IgsReference(PROVIDER, SOLUTION, ReferencePoint.APC)

# G05's corrections as sent (IDF013 to IDF021 integers) and in metres
ORBIT = dict(radial=5000, along=-625, cross=310, dot_radial=1000, dot_along=-100, dot_cross=50)
CLOCK = dict(c0=12340, c1=500, c2=50)
D_O = (0.5, -0.25, 0.124)
RATE = (0.001, -0.0004, 0.0002)
# the engine's δC, C1, C2: the IGS values negated (IG5, Δt = Δt_b + δC_IGS / c)
C012 = (-1.234, -0.0005, -1e-6)


def igs(cls, subtype: int, sats, *, epoch: GnssTime | float = EPOCH, interval=INTERVAL, iod_ssr=1,
        provider=PROVIDER, solution=SOLUTION, **extra) -> bytes:
    """[spec: IGS SSR v1.00 8.1 to 8.4] a 4076 frame of ``cls``."""
    sow = epoch.second if isinstance(epoch, GnssTime) else int(epoch)
    values = dict(number=4076, version=0, subtype=subtype, epoch=sow, interval=interval, multiple=0,
                  iod_ssr=iod_ssr, provider=provider, solution=solution)
    names = {f.name for f in dataclasses.fields(cls)}
    if "regional" in names:
        values["regional"] = 0
    if "dispersive" in names:
        values.update(dispersive=0, mw=0)
    if "satellites" in names:
        values.update(satellite_count=len(sats), satellites=tuple(sats))
    values.update(extra)
    return rtcm3.encode(cls(**values))


def orbit_msg(iod=33, sat=5, subtype=21, **kw) -> bytes:
    return igs(rtcm3.IgsOrbit, subtype, [rtcm3.IgsOrbitSat(sat, iod, **ORBIT)], **kw)


def clock_msg(sat=5, subtype=22, **kw) -> bytes:
    return igs(rtcm3.IgsClock, subtype, [rtcm3.IgsClockSat(sat, **CLOCK)], **kw)


def code_msg(biases=((0, 123), (11, -40)), sat=5, subtype=25, **kw) -> bytes:
    return igs(rtcm3.IgsCodeBias, subtype,
               [rtcm3.IgsCodeBiasSat(sat, len(biases), tuple(rtcm3.IgsCodeBiasSignal(s, b) for s, b in biases))],
               **kw)


def phase_msg(biases=((0, 1, 2, 1234), (11, 0, 0, -567)), sat=5, yaw=64, yaw_rate=-8, subtype=26,
              **kw) -> bytes:
    """biases: (signal, integer, discontinuity, IDF028 integer)."""
    signals = tuple(rtcm3.IgsPhaseBiasSignal(s, i, 0, d, b) for s, i, d, b in biases)
    return igs(rtcm3.IgsPhaseBias, subtype, [rtcm3.IgsPhaseBiasSat(sat, len(signals), yaw, yaw_rate, signals)],
               **kw)


def scene(*frames: bytes, at: GnssTime = EPOCH + 1.0, reference=True) -> list:
    """The ephemerides, then ``frames`` received one by one from ``at``."""
    pieces = [(frame_1019(), T - 30.0), (frame_1020(), T - 29.0)]
    return pieces + [(f, at + 0.01 * k) for k, f in enumerate(frames)]


def engine(*frames: bytes, config: Config | None = None, providers=None, reference=True) -> Engine:
    e = Engine(config, providers)
    assert e.feed(LeapSeconds(LEAP)).outcome == FeedOutcome.STORED
    if reference:
        assert e.feed(REFERENCE).outcome == FeedOutcome.STORED
    feed(e, scene(*frames))
    return e


def adapted(*frames: bytes, ephemerides: EphemerisStore | None = None, adapter: IgsSsr | None = None):
    """Each frame through a fresh (or given) adapter of a declared stream; the results."""
    a = adapter or IgsSsr()
    a.declare(REFERENCE)
    out = []
    for f in frames:
        message = rtcm3.decode(rtcm3.parse(f))
        out.append(a.input(message, EPOCH + 1.0, ephemerides or EphemerisStore()))
    return out


FULL = (orbit_msg(), clock_msg(), code_msg(), phase_msg())


# ---- IG1 stream identity ---------------------------------------------------------------


def test_ig1_the_stream_is_provider_and_solution_and_never_a_spartn_stream():
    e = engine(*FULL)
    assert e.streams.streams() == (IGS,)
    assert IGS == IgsStreamId(2, 300) and IGS != StreamId(2, 300)
    assert {IGS: 1}.get(StreamId(2, 300)) is None
    assert igs_ssr.IgsStreamId.of(rtcm3.decode(rtcm3.parse(orbit_msg()))) == IGS


# ---- IG2 issue of update, sets, repeats ------------------------------------------------


def test_ig2_iod_ssr_is_the_issue_so_parts_of_two_configurations_do_not_combine():
    e = engine(orbit_msg(), clock_msg(), code_msg(iod_ssr=2), phase_msg(iod_ssr=2))
    choice = e._choice(G05, T)
    assert choice.selection is not None and choice.selection.issue.value == 1
    assert choice.selection.biases is None and choice.selection.bias_missing == Missing.SIOU_MISMATCH


def test_ig2_no_end_of_set_is_claimed_and_a_newer_epoch_ends_the_set():
    e = engine(*FULL, orbit_msg(epoch=EPOCH + 5.0))
    st = e.streams.status(IGS)
    assert st.last_set is not None and st.last_set.time == EPOCH and not st.last_set.end_of_set
    assert st.open_set == EPOCH + 5.0


def test_ig2_an_identical_repeat_is_a_duplicate_and_does_not_undo_a_completed_block():
    e = Engine()
    e.feed(LeapSeconds(LEAP)), e.feed(REFERENCE)
    results = feed(e, scene(code_msg(), phase_msg(), code_msg(), phase_msg()))[2:]
    assert [r.outcome for r in results] == [FeedOutcome.STORED, FeedOutcome.STORED,
                                            FeedOutcome.DUPLICATE, FeedOutcome.DUPLICATE]
    (bias_blocks,) = [b[2] for b in e.corrections.blocks(IGS, G05).values()]
    assert len(bias_blocks) == 1 and bias_blocks[0].value.phase and bias_blocks[0].value.code


# ---- IG3 epoch and reference time --------------------------------------------------------


def test_ig3_t0_is_epoch_plus_half_the_interval_and_the_epoch_for_interval_0():
    (_, a) = adapted(orbit_msg(), clock_msg())
    (_, b) = adapted(orbit_msg(interval=0), clock_msg(interval=0))
    assert a.item.satellites[0].orbit.reference == T0
    assert b.item.satellites[0].orbit.reference == EPOCH
    assert igs_ssr.update_interval(rtcm3.decode(rtcm3.parse(orbit_msg(interval=6)))) == 60.0


def test_ig3_a_beidou_block_is_stored_at_its_gps_time():
    bdt = EPOCH.second - 14
    e = engine(orbit_msg(sat=19, subtype=101, epoch=bdt), clock_msg(sat=19, subtype=102, epoch=bdt))
    (blocks,) = e.corrections.blocks(IGS, Sat.beidou(19)).values()
    assert blocks[0][0].time == EPOCH and blocks[1][0].time == EPOCH
    assert blocks[0][0].value.reference == T0


def test_ig3_the_epoch_is_gps_time_and_a_beidou_epoch_is_bdt():
    gps, bds = (rtcm3.decode(rtcm3.parse(f)) for f in (orbit_msg(), orbit_msg(sat=19, subtype=101)))
    assert igs_ssr.epoch_time(gps, System.GPS, T) == EPOCH
    # IDF003 "-14 s offset from BDS": the same instant is sent 14 s lower
    week, sow = timescale.bdt_to_gps(*timescale.gps_to_bdt(EPOCH.week, EPOCH.seconds))
    assert igs_ssr.epoch_time(bds, System.BEIDOU, T) == EPOCH + 14.0
    assert GnssTime(week, sow) == EPOCH
    # the week is the one nearest the receipt time, across a week boundary
    late = rtcm3.decode(rtcm3.parse(orbit_msg(epoch=604795)))
    assert igs_ssr.epoch_time(late, System.GPS, GnssTime(2296, 3.0)) == GnssTime(2295, 604795.0)


# ---- IG4 orbit, IG5 clock: units and signs -----------------------------------------------


def test_ig4_orbit_units_and_rates():
    (_, a) = adapted(orbit_msg(), clock_msg())
    o = a.item.satellites[0].orbit
    assert (o.radial, o.along, o.cross) == pytest.approx(D_O, abs=1e-12)
    assert o.rate == pytest.approx(RATE, abs=1e-15)
    assert (o.iode, o.ephemeris_type) == (33, 0)
    assert a.item.siou == 1 and a.item.datum == 0 and a.item.end_of_set is False


def test_ig4_an_orbit_of_a_new_iod_waits_for_the_clock_of_its_epoch():
    a = IgsSsr()
    first, c0 = adapted(orbit_msg(epoch=EPOCH - 5.0), clock_msg(epoch=EPOCH - 5.0), adapter=a)
    # the first orbit is held, then sent in the block of its clock
    assert isinstance(first, igs_ssr.Kept) and not first.duplicate
    assert c0.item.satellites[0].orbit.iode == 33 and c0.item.tag == igs_ssr.tag_of(EPOCH - 5.0)
    # the same IOD goes at once
    (same,) = adapted(orbit_msg(epoch=EPOCH), adapter=a)
    assert same.item.satellites[0].orbit.iode == 33
    # a new IOD waits; its clock of the same epoch is lost, the next one takes it
    new, later = adapted(orbit_msg(iod=34, epoch=EPOCH + 5.0), clock_msg(epoch=EPOCH + 10.0), adapter=a)
    assert isinstance(new, igs_ssr.Kept)
    (block,) = later.item.satellites
    assert block.orbit.iode == 34 and block.clock is not None
    assert block.orbit.reference == EPOCH + 7.5 and later.item.tag == igs_ssr.tag_of(EPOCH + 10.0)


def test_ig4_a_held_orbit_is_dropped_without_a_clock_within_its_age():
    a = IgsSsr()
    adapted(orbit_msg(epoch=EPOCH - 5.0), clock_msg(epoch=EPOCH - 5.0), adapter=a)
    (_, late) = adapted(orbit_msg(iod=34, epoch=EPOCH), clock_msg(epoch=EPOCH + 200.0), adapter=a)
    assert late.item.satellites[0].orbit is None
    assert a.left_out["orbit of a new IOD without a clock"] == 1


def test_ig4_the_old_orbit_and_clock_stay_paired_until_the_new_clock():
    # the query between the orbit of a new IOD and the clock of its epoch uses the
    # old IOD with the old clock, never the new IOD with the old clock
    old = orbit_msg(epoch=EPOCH - 5.0), clock_msg(epoch=EPOCH - 5.0)
    e = engine(*old, orbit_msg(iod=34, epoch=EPOCH))
    (blocks,) = e.corrections.blocks(IGS, G05).values()
    assert [b.value.iode for b in blocks[0]] == [33]


def test_ig5_clock_units_and_iode_continuity_from_the_orbit_interval():
    (_, c) = adapted(orbit_msg(interval=4), clock_msg(epoch=EPOCH, interval=2))
    k = c.item.satellites[0].clock
    assert (k.correction, k.rate, k.acceleration) == pytest.approx(C012, abs=1e-15)
    assert k.reference == T0
    # 87 395 s is 5 s after the start of its 15 s orbit interval (87 390 s): C4 pairs the
    # clock with an orbit of that interval or later
    assert k.iode_continuity == 87395 % 86400 % 15 == 5


def test_ig5_a_high_rate_clock_is_added_to_the_polynomial_of_its_interval():
    hr = igs(rtcm3.IgsHighRateClock, 24, [rtcm3.IgsHighRateClockSat(5, -250)], epoch=EPOCH + 2.0, interval=0)
    late = igs(rtcm3.IgsHighRateClock, 24, [rtcm3.IgsHighRateClockSat(5, 1)], epoch=EPOCH + 6.0, interval=0)
    a = IgsSsr()
    _, c, x = adapted(orbit_msg(), clock_msg(), hr, adapter=a)
    (gone,) = adapted(late, adapter=a)
    k = x.item.satellites[0].clock
    assert k.correction == pytest.approx(-(1.234 - 0.025), abs=1e-12)
    assert (k.rate, k.acceleration, k.reference) == (c.item.satellites[0].clock.rate, -2e-8 * 50, T0)
    assert x.completes
    assert isinstance(gone, igs_ssr.Refused) and a.left_out["high-rate clock without its polynomial"] == 1


# ---- IG6 biases --------------------------------------------------------------------------


def test_ig6_bias_signs_are_the_engines_and_the_fix_flag_is_idf029():
    _, p = adapted(code_msg(), phase_msg())
    sat = p.item.satellites[0]
    code = {b.signal: b.value for b in sat.code_biases}
    phase = {b.signal: (b.value, b.fixed) for b in sat.phase_biases}
    # IGS 4.3: added to the measurement to correct it, so the measurement holds −value
    assert code == pytest.approx({L1C: -1.23, L2W: 0.40})
    assert phase[L1C] == (pytest.approx(-0.1234), True) and phase[L2W] == (pytest.approx(0.0567), False)


def test_ig6_a_code_message_leaves_out_the_phase_biases_when_a_newer_phase_message_is_due():
    a = IgsSsr()
    _, c1, c2 = adapted(phase_msg(), code_msg(), code_msg(epoch=EPOCH + 5.0), adapter=a)
    assert c1.item.satellites[0].phase_biases          # same epoch: not due
    assert c2.item.satellites[0].phase_biases is None  # phase of EPOCH + 5 due


def test_ig6_a_discontinuity_counter_change_ends_the_signals_arc():
    first = phase_msg()
    same = phase_msg(epoch=EPOCH + 5.0)
    changed = phase_msg(biases=((0, 1, 3, 1234), (11, 0, 0, -567)), epoch=EPOCH + 10.0)
    e = engine(orbit_msg(), clock_msg(), code_msg(), first, same, code_msg(epoch=EPOCH + 10.0), changed)
    assert e.corrections.phase_discontinuities(IGS, G05, L1C) == 1
    assert e.corrections.phase_discontinuities(IGS, G05, L2W) == 0
    rec = e.query(G05, EPOCH + 11.0, RECEIVER, (L1C, L2W))
    assert rec.phase_bias[L1C].discontinuities == 1 and rec.phase_bias[L1C].continuity == 0.0
    assert rec.phase_bias[L2W].continuity == 10.0


# ---- IG7 signals -------------------------------------------------------------------------


def test_ig7_signal_ids_map_to_rinex_codes_and_bands_without_a_frequency_are_left_out():
    assert igs_ssr.signal(System.GPS, 2) == Signal(System.GPS, "1W")
    assert igs_ssr.signal(System.GALILEO, 9) == Signal(System.GALILEO, "7Q")
    assert igs_ssr.signal(System.BEIDOU, 13) == Signal(System.BEIDOU, "5P")
    assert igs_ssr.signal(System.QZSS, 0) == Signal(System.QZSS, "1C")
    assert igs_ssr.signal(System.GLONASS, 4) is None          # G1a: no band of the engine
    assert igs_ssr.signal(System.GPS, 9) == Signal(System.GPS, "2X")          # RTCM DF380
    assert igs_ssr.signal(System.GALILEO, 12) == Signal(System.GALILEO, "8Q")  # RTCM DF382
    assert igs_ssr.signal(System.GPS, 16) is None             # defined by neither
    assert igs_ssr.signal(System.GALILEO, 3) is None
    a = IgsSsr()
    (r,) = adapted(code_msg(biases=((0, 10), (4, 3)), sat=1, subtype=45), adapter=a)
    assert [b.signal for b in r.item.satellites[0].code_biases] == [Signal.GLONASS_L1C]
    assert a.left_out["code bias signal"] == 1


def test_ig1_satellite_ids():
    assert igs_ssr.satellite(System.GPS, 0) == Sat.gps(64)
    assert igs_ssr.satellite(System.QZSS, 1) == Sat.qzss(1)
    assert igs_ssr.satellite(System.GALILEO, 37) is None
    assert igs_ssr.satellite(System.GLONASS, 0) is None      # 64 has no GLONASS slot


# ---- IG8 issue of data -----------------------------------------------------------------


def test_ig8_the_iod_names_the_stored_issue_of_each_system():
    store = EphemerisStore()
    e = engine()
    gps = e.ephemerides.issues(G05)[0].ephemeris
    gal = dataclasses.replace(gps, sat=Sat.galileo(11), message=1046, iod=300)
    bds_toe = GnssTime(*timescale.bdt_to_gps(gps.toe.week - 1356, 7200.0 * 13))
    bds = dataclasses.replace(gps, sat=Sat.beidou(19), message=1042, iod=77, toe=bds_toe)
    for eph in (gal, bds):
        store.add(eph)
    assert igs_ssr.ephemeris_iod(Sat.galileo(11), 300 & 0xFF, store) == 300
    assert igs_ssr.ephemeris_iod(Sat.galileo(11), 12, store) == 12            # none stored: as sent
    assert igs_ssr.ephemeris_iod(Sat.beidou(19), (7200 * 13 // 720) % 240, store) == 77
    assert igs_ssr.ephemeris_iod(R01, 13, store) == 13
    assert igs_ssr.ephemeris_iod(R01, 0x80 | 13, store) is None
    assert igs_ssr.ephemeris_type(Sat.galileo(11)) == 1                      # I/NAV, IGS Table 1
    assert igs_ssr.ephemeris_type(Sat.beidou(1)) == 1 and igs_ssr.ephemeris_type(Sat.beidou(19)) == 0


def test_ig4_a_beidou_geo_orbit_is_left_out():
    a = IgsSsr()
    (r,) = adapted(orbit_msg(sat=1, subtype=101), adapter=a)
    assert isinstance(r, igs_ssr.Refused) and a.left_out["BeiDou GEO orbit"] == 1


# ---- IG9 URA ---------------------------------------------------------------------------


# IDF034 → URA bound 3^class (1 + value/4) − 1 mm → the first SF024 step (1, 2, 5, 10, 30,
# 100 cm) not below it
@pytest.mark.parametrize("ura, index", [
    (0, 0), (63, 7),                    # unknown; above 5466.5 mm
    (0b000001, 1), (0b001000, 1), (0b001100, 1),     # 0.25, 2, 5 mm
    (0b010010, 2), (0b011011, 3), (0b011111, 4),     # 12.5, 46.25, 73.25 mm
    (0b100100, 5), (0b101100, 6), (0b110000, 6),     # 161, 485, 728 mm
    (0b110111, 7),                                   # 2003.75 mm
])
def test_ig9_ura_to_the_sf024_index_whose_bound_holds(ura, index):
    assert igs_ssr.ure_index(ura) == index


def test_ig9_the_ura_goes_with_the_clocks_that_follow():
    ura = igs(rtcm3.IgsUra, 27, [rtcm3.IgsUraSat(5, 0b011011)])
    e = engine(ura, *FULL)
    rec = e.query(G05, T, RECEIVER)
    # class 3, value 3: ≤ 46.25 mm, within SF024's 5 cm
    assert rec.clock_correction.ure_index == igs_ssr.ure_index(0b011011) == 3


# ---- IG10 yaw --------------------------------------------------------------------------


def test_ig10_yaw_of_the_phase_biases_at_the_orbits_t0():
    a = IgsSsr()
    _, _, o = adapted(phase_msg(epoch=EPOCH - 5.0, yaw=64, yaw_rate=-8), orbit_msg(), clock_msg(), adapter=a)
    # 64/256 semicircle = 45°; −8/8192 semicircle/s; reference epoch + ½ interval
    expected = (45.0 - 8 * 180 / 8192 * (T0 - (EPOCH - 2.5))) % 360.0
    assert o.item.satellites[0].orbit.yaw == pytest.approx(expected)
    assert o.item.yaw_present


# ---- IG11 reference point, IG12 what is refused ------------------------------------------


def test_ig11_satellite_messages_need_a_declared_antenna_phase_centre():
    e = engine(*FULL, reference=False)
    assert e.streams.streams() == ()
    e2 = Engine()
    e2.feed(LeapSeconds(LEAP))
    e2.feed(IgsReference(PROVIDER, SOLUTION, ReferencePoint.COM))
    results = feed(e2, scene(*FULL))[2:]
    assert all(r.outcome == FeedOutcome.REJECTED and r.reason == FeedReason.NOT_AN_INPUT for r in results)
    assert "centre-of-mass" in results[0].detail


def test_ig12_sbas_is_refused_as_an_unsupported_system():
    e = engine()
    r = feed(e, [(orbit_msg(sat=1, subtype=121), T)])[0]
    assert (r.outcome, r.reason) == (FeedOutcome.REJECTED, FeedReason.UNSUPPORTED_SYSTEM)


def test_ig12_a_reserved_sub_type_is_counted_with_its_number():
    # 4076 with sub-type 202: DF002, IDF001 0, IDF002 202, then padding (IGS Table 5)
    bits = ((4076 << 11) | (0 << 8) | 202) << 1
    payload = bits.to_bytes(3, "big") + bytes(5)
    (failure,) = Framer().push("ip", rtcm3.encode_frame(payload), T)
    assert failure == FrameFailure(FailureKind.UNSUPPORTED, Protocol.RTCM3, 4076, 202, T)


def test_igs_messages_without_a_time_are_held_and_released():
    e = Engine()
    e.feed(LeapSeconds(LEAP)), e.feed(REFERENCE)
    framer = Framer()
    held = [e.feed(m) for f in FULL for m in framer.push("ip", f)]
    assert all(r.outcome == FeedOutcome.HELD and r.reason == FeedReason.NO_TIME for r in held)
    feed(e, [(frame_1019(), T - 30.0), (frame_1020(), T - 29.0)])
    assert e.timebase.held == 0
    timed = engine(*FULL)
    assert e.query(G05, T, RECEIVER) == timed.query(G05, T, RECEIVER)


# ---- the engine end to end ---------------------------------------------------------------


def test_the_record_is_the_toolkits_value_with_rates_and_polynomial_at_t_tx():
    # the phase biases of the epoch before give the orbit its yaw (IG10)
    e = engine(phase_msg(epoch=EPOCH - 5.0), *FULL)
    rec = e.query(G05, T, RECEIVER, (L1C, L2W))
    # the default atmosphere is SPARTN's, and this scene has none (A9)
    assert set(rec.missing) == {("troposphere", None), ("ionosphere", None)}, rec.missing
    assert rec.corrected and rec.quality.stream.processor_id == PROVIDER
    eph = e.ephemerides.issues(G05)[0].ephemeris
    state, clock = broadcast_state(eph, rec.t_tx)
    dt = rec.t_tx - T0
    assert 2.4 < dt < 2.45                                  # the query is 2.5 s after t0, less τ
    d_o = ssr.orbit_at(D_O, RATE, dt)
    assert d_o != pytest.approx(D_O, abs=1e-4)              # the rate moved it
    assert rec.satellite.position == ssr.corrected_position(state.position, state.velocity, d_o)
    d_c = ssr.clock_at(*C012, dt)
    assert d_c - C012[0] == pytest.approx(-(0.0005 * dt + 1e-6 * dt * dt))
    assert rec.satellite.clock == ssr.corrected_clock(clock, d_c)
    # IGS SSR 4.2 with t = t_sv − Δt_sv: the IGS δC is added to the broadcast clock offset
    igs_dc = 1.234 + 0.0005 * dt + 1e-6 * dt * dt
    assert rec.satellite.clock == pytest.approx(clock + igs_dc / 299792458.0, abs=1e-18)
    drift = e._choice(G05, T).broadcast.drift(rec.t_tx)
    assert rec.satellite.clock_drift == ssr.corrected_clock(drift, ssr.clock_rate_at(C012[1], C012[2], dt))
    assert (rec.orbit.radial, rec.orbit.along, rec.orbit.cross) == d_o
    assert rec.orbit.reference == T0 and rec.clock_correction.reference == T0
    assert rec.code_bias[L1C].value == pytest.approx(-1.23) and rec.phase_bias[L2W].value == pytest.approx(0.0567)
    # wind-up with the yaw of IG10, in the frame of SPARTN 8.14
    geodetic = geodesy.ecef_to_geodetic(RECEIVER)
    east, north, _ = geodesy.enu_basis(geodetic[0], geodetic[1])
    assert rec.orbit.yaw == pytest.approx((45.0 - 8 * 180 / 8192 * 5.0) % 360.0)
    fx, fy, _ = wu.satellite_frame(rec.satellite.position, state.velocity, math.radians(rec.orbit.yaw))
    assert rec.modelled.windup == wu.windup(rec.satellite.position, RECEIVER, east, north, fx, fy)


def test_a_correction_of_an_issue_not_yet_received_is_matched_once_it_is():
    e = engine(orbit_msg(iod=34), clock_msg())
    assert e.query(G05, T, RECEIVER).missing[("orbit", None)] == Missing.IODE_MISMATCH
    feed(e, [(frame_1019(iode=34), T - 1.0)])
    assert e.query(G05, T, RECEIVER).corrected


# ---- a satellite stream of IGS SSR with a SPARTN stream's atmosphere ------------------------


def test_igs_satellite_corrections_with_the_atmosphere_of_a_spartn_stream():
    config = Config(satellite_stream=IGS, atmosphere=(SpartnSource(STREAM),))
    e = Engine(config)
    e.feed(LeapSeconds(LEAP)), e.feed(REFERENCE)
    pieces = stream() + [(f, EPOCH + 1.0 + 0.01 * k) for k, f in enumerate(FULL)]
    feed(e, sorted(pieces, key=lambda p: p[1]))
    rec = e.query(G05, T, RECEIVER, (L1C,))
    assert rec.quality.stream.processor_id == PROVIDER                # S10: the IGS stream
    assert rec.clock_correction.reference == T0
    assert rec.troposphere is not None and rec.troposphere.origin.stream == STREAM   # A10
    io = rec.ionosphere
    assert io is not None and io.origin.stream == STREAM and io.with_biases is False  # A11
    # the SPARTN-only engine has the same atmosphere for the receiver
    spartn_only = Engine()
    spartn_only.feed(LeapSeconds(LEAP))
    feed(spartn_only, stream())
    ref = spartn_only.query(G05, T, RECEIVER, (L1C,))
    assert rec.troposphere.slant == pytest.approx(ref.troposphere.slant, abs=1e-6)


def test_a_spartn_replay_is_unchanged_by_igs_messages_it_does_not_use():
    plain = Engine()
    plain.feed(LeapSeconds(LEAP))
    feed(plain, stream())
    mixed = Engine()
    mixed.feed(LeapSeconds(LEAP))
    results = feed(mixed, stream() + [(f, T - 1.0) for f in FULL])
    assert all(r.reason == FeedReason.NOT_AN_INPUT for r in results[-len(FULL):])
    for sat_ in (G05, R01):
        assert mixed.query(sat_, T, RECEIVER) == plain.query(sat_, T, RECEIVER)


# ---- IM201 through the VTEC provider -----------------------------------------------------


def vtec_msg(cosine, sine, *, height=45, degree=2, order=2, quality=20, epoch=EPOCH, iod_ssr=1) -> bytes:
    layer = rtcm3.IgsVtecLayer(height, degree - 1, order - 1, tuple(cosine), tuple(sine))
    return igs(rtcm3.IgsVtec, 201, None, epoch=epoch, interval=6, iod_ssr=iod_ssr, quality=quality,
               layers_minus_1=0, layers=(layer,))


COS = (2000, -300, 150, 400, -120, 60)       # 0.005 TECU: 10, −1.5, 0.75, 2, −0.6, 0.3
SIN = (200, -80, 40)


def test_im201_through_the_provider_is_the_documents_slant_tec():
    provider = IgsVtecIonosphere()
    config = Config(atmosphere=(CallerSource("igs_vtec"),))
    e = engine(*FULL, vtec_msg(COS, SIN), config=config, providers={"igs_vtec": provider})
    rec = e.query(G05, T, RECEIVER, (L1C,))
    io = rec.ionosphere
    assert io is not None and io.origin.provider == "igs_vtec" and io.origin.stream == IGS
    assert io.with_biases is True                      # A11: the code biases' own stream
    # the steps of IGS 4.5.1 through the toolkit, from the record's line of sight
    geodetic = geodesy.ecef_to_geodetic(RECEIVER)
    g = rec.geometry
    enu = (math.cos(g.elevation) * math.sin(g.azimuth), math.cos(g.elevation) * math.cos(g.azimuth),
           math.sin(g.elevation))
    look = ionosphere.spherical(geodesy.geodetic_to_ecef(geodetic),
                                geodesy.enu_to_ecef(enu, geodetic[0], geodetic[1]))
    p = ionosphere.pierce(look, 450e3)
    v = ionosphere.vtec(2, 2, [c * 0.005 for c in COS], [s * 0.005 for s in SIN], p.latitude, p.longitude,
                        T.second % 86400)
    assert io.slant_tec == pytest.approx(max(v, 0.0) * p.slant_factor, rel=1e-12)
    assert io.code_delay[L1C] == ionosphere.delay(io.slant_tec, L1C.frequency())[0]
    assert io.sigma == pytest.approx(20 * 0.05 * p.slant_factor)
    assert io.age == T - EPOCH


def test_im201_negative_vtec_is_zero_and_unavailable_coefficients_give_no_value():
    provider = IgsVtecIonosphere()
    config = Config(atmosphere=(CallerSource("igs_vtec"),))
    negative = engine(*FULL, vtec_msg((-2000, 0, 0, 0, 0, 0), (0, 0, 0)), config=config,
                      providers={"igs_vtec": provider})
    assert negative.query(G05, T, RECEIVER).ionosphere.slant_tec == 0.0
    provider2 = IgsVtecIonosphere()
    missing = engine(*FULL, vtec_msg((-32768, 0, 0, 0, 0, 0), (0, 0, 0)), config=config,
                     providers={"igs_vtec": provider2})
    assert missing.query(G05, T, RECEIVER).missing[("ionosphere", None)] == Missing.NO_ATMOSPHERE


def test_im201_of_another_iod_ssr_is_not_the_satellite_streams():
    provider = IgsVtecIonosphere()
    config = Config(atmosphere=(CallerSource("igs_vtec"),))
    e = engine(*FULL, vtec_msg(COS, SIN, iod_ssr=3), config=config, providers={"igs_vtec": provider})
    assert e.query(G05, T, RECEIVER).missing[("ionosphere", None)] == Missing.SIOU_MISMATCH


def test_im201_without_a_vtec_provider_is_no_input():
    e = engine()
    r = feed(e, [(vtec_msg(COS, SIN), T)])[0]
    assert (r.outcome, r.reason) == (FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT)


def test_ig2_the_adapter_gives_typed_inputs():
    (_, a) = adapted(orbit_msg(), clock_msg())
    assert isinstance(a.item, OcbInput) and a.item.stream == IGS
