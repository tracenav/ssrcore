"""Shared types and the line-of-sight record (DESIGN 2.4, 2.5).

Exit test: the types construct, are frozen, and are hashable where they serve
as keys. Values are checked against the documents by a second route where
one exists: frequencies as multiples of the fundamental of IS-GPS-200N with IRN-IS-200N-003 and -004
3.3.1.1 and against the GLONASS ICD Table 3.1 rows; GPS times against the calendar;
the record's two sums against the geometry-free combination.
"""

from __future__ import annotations

import dataclasses
import datetime as dt
import math

import pytest

from ssrcore.model.types import (SPARTN_SIGNALS, Protocol, SPEED_OF_LIGHT, WEEK_SECONDS, AreaDefinition, CallerSource,
                                 ClockCorrection, CodeBias, Config, ModelSource, SpartnSource, FailureKind, FeedOutcome, FeedReason, FeedResult, FrameFailure,
                                 GadInput, GlonassEphemeris, GnssTime, HpacArea, HpacInput, IonosphereBlock,
                                 IonosphereSatellite, KeplerEphemeris, LeapSeconds, MaxAge, MsmType, OcbInput,
                                 OcbSatellite, OrbitCorrection, PhaseBias, RoverPosition, Sat, Signal, SpartnTag,
                                 StreamId, System, TroposphereBlock, UtcTime, VrsConfig)
from ssrcore.model.providers import Kind, Origin
from ssrcore.record import (Broadcast, ClockTerm, CodeBiasTerm, Geometry, IonosphereTerm, LineOfSight, Missing,
                            MissingTermError, Modelled, OrbitTerm, PhaseBiasTerm, SatelliteState, Term,
                            TroposphereTerm)

GPS_EPOCH = dt.datetime(1980, 1, 6)


def gps_time(when: dt.datetime) -> GnssTime:
    """A calendar instant (on the GPS scale) as week and seconds, by the calendar."""
    delta = when - GPS_EPOCH
    week, day = divmod(delta.days, 7)
    return GnssTime(week, day * 86400 + delta.seconds + delta.microseconds * 1e-6)


# ---- Sat --------------------------------------------------------------------------------


def test_sat_constructs_by_system():
    assert Sat.gps(5) == Sat(System.GPS, 5) == Sat(0, 5)
    assert Sat.glonass(24).system is System.GLONASS
    assert [s.system for s in (Sat.galileo(36), Sat.beidou(46), Sat.qzss(2))] \
        == [System.GALILEO, System.BEIDOU, System.QZSS]


def test_sat_numbers_are_those_the_spartn_masks_can_name():
    # SF011 (64), SF012 (63), SF093 to SF095 (64): position + 1.
    assert Sat.glonass(63) and Sat.gps(64)
    for bad in (lambda: Sat.gps(0), lambda: Sat.glonass(64), lambda: Sat.qzss(65)):
        with pytest.raises(ValueError):
            bad()
    with pytest.raises(TypeError):
        Sat.gps(5.0)


def test_sat_is_frozen_hashable_and_ordered():
    s = Sat.gps(5)
    with pytest.raises(dataclasses.FrozenInstanceError):
        s.number = 6
    assert {s: 1}[Sat.gps(5)] == 1
    assert sorted([Sat.galileo(1), Sat.gps(7), Sat.gps(2)]) == [Sat.gps(2), Sat.gps(7), Sat.galileo(1)]


# ---- GnssTime ---------------------------------------------------------------------------


def test_spartn_tag_epoch_by_the_calendar():
    # SPARTN TF009: the 32-bit tag counts from 2010-01-01 00:00:00. By the
    # calendar that is a Friday, 1564 weeks and 5 days after 1980-01-06.
    t = gps_time(dt.datetime(2010, 1, 1))
    assert t == GnssTime(1564, 5 * 86400.0)
    assert dt.datetime(2010, 1, 1).weekday() == 4  # Friday: day 5 of a GPS week starting Sunday


def test_week_crossover_both_ways():
    # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1: differences account for the end of
    # week.
    end = GnssTime(2000, WEEK_SECONDS - 0.5)
    start = end + 1.0
    assert start == GnssTime(2001, 0.5)
    assert start - end == 1.0 and end - start == -1.0
    assert start - 1.0 == end
    assert GnssTime(2001, 10.0) - GnssTime(2000, 604790.0) == 20.0
    assert end + 3 * WEEK_SECONDS == GnssTime(2003, WEEK_SECONDS - 0.5)
    # A remainder that rounds up to a whole week is the start of the next week.
    assert GnssTime.normalised(5, -1e-20) == GnssTime(5, 0.0)


def test_differences_keep_the_resolution_of_seconds_of_week():
    # DESIGN 4: a double second of week resolves better than 1e-10 s.
    t = GnssTime(2300, 604799.0)
    for step in (1e-9, 3.3e-10, 0.075):
        assert abs((t + step) - t - step) < 1.2e-10
    assert (GnssTime(2300, 1.0) - GnssTime(1, 1.0)) == 2299 * WEEK_SECONDS


def test_leap_second_change_against_the_calendar():
    # The leap second at the end of 2016: GPS - UTC was 17 s and became 18 s
    # at 2017-01-01 00:00:00 UTC, which is 00:00:18 GPS time.
    effective = gps_time(dt.datetime(2017, 1, 1, 0, 0, 18))
    assert effective == GnssTime(1930, 18.0)
    leap = LeapSeconds(17, 18, effective)
    assert leap.at(effective - 1e-3) == 17 and leap.at(effective) == 18
    assert LeapSeconds(18).at(effective) == 18
    with pytest.raises(ValueError):
        LeapSeconds(17, 18)


def test_gnss_time_validates_and_orders():
    for bad in ((-1, 0.0), (0, WEEK_SECONDS), (0, -0.1)):
        with pytest.raises(ValueError):
            GnssTime(*bad)
    with pytest.raises(TypeError):
        GnssTime(1.0, 0.0)
    assert GnssTime(1, 604799.0) < GnssTime(2, 0.0)
    assert hash(GnssTime(3, 1)) == hash(GnssTime(3, 1.0))
    assert GnssTime(3, 10.25) == GnssTime(3, 10, 0.25) and GnssTime(3, 10.25).seconds == 10.25


def test_gnss_time_keeps_light_travel_resolution_at_the_end_of_a_week():
    """DESIGN 4, time: a double second of week near 604800 is spaced 2^-33 s
    (35 mm of light travel); GnssTime keeps whole seconds apart, so a light
    time taken off a reception time comes back to about 1e-16 s (well under
    0.1 mm), across a week boundary too."""
    assert math.ulp(604799.9) == 2.0 ** -33
    tau = 0.0712345678901234
    for t_rx in (GnssTime(2300, 604799.97), GnssTime(2300, 604799, 0.123456789012345678)):
        t_tx = t_rx - tau
        assert abs((t_rx - t_tx) - tau) < 1e-16
        assert abs(SPEED_OF_LIGHT * ((t_rx - t_tx) - tau)) < 1e-7          # 0.1 µm
        across = t_rx + 0.95                                                 # into week 2301
        assert across.week == 2301 and abs((across - t_tx) - (tau + 0.95)) < 1e-15
    # The same in one double second of week loses up to half a spacing, 17 mm:
    rx = 604799.97
    worst = max(abs((rx - (rx - tau * k / 97)) - tau * k / 97) for k in range(1, 98))
    assert SPEED_OF_LIGHT * worst > 1e-3
    # Small steps are kept to the spacing of the fraction, 2^-53 s at most.
    t = GnssTime(2300, 604799.5)
    assert abs((t + 1e-12) - t - 1e-12) <= 2.0 ** -53


# ---- Signal -----------------------------------------------------------------------------

# IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.1.1: the nominal frequency of the SV's common source
F0 = 10.23e6


def test_gps_and_galileo_frequencies_as_multiples_of_the_fundamental():
    # L1 = 154 f0, L2 = 120 f0, L5 = 115 f0; Galileo E1 = 154 f0, E5a = 115 f0,
    # E5b = 118 f0, E6 = 125 f0; E5, 1191.795 MHz, is midway between E5a and E5b (the
    # values of Galileo OS SIS ICD 2.2 Table 2; its note: E5a and E5b are part of E5).
    assert Signal.GPS_L1C.frequency() == 154 * F0
    assert Signal.GPS_L2W.frequency() == 120 * F0 == Signal.GPS_L2L.frequency()
    assert Signal.GPS_L5Q.frequency() == 115 * F0
    assert Signal.GALILEO_L1C.frequency() == 154 * F0
    assert Signal.GALILEO_L5Q.frequency() == 115 * F0
    assert Signal.GALILEO_L7Q.frequency() == 118 * F0
    assert Signal(System.GALILEO, "6C").frequency() == 125 * F0
    assert Signal(System.GALILEO, "8X").frequency() == (115 + 118) / 2 * F0


def test_beidou_and_qzss_frequencies():
    # BDS ICD 4.2.1: B1I 1561.098 MHz = 152.6 f0; B3 1268.52 MHz = 124 f0; B1C
    # and B2a share L1 and E5a, B2b shares E5b, B2a+b is E5's centre.
    assert Signal.BEIDOU_L2I.frequency() == pytest.approx(152.6 * F0, abs=1e-3)
    assert Signal.BEIDOU_L6I.frequency() == 124 * F0
    assert Signal.BEIDOU_L1P.frequency() == Signal.GPS_L1C.frequency()
    assert Signal.BEIDOU_L5P.frequency() == Signal.GALILEO_L5Q.frequency()
    assert Signal.BEIDOU_L7I.frequency() == Signal.BEIDOU_L7P.frequency() == Signal.GALILEO_L7Q.frequency()
    assert Signal.BEIDOU_L8P.frequency() == Signal(System.GALILEO, "8Q").frequency()
    assert [Signal.QZSS_L1C.frequency(), Signal.QZSS_L2L.frequency(), Signal.QZSS_L5Q.frequency()] \
        == [154 * F0, 120 * F0, 115 * F0]


# GLONASS ICD Table 3.1, rows K = 6, 1, 0, -7 (MHz).
GLONASS_TABLE_3_1 = {6: (1605.375, 1248.625), 1: (1602.5625, 1246.4375), 0: (1602.0, 1246.0),
                     -7: (1598.0625, 1242.9375)}


def test_glonass_frequencies_against_table_3_1_and_tm23():
    for k, (l1, l2) in GLONASS_TABLE_3_1.items():
        assert Signal.GLONASS_L1C.frequency(k) == l1 * 1e6
        assert Signal.GLONASS_L2C.frequency(k) == l2 * 1e6
    for k in range(-7, 7):
        # TM-23 2.2.2: f1 = (2848 + k) 9/16 MHz, f2 = (2848 + k) 7/16 MHz; ICD 3.3.1.1: f2/f1 = 7/9.
        f1, f2 = Signal.GLONASS_L1C.frequency(k), Signal.GLONASS_L2C.frequency(k)
        assert f1 == (2848 + k) * 9 / 16 * 1e6 and f2 == (2848 + k) * 7 / 16 * 1e6
        assert f2 / f1 == pytest.approx(7 / 9, rel=1e-15)
    with pytest.raises(ValueError):
        Signal.GLONASS_L1C.frequency()
    with pytest.raises(ValueError):
        Signal.GLONASS_L1C.frequency(7)


def test_wavelengths_against_tm23():
    # TM-23 Table 4.1 lists λL1 = 0.190, λL2 = 0.244, λL5 = 0.255,
    # λG1 = 0.187 (k = 0), λE5b = 0.248 m, to three decimals.
    for signal, metres, k in ((Signal.GPS_L1C, 0.190, None), (Signal.GPS_L2W, 0.244, None),
                              (Signal.GPS_L5Q, 0.255, None), (Signal.GLONASS_L1C, 0.187, 0),
                              (Signal.GALILEO_L7Q, 0.248, None)):
        assert round(signal.wavelength(k), 3) == metres
    assert Signal.GPS_L1C.wavelength() == SPEED_OF_LIGHT / (154 * F0)


def test_spartn_bias_masks():
    # SPARTN Table 6.2 SF025 to SF028, SF102 to SF107.
    assert [Signal.from_spartn(System.GPS, k) for k in range(5)] \
        == [Signal.GPS_L1C, Signal.GPS_L2W, Signal.GPS_L2L, Signal.GPS_L5Q, None]
    assert Signal.from_spartn(System.BEIDOU, 6) == Signal.BEIDOU_L8P
    # SPARTN v2.0.3 SF102, SF105: Galileo 3 is L8Q/C8Q, 4 is L6C/C6C (spare in v2.0.2).
    assert [Signal.from_spartn(System.GALILEO, k) for k in range(3, 6)] \
        == [Signal.GALILEO_L8Q, Signal.GALILEO_L6C, None]
    assert Signal.from_spartn(System.GLONASS, 2) is None and Signal.from_spartn(System.QZSS, -1) is None
    assert Signal.BEIDOU_L6I.spartn_position == 3 and Signal(System.GPS, "2X").spartn_position is None
    assert (Signal.GPS_L2L.code_observation, Signal.GPS_L2L.phase_observation) == ("C2L", "L2L")
    # The small masks (6, 5, 8, 8, 6 bits) hold every defined signal.
    small = {System.GPS: 6, System.GLONASS: 5, System.GALILEO: 8, System.BEIDOU: 8, System.QZSS: 6}
    assert all(len(SPARTN_SIGNALS[s]) <= small[s] for s in System)


def test_signal_other_codes_are_named_but_not_substituted():
    # Rule C9: C2X is a signal of its own, with a frequency and no SPARTN bias.
    x = Signal(System.GPS, "2X")
    assert x != Signal.GPS_L2L and x.frequency() == Signal.GPS_L2L.frequency()
    for bad in ((System.GPS, "3C"), (System.GLONASS, "5Q"), (System.GPS, "1c"), (System.GPS, "1CC")):
        with pytest.raises(ValueError):
            Signal(*bad)


def test_signal_is_a_frozen_key():
    with pytest.raises(dataclasses.FrozenInstanceError):
        Signal.GPS_L1C.code = "2W"
    assert {Signal(System.GPS, "1C"): "L1"}[Signal.GPS_L1C] == "L1"


# ---- configuration ----------------------------------------------------------------------


def test_config_defaults_are_those_of_design():
    c = Config()
    assert (c.max_age.orbit, c.max_age.clock, c.max_age.bias, c.max_age.atmosphere) == (120.0, 30.0, 300.0, 120.0)
    assert c.elevation_cutoff == 0.0 and c.troposphere_height_reduction
    assert (c.stream_switch, c.hold_limit) == (5.0, 256)
    v = c.vrs
    assert (v.move_distance, v.move_confirm, v.rate, v.station_interval, v.msm, v.station_id) \
        == (5000.0, 10.0, 1.0, 10.0, MsmType.MSM7, 0)
    assert v.elevation_mask == pytest.approx(math.radians(5.0)) and v.float_phase and v.pinned_position is None
    assert hash(c) == hash(Config())
    with pytest.raises(dataclasses.FrozenInstanceError):
        c.hold_limit = 1


def test_config_refuses_values_the_specifications_exclude():
    for bad in (lambda: MaxAge(orbit=600.0), lambda: MaxAge(atmosphere=450.0), lambda: VrsConfig(rate=5.0),
                lambda: VrsConfig(station_id=4096), lambda: VrsConfig(station_interval=20.0),
                lambda: VrsConfig(msm=5), lambda: Config(hold_limit=0)):
        with pytest.raises(ValueError):
            bad()
    assert VrsConfig(pinned_position=(1, 2, 3)).pinned_position == (1.0, 2.0, 3.0)


def test_config_names_the_satellite_stream_and_the_atmosphere_sources():
    """DESIGN 2.5, rules S10 and A9: the defaults are the engine before the seam."""
    c = Config()
    assert c.satellite_stream is None and c.atmosphere == (SpartnSource(),)
    a, b = StreamId(1, 0), StreamId(2, 7)
    two = Config(satellite_stream=a, atmosphere=[SpartnSource(), SpartnSource(b), ModelSource(),
                                                 CallerSource("sensor")])
    assert two.atmosphere == (SpartnSource(), SpartnSource(b), ModelSource("unb3m"), CallerSource("sensor"))
    assert hash(two) == hash(Config(satellite_stream=a, atmosphere=two.atmosphere))
    assert Config(atmosphere=()).atmosphere == ()
    for bad, error in ((lambda: Config(satellite_stream=(1, 0)), TypeError),
                       (lambda: Config(atmosphere=("spartn",)), TypeError),
                       (lambda: Config(atmosphere=(CallerSource("x"), CallerSource("x"))), ValueError),
                       (lambda: ModelSource("saastamoinen"), ValueError),
                       (lambda: CallerSource(""), ValueError),
                       (lambda: SpartnSource((2, 7)), TypeError)):
        with pytest.raises(error):
            bad()


def test_corrections_with_rates_carry_their_reference_time():
    """Rule C11: rates and polynomial terms come with t0; SPARTN's carry none."""
    t0 = GnssTime(2300, 100.0)
    o = OrbitCorrection(40, 0.1, -0.2, 0.05, rate=[0.001, -0.004, 0], reference=t0)
    assert o.rate == (0.001, -0.004, 0.0) and o.reference == t0
    assert OrbitCorrection(40, 0.1, -0.2, 0.05).rate is None
    c = ClockCorrection(0.0, 0.3, 2, rate=0.001, acceleration=None, reference=t0)
    assert c.polynomial and not ClockCorrection(320.0, 0.3, 2).polynomial
    for bad, error in ((lambda: OrbitCorrection(40, 0.1, 0.0, 0.0, rate=(0.1, 0.2)), ValueError),
                       (lambda: OrbitCorrection(40, 0.1, 0.0, 0.0, rate=(0.1, 0.2, 0.3)), ValueError),
                       (lambda: OrbitCorrection(40, 0.1, 0.0, 0.0, reference=t0), ValueError),
                       (lambda: OrbitCorrection(40, 0.1, 0.0, 0.0, rate=(0, 0, 0), reference=100.0), TypeError),
                       (lambda: ClockCorrection(0.0, 0.3, 2, acceleration=1e-5), ValueError),
                       (lambda: ClockCorrection(0.0, 0.3, 2, reference=t0), ValueError)):
        with pytest.raises(error):
            bad()


# ---- inputs -----------------------------------------------------------------------------


def ocb() -> OcbInput:
    return OcbInput(
        system=System.GPS, stream=StreamId(1, 2), tag=SpartnTag(12345, full=False), siou=17,
        end_of_set=True, yaw_present=False, datum=0, ephemeris_type=0,
        satellites=(
            OcbSatellite(Sat.gps(5), dnu=False, continuity=60.0,
                         orbit=OrbitCorrection(iode=40, radial=0.1, along=-0.2, cross=0.05),
                         clock=ClockCorrection(iode_continuity=320.0, correction=0.3, ure=2),
                         phase_biases=(PhaseBias(Signal.GPS_L1C, 0.01, True, 120.0),),
                         code_biases=(CodeBias(Signal.GPS_L1C, -0.4), CodeBias(Signal.GPS_L2W, 0.2))),
            OcbSatellite(Sat.gps(9), dnu=True),
        ))


def test_inputs_are_frozen_and_comparable_for_duplicates():
    # Rule S8 compares whole messages: equal content, equal hash.
    a, b = ocb(), ocb()
    assert a == b and hash(a) == hash(b) and len({a, b}) == 1
    with pytest.raises(dataclasses.FrozenInstanceError):
        a.siou = 3
    assert dataclasses.replace(a, siou=18) != a


def test_dnu_block_carries_nothing_else():
    # SPARTN SF013: with DNU set the remaining fields are not present.
    with pytest.raises(ValueError):
        OcbSatellite(Sat.gps(9), dnu=True, continuity=5.0)


def test_spartn_time_tags():
    # TF009: 16-bit tags are seconds of the half day; 32-bit tags fill 32 bits.
    assert SpartnTag(43199, full=False) and SpartnTag(2 ** 32 - 1, full=True)
    with pytest.raises(ValueError):
        SpartnTag(43200, full=False)


def test_atmosphere_inputs_construct():
    hpac = HpacInput(System.GALILEO, StreamId(1, 2), SpartnTag(500, False), siou=17, aiou=3, areas=(
        HpacArea(area_id=7, grid_points=4,
                 troposphere=TroposphereBlock(1, 3, 2.31, (0.12, 0.001, -0.002), (0.01, None, 0.0, -0.01)),
                 ionosphere=IonosphereBlock(0, (IonosphereSatellite(Sat.galileo(11), 6, (12.5,), None),))),))
    gad = GadInput(StreamId(1, 2), SpartnTag(500, False), siou=17, aiou=3, areas=(
        AreaDefinition(7, latitude=52.0, longitude=4.0, latitude_nodes=2, longitude_nodes=2,
                       latitude_spacing=2.5, longitude_spacing=2.5),))
    assert hash(hpac) and hash(gad)


def test_ephemerides_construct_for_their_systems():
    t = GnssTime(2300, 7200.0)
    k = KeplerEphemeris(Sat.gps(5), 1019, 40, t, t, 5153.6, 0.01, 0.1, 1e-9, 0.2, 0.3, -2.6e-9, 0.3, 1e-10,
                        1e-6, 1e-6, 200.0, 20.0, 1e-7, 1e-7, 1e-4, 1e-12, 0.0, (5e-9,), 0, 1)
    g = GlonassEphemeris(Sat.glonass(3), 5, 30, t, (1e7, 2e7, 1e6), (100.0, -200.0, 3000.0), (0.0, 0.0, 1e-6),
                         1e-5, 1e-12, 2e-9, 0, 1800.0)
    assert hash(k) and hash(g)
    with pytest.raises(ValueError):
        dataclasses.replace(k, sat=Sat.glonass(3))
    with pytest.raises(ValueError):
        dataclasses.replace(g, frequency_channel=7)


def test_other_inputs_and_results():
    assert UtcTime(2024, 2, 29, 23, 59, 59.5) and RoverPosition(0.9, 0.1, 45.0, 3600.0)
    f = FrameFailure("decode", "spartn", 0, 3)
    assert f.kind is FailureKind.DECODE and f.protocol is Protocol.SPARTN and hash(f)
    with pytest.raises(ValueError):
        FrameFailure("decode", "can")
    r = FeedResult("stored", sats=(Sat.gps(5),))
    assert r.outcome is FeedOutcome.STORED and r.reason is None
    assert FeedResult(FeedOutcome.HELD, FeedReason.NO_LEAP_SECONDS).reason is FeedReason.NO_LEAP_SECONDS
    with pytest.raises(ValueError):
        FeedResult(FeedOutcome.REJECTED)


# ---- the record -------------------------------------------------------------------------

T_RX = GnssTime(2300, 100000.0)
T_TX = T_RX - 0.075
L1, L2 = Signal.GPS_L1C, Signal.GPS_L2W


def record(**changes) -> LineOfSight:
    sat = SatelliteState(position=(2e7, 1e7, 1.5e7), velocity=(1000.0, -2000.0, 500.0), clock=1.5e-4,
                         clock_drift=1e-11, broadcast=Broadcast((2e7, 1e7, 1.5e7), 1.5e-4, 40, 1019, T_RX - 3600.0,
                                                               3600.0))
    tec = 20.0
    fields = dict(
        sat=Sat.gps(5), t_rx=T_RX, receiver=(3.9e6, 3e5, 5e6), signals=(L1, L2), t_tx=T_TX,
        geometry=Geometry(22_000_000.0, -150.0, (0.6, 0.48, 0.64), 1.0, 0.7),
        satellite=sat,
        orbit=OrbitTerm(0.1, -0.2, 0.05, (0.01, 0.02, 0.03), 40, None, T_RX - 5.0, 5.0),
        clock_correction=ClockTerm(0.3, 320.0, 2, 0.02, T_RX - 2.0, 2.0),
        modelled=Modelled(earth_rotation=12.5, shapiro=0.018, relativity_clock=-2e-8, solid_tide=-0.07,
                          tide_displacement=(0.01, 0.02, 0.1), windup=0.2),
        code_bias={L1: CodeBiasTerm(-0.4, T_RX - 2.0, 2.0), L2: CodeBiasTerm(0.2, T_RX - 2.0, 2.0)},
        phase_bias={L1: PhaseBiasTerm(0.01, True, 120.0, 0, T_RX - 2.0, 2.0)},
        troposphere=TroposphereTerm(4.6, 2.3, 2.0, 2.01, 0.05, 0.002, 3, 7, 10.0),
        ionosphere=IonosphereTerm(tec, tec - 0.5, 0.5,
                                  {s: 40.3e16 * tec / s.frequency() ** 2 for s in (L1, L2)}, 6, 7, 10.0),
        frequency={L1: L1.frequency(), L2: L2.frequency()},
        missing={(Term.PHASE_BIAS, L2): Missing.NO_BIAS},
    )
    fields.update(changes)
    return LineOfSight(**fields)


def test_the_seams_fields_default_to_none_and_a_models_terms_make_the_sums():
    """DESIGN 2.4 and 3.6: every field the providers added has a default,
    so records built as before are unchanged; a model's troposphere and a caller's
    ionosphere, without the SPARTN parts, still make the sums."""
    r = record()
    assert (r.troposphere.wet_zenith, r.troposphere.origin, r.troposphere.sigma) == (None, None, None)
    assert (r.ionosphere.origin, r.ionosphere.sigma, r.ionosphere.with_biases) == (None, None, None)
    assert (r.orbit.rate, r.orbit.reference) == (None, None)
    assert (r.clock_correction.rate, r.clock_correction.acceleration, r.clock_correction.reference) == (None,) * 3
    model = TroposphereTerm(r.troposphere.slant, 2.3, 2.0, 2.01, None, None, None, None, None, 0.1,
                            Origin("unb3m", Kind.MODEL), 0.049)
    shell = IonosphereTerm(r.ionosphere.slant_tec, None, None, r.ionosphere.code_delay, None, None, 3.0,
                           Origin("shell", Kind.CORRECTION), 2.0, False)
    m = record(troposphere=model, ionosphere=shell)
    assert m.pseudorange_model(L1) == r.pseudorange_model(L1) and m.phase_model(L1) == r.phase_model(L1)


def test_pseudorange_model_is_the_sum_of_design_2_4():
    r = record()
    m = r.modelled
    expected = (r.geometry.range + m.earth_rotation + m.shapiro
                - SPEED_OF_LIGHT * (r.satellite.clock + m.relativity_clock)
                + r.troposphere.slant + r.ionosphere.code_delay[L1] + r.code_bias[L1].value + m.solid_tide)
    assert r.pseudorange_model(L1) == pytest.approx(expected, abs=1e-9)
    assert r.available


def test_code_minus_phase_is_twice_the_ionosphere_plus_biases_and_windup():
    # Independent route: in P - Φ every non-dispersive term cancels, leaving
    # 2·I + b − B − λ·w (SPARTN 8.11: the phase delay is the code delay negated).
    r = record()
    lam = SPEED_OF_LIGHT / (154 * F0)
    difference = r.pseudorange_model(L1) - r.phase_model(L1)
    expected = 2 * r.ionosphere.code_delay[L1] + r.code_bias[L1].value - r.phase_bias[L1].value - lam * 0.2
    assert difference == pytest.approx(expected, abs=1e-6)
    assert r.ionosphere.phase_delay(L1) == -r.ionosphere.code_delay[L1]
    # TEC to metres: 1 TECU on L1 is 0.162 m of code delay.
    assert 40.3e16 / L1.frequency() ** 2 == pytest.approx(0.1624, abs=1e-4)


def test_a_missing_term_fails_the_sum_and_says_why():
    r = record()
    with pytest.raises(MissingTermError) as e:
        r.phase_model(L2)
    assert (e.value.term, e.value.signal, e.value.reason) == (Term.PHASE_BIAS, L2, Missing.NO_BIAS)
    no_tide = record(modelled=dataclasses.replace(r.modelled, solid_tide=None),
                     missing={(Term.SOLID_TIDE, None): Missing.NO_EPHEMERIS})
    with pytest.raises(MissingTermError, match="solid_tide"):
        no_tide.pseudorange_model(L1)
    bare = LineOfSight(Sat.gps(5), T_RX, (0.0, 0.0, 0.0), missing={("satellite", None): "no_ephemeris"})
    assert not bare.available and bare.why(Term.SATELLITE) is Missing.NO_EPHEMERIS
    with pytest.raises(MissingTermError):
        bare.pseudorange_model(L1)


def test_record_is_immutable_throughout():
    r = record()
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.t_tx = T_RX
    with pytest.raises(TypeError):
        r.code_bias[L2] = CodeBiasTerm(0.0, T_RX, 0.0)
    with pytest.raises(TypeError):
        r.ionosphere.code_delay[L1] = 0.0
    with pytest.raises(ValueError):
        record(missing={(Term.ORBIT, None): "not a reason"})
