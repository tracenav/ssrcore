"""The observations of the virtual base station (DESIGN 2.6.3).

The synthetic input of examples/vrs_stream.py (seven GPS satellites and
GLONASS R01; SPARTN updates every 5 s; an ephemeris issue change, a
phase-bias arc and one station move) runs through Framer, Engine and
VirtualBase; every message goes through gnsswire's encode and decode, and
the ranges restored from the wire (RTCM 10403.3 3.5.16.3.1) are compared
with the engine's ``pseudorange_model`` and ``phase_model`` at the station
at the same epoch (DESIGN 2.6.3).
"""

from __future__ import annotations

import dataclasses
import math

import pytest
from gnsswire import rtcm3
from vrs_run import decoded_epochs, example, msm_cells, run

from ssrcore.engine import Engine
from ssrcore.model.providers import Kind, Origin
from ssrcore.model.types import (Config, GnssTime, Missing, ModelSource, Sat, Signal, SpartnSource, StreamId, System,
                                 Term, VrsConfig)
from ssrcore.vrs import msm
from ssrcore.vrs.observations import (IONOSPHERE_NOT_WITH_BIASES, MAX_CELLS, TROPOSPHERE_NOT_A_CORRECTION, Observation,
                                      Observations)
from ssrcore.vrs.placement import Station

G05, G09, G12, G17, R01 = Sat.gps(5), Sat.gps(9), Sat.gps(12), Sat.gps(17), Sat.glonass(1)
L2W = Signal.GPS_L2W


def arcs_of(result):
    """Per (sat, signal) the decoded cells by epoch, and the epochs at which the
    decoded lock time says the phase restarted (RTCM 10403.3 3.5.16.3.3)."""
    level = None
    series: dict = {}
    for t, station, messages in decoded_epochs(result):
        for key, cell in msm_cells(messages).items():
            series.setdefault(key, []).append((t, station, cell))
        for m in messages:
            if isinstance(m, rtcm3.Msm):
                level = m.number % 10
    return series, level


def loss_of_lock(p: int, n: int, a: int, b: int, dt: int, level: int) -> str:
    """[spec: RTCM 10403.3 3.5.16.3.3] the suggested loss-of-lock detection of a decoder:
    'loc', 'possible' or 'ok' from the minimum lock times p, n (ms), the supplementary
    coefficients a, b (DF407) and the interval dt (ms)."""
    if not msm.extended(level):
        if p > n or (p == n and dt >= p):
            return "loc"
        if p == n or dt <= n:
            return "ok"
        return "loc" if dt >= 2 * n - p else "possible"
    if p > n or (p == n and dt >= a):
        return "loc"
    if p == n:
        return "ok"
    if b > p:
        return "ok" if dt <= n else ("loc" if dt >= n + b - p else "possible")
    return "loc" if dt > n else "ok"


def coefficient(i: int) -> int:
    """[spec: RTCM 10403.3 Table 3.5-75] the supplementary coefficient k of a DF407 indicator."""
    return 1 if i < 64 else 2 << ((i - 64) // 32)


# ---- the exit test: through the wire and back -------------------------------------------


@pytest.mark.parametrize("level", [7, 4])
def test_decoded_ranges_equal_the_engine_models_within_the_msm_resolution(level):
    result = run(msm=level, seconds=300 if level == 7 else 40)
    res_p, res_l = msm.resolution(level)
    compared = phases = 0
    for t, station, messages in decoded_epochs(result):
        msms = [m for m in messages if isinstance(m, rtcm3.Msm)]
        assert msms and all(m.number % 10 == level for m in msms)
        for (sat, signal), cell in msm_cells(messages).items():
            p_model, l_model, wavelength, cycles = result.models[(t, sat, signal)]
            assert abs(cell.pseudorange - p_model) <= res_p / 2 + 1e-6
            compared += 1
            if l_model is None:
                assert cell.phase is None
                continue
            # Φ_s = phase_model(s) + λs·(k + N_s): the model and whole cycles, to the resolution
            assert abs(cell.phase - l_model - wavelength * cycles) <= res_l / 2 + 1e-6
            phases += 1
    assert phases == compared > 7 * 3 * 30


def test_every_cell_field_is_what_design_says():
    result = run()
    for _, _, messages in decoded_epochs(result)[:20]:
        for m in messages:
            if not isinstance(m, rtcm3.Msm):
                continue
            assert (m.clock_steering, m.external_clock, m.smoothing, m.smoothing_interval, m.iods) == (1, 0, 0, 0, 0)
            assert set(m.cnr) == {0} and set(m.half_cycle) == {0}            # DF403 0: not measured; DF420 0
        cells = msm_cells(messages)
        assert {c.channel for (s, _), c in cells.items() if s == R01} == {1}  # DF419 = K + 7
        assert {sat for sat, _ in cells} == {Sat.gps(n) for n in example().GPS_ORBITS} | {R01}
        # signals with both biases, under their own codes; L5Q only where the stream has it
        assert {sig for sat, sig in cells if sat == G05} == {Signal.GPS_L1C, L2W, Signal.GPS_L2L, Signal.GPS_L5Q}
        assert {sig for sat, sig in cells if sat == G09} == {Signal.GPS_L1C, L2W, Signal.GPS_L2L}


def test_rates_are_the_derivative_of_the_sent_phase_range():
    """DF399 + DF404 against the central difference of the decoded phase range."""
    series, _ = arcs_of(run())
    worst = 0.0
    for key, cells in series.items():
        for (t0, _, a), (t1, _, b), (t2, _, c) in zip(cells, cells[1:], cells[2:]):
            if a.phase is None or c.phase is None or c.lock < a.lock or t2 - t0 != 2.0:
                continue
            worst = max(worst, abs(b.rate - (c.phase - a.phase) / 2.0))
    print(f"largest rate difference {worst * 1000:.2f} mm/s")
    assert worst < 0.02, worst                      # m/s: the terms whose rate is not in the record


# ---- continuity -----------------------------------------------------------------------


def resets(series, level):
    """{key: [seconds of each epoch at which the decoder must take the phase as new]}:
    a first appearance, a lower lock time, or the RTCM loss-of-lock algorithm saying 'loc'."""
    out = {}
    for key, cells in series.items():
        times = []
        previous = None
        for t, _, cell in cells:
            if cell.phase is None:
                previous = None
                continue
            if previous is None:
                times.append(t)
            else:
                tp, cp = previous
                p, n = msm.lock_minimum(cp.lock, level), msm.lock_minimum(cell.lock, level)
                verdict = loss_of_lock(p, n, coefficient(cp.lock), coefficient(cell.lock),
                                       round((t - tp) * 1000), level)
                if verdict == "loc":
                    times.append(t)
            previous = (t, cell)
        out[key] = [t - example().START for t in times]
    return out


def test_phase_is_continuous_between_the_announced_resets():
    """Several minutes in which corrections update every 5 s, G05 changes its
    ephemeris issue, G09's L2W phase bias starts a new arc and the station moves:
    every phase restarts exactly at the start, at the move and (G09 L2W) at the
    first epoch after the bias block that starts a new arc (sent at 200 s,
    received at 200.5 s), and between them it never jumps."""
    s = example().Scenario()
    result = run()
    series, level = arcs_of(result)
    found = resets(series, level)
    move = s.move + 9.0                         # the 10th epoch beyond 5 km confirms (10 s at 1 Hz)
    for key, times in found.items():
        expected = [0.0, move] + ([s.bias_arc + 1.0] if key == (G09, L2W) else [])
        assert times == expected, (key, times)

    # between resets: the change of the phase is the change the rate predicts, to 1 cm
    worst = 0.0
    for key, cells in series.items():
        for (t0, _, a), (t1, _, b) in zip(cells, cells[1:]):
            if b.phase is None or a.phase is None or (t1 - example().START) in found[key]:
                continue
            worst = max(worst, abs((b.phase - a.phase) - (a.rate + b.rate) / 2.0 * (t1 - t0)))
    assert worst < 0.01, worst

    # and Φ − phase_model is the same whole number of cycles throughout each arc (the
    # wind-up does not cross ±½ cycle in these minutes, so k + N_s is constant)
    arcs: dict = {}
    for (t, sat, signal), (_, _, _, cycles) in sorted(result.models.items(), key=lambda kv: kv[0][0]):
        at = t - example().START
        if at in found[(sat, signal)]:
            arcs.setdefault((sat, signal), []).append(set())
        arcs[(sat, signal)][-1].add(cycles)
    assert all(len(arc) == 1 for runs in arcs.values() for arc in runs)
    print(f"largest unpredicted step {worst * 1000:.2f} mm")


def test_lock_time_grows_and_the_decoder_sees_no_loss_inside_an_arc():
    series, level = arcs_of(run())
    for key, cells in series.items():
        for (t0, _, a), (t1, _, b) in zip(cells, cells[1:]):
            if b.lock < a.lock:
                continue
            p, n = msm.lock_minimum(a.lock, level), msm.lock_minimum(b.lock, level)
            verdict = loss_of_lock(p, n, coefficient(a.lock), coefficient(b.lock), round((t1 - t0) * 1000), level)
            # the second epoch of an arc may read as a possible loss: its lock time was 0 at the first
            assert verdict == "ok" or (verdict == "possible" and a.lock == 0), (key, t1, a.lock, b.lock)


def test_an_ephemeris_issue_change_keeps_the_arc():
    """DESIGN 2.6.3: G05 moves to IODE 34 (same orbit, corrections naming it) without a reset."""
    result = run()
    assert result.engine.query(G05, result.epochs[-1][0], result.epochs[-1][2].position).orbit.iode == 34
    issue = example().Scenario().issue_change + 20.0
    assert issue not in sum(resets(*arcs_of(result)).values(), [])


def test_the_station_move_resets_every_lock_and_changes_the_id():
    s = example().Scenario()
    result = run()
    move = s.move + 9.0
    for t, station, messages in decoded_epochs(result):
        at = t - example().START
        ids = {m.station for m in messages}
        assert ids == {0 if at < move else 1}
        numbers = [m.number for m in messages]
        if at == move:
            assert numbers[:3] == [1005, 1033, 1230]
            m1005 = messages[0]
            assert (m1005.x, m1005.y, m1005.z) == station.arp
            assert set(sum((m.lock for m in messages if isinstance(m, rtcm3.Msm)), ())) == {0}


def test_message_order_and_rates_match_design_2_6_2():
    result = run()
    station_epochs = []
    for t, _, messages in decoded_epochs(result):
        numbers = [m.number for m in messages]
        msms = [n for n in numbers if n > 1070 and n < 1130]
        assert msms == [1077, 1087]                                  # GPS first, then GLONASS
        assert [m.multiple for m in messages if isinstance(m, rtcm3.Msm)] == [1, 0]
        if 1005 in numbers:
            assert numbers == [1005, 1033, 1230, 1077, 1087]
            station_epochs.append(t - example().START)
    gaps = [b - a for a, b in zip(station_epochs, station_epochs[1:])]
    assert station_epochs[0] == 0.0 and max(gaps) <= 10.0 and 159.0 in station_epochs


def test_the_same_input_replayed_gives_the_same_bytes():
    ex = example()
    assert ex.run(ex.Scenario(seconds=60)).wire == ex.run(ex.Scenario(seconds=60)).wire


# ---- the rules of the arcs (DESIGN 2.6.3) ------------------------------------------------


def test_a_siou_change_resets_every_arc():
    """From 50 s every message carries SIOU 4 (received from 50.1 s): every arc
    restarts at the next epoch. (A SIOU change between atmosphere messages leaves
    every satellite without atmosphere until the next HPAC of the new SIOU, rule
    A1, and so out of the observations; the change here comes with the GAD.)"""
    series, level = arcs_of(run(siou_change=50.0, seconds=80, issue_change=None, bias_arc=None, move=None))
    for key, times in resets(series, level).items():
        assert times == [0.0, 51.0], (key, times)


def test_a_fix_flag_change_resets_and_float_phases_follow_the_configuration():
    events = dict(float_from=40.0, seconds=60, issue_change=None, bias_arc=None, move=None)
    series, level = arcs_of(run(**events))
    for key, times in resets(series, level).items():
        assert times == ([0.0, 41.0] if key[0] == G12 else [0.0]), (key, times)
    # without float_phase the phase of a float bias is withheld; the pseudorange stays
    series, level = arcs_of(run(float_phase=False, **events))
    for t, _, cell in series[(G12, Signal.GPS_L1C)]:
        at = t - example().START
        assert (cell.phase is None) == (at >= 41.0) and cell.pseudorange > 0
        assert cell.lock == 0 or at < 41.0


def test_a_satellite_flagged_dnu_leaves_and_returns_with_a_new_arc():
    result = run(dnu=(40.0, 60.0), seconds=90, issue_change=None, bias_arc=None, move=None)
    series, level = arcs_of(result)
    present = [t - example().START for t, _, _ in series[(G17, Signal.GPS_L1C)]]
    # the DNU block of 40 s arrives at 40.5 s; the first block without it at 60.5 s
    assert 40.0 in present and not [a for a in present if 41.0 <= a < 61.0] and 61.0 in present
    assert resets(series, level)[(G17, Signal.GPS_L1C)] == [0.0, 61.0]


def test_a_phase_bias_discontinuity_resets_that_signal_only():
    series, level = arcs_of(run(seconds=220, issue_change=None, move=None))
    found = resets(series, level)
    assert found[(G09, L2W)] == [0.0, 201.0]
    assert all(times == [0.0] for key, times in found.items() if key != (G09, L2W))


# ---- what is left out (DESIGN 2.6.3) -----------------------------------------------------


class Altered:
    """The engine of a run, with records changed by ``change(los)`` (None to keep)."""

    def __init__(self, engine, change):
        self._engine, self._change = engine, change
        self.ephemerides, self.timebase = engine.ephemerides, engine.timebase

    def query(self, sat, t, receiver, signals=()):
        los = self._engine.query(sat, t, receiver, signals)
        return self._change(los) or los


def observe_once(change, config=VrsConfig()):
    """The observations of one epoch of the run's final state, with altered records."""
    result = run(seconds=20, issue_change=None, bias_arc=None, move=None)
    t, _, station, _ = result.epochs[-1]
    source = Observations(Altered(result.engine, change), config)
    here = Station(station.arp, 0, 1, True)
    return source.observations(t, here), source


def by_key(observations):
    return {(o.sat, o.signal): o for o in observations}


def test_terms_missing_leave_out_what_design_says():
    full = by_key(observe_once(lambda los: None)[0])
    assert len({sat for sat, _ in full}) == 8

    def drop(sat, **fields):
        return lambda los: dataclasses.replace(los, **fields) if los.sat == sat else None

    # no ionosphere, no troposphere, no tide, no corrections: the satellite is left out whole
    for fields in (dict(ionosphere=None), dict(troposphere=None), dict(orbit=None),
                   dict(modelled=dataclasses.replace(full_record_modelled(), solid_tide=None))):
        out = by_key(observe_once(drop(G09, **fields))[0])
        assert not [k for k in out if k[0] == G09] and len(out) == len(full) - 3, fields

    # a signal without its code bias is not sent; without its phase bias only the pseudorange is
    def bias(los, name, signal):
        if los.sat != G09:
            return None
        kept = {s: v for s, v in getattr(los, name).items() if s != signal}
        return dataclasses.replace(los, **{name: kept})

    out = by_key(observe_once(lambda los: bias(los, "code_bias", L2W))[0])
    assert (G09, L2W) not in out and (G09, Signal.GPS_L1C) in out
    out = by_key(observe_once(lambda los: bias(los, "phase_bias", L2W))[0])
    assert out[(G09, L2W)].phase is None and out[(G09, L2W)].lock == 0
    assert out[(G09, Signal.GPS_L1C)].phase is not None

    # no wind-up: no phase for any signal of the satellite
    out = by_key(observe_once(lambda los: dataclasses.replace(
        los, modelled=dataclasses.replace(los.modelled, windup=None)) if los.sat == G09 else None)[0])
    assert all(o.phase is None for k, o in out.items() if k[0] == G09)

    # the elevation mask
    high = by_key(observe_once(lambda los: None, VrsConfig(elevation_mask=math.radians(50.0)))[0])
    assert {sat for sat, _ in high} < {sat for sat, _ in full}


def full_record_modelled():
    result = run(seconds=20, issue_change=None, bias_arc=None, move=None)
    t, _, station, _ = result.epochs[-1]
    return result.engine.query(G09, t, station.position).modelled


# ---- the atmosphere the station may use, and origin changes (DESIGN 2.6.3, rule A9) ---------


def test_a_missing_term_is_left_out_with_the_records_reason():
    def drop(los):
        if los.sat != G09:
            return None
        missing = dict(los.missing)
        missing[(Term.TROPOSPHERE, None)] = Missing.STALE
        return dataclasses.replace(los, troposphere=None, missing=missing)

    out, source = observe_once(drop)
    assert G09 not in {k[0] for k in by_key(out)} and source.left == {G09: Missing.STALE}
    assert observe_once(lambda los: None)[1].left == {}


def test_a_model_troposphere_leaves_every_satellite_out(monkeypatch):
    """DESIGN 2.6.3: the station uses a troposphere only when it is
    a correction. An engine whose troposphere comes from UNB3m first observes nothing."""
    ex = example()
    model_first = Config(atmosphere=(ModelSource(), SpartnSource()))
    monkeypatch.setattr(ex, "Engine", lambda: Engine(model_first))
    result = ex.run(ex.Scenario(seconds=10, issue_change=None, bias_arc=None, move=None))
    t, _, station, _ = result.epochs[-1]
    los = result.engine.query(G09, t, station.position)
    assert los.corrected and los.troposphere.origin == Origin("unb3m", Kind.MODEL)
    assert los.ionosphere.origin.kind == Kind.CORRECTION and los.ionosphere.with_biases
    source = result.station.source
    assert source.last == [] and len(source.left) == 8
    assert set(source.left.values()) == {TROPOSPHERE_NOT_A_CORRECTION}
    assert not [m for m in result.epochs[-1][3] if isinstance(m, rtcm3.Msm) and m.sat_mask]


def test_an_ionosphere_not_with_the_biases_leaves_the_satellite_out():
    def other(los):
        if los.sat != G09:
            return None
        return dataclasses.replace(los, ionosphere=dataclasses.replace(
            los.ionosphere, origin=Origin("spartn:2/7", Kind.CORRECTION, StreamId(2, 7)), with_biases=False))

    out, source = observe_once(other)
    assert G09 not in {k[0] for k in by_key(out)} and source.left == {G09: IONOSPHERE_NOT_WITH_BIASES}
    out, source = observe_once(lambda los: dataclasses.replace(los, ionosphere=dataclasses.replace(
        los.ionosphere, origin=Origin("shell", Kind.MODEL))) if los.sat == G09 else None)
    assert source.left == {G09: IONOSPHERE_NOT_WITH_BIASES}


def test_an_origin_change_ends_the_arcs_of_the_satellite():
    """DESIGN 2.6.3 *Arcs* item 7: a troposphere from another provider at the next epoch."""
    result = run(seconds=20, issue_change=None, bias_arc=None, move=None)
    (t1, _, station, _), (t2, _, _, _) = result.epochs[-2], result.epochs[-1]

    def changed(los):
        if los.sat != G09 or los.t_rx != t2:
            return None
        return dataclasses.replace(los, troposphere=dataclasses.replace(
            los.troposphere, origin=Origin("sensor", Kind.CORRECTION)))

    for change, restarted in ((lambda los: None, False), (changed, True)):
        source = Observations(Altered(result.engine, change), VrsConfig())
        source.observations(t1, Station(station.arp, 0, 1, True))
        second = by_key(source.observations(t2, Station(station.arp, 0, 1, False)))
        g09 = [o for k, o in second.items() if k[0] == G09 and o.phase is not None]
        others = [o for k, o in second.items() if k[0] != G09 and o.phase is not None]
        assert g09 and others and all(o.lock > 0 for o in others)
        assert all((o.lock == 0) == restarted for o in g09)


def test_cells_beyond_64_are_split_by_satellite():
    signals = (Signal.GPS_L1C, L2W, Signal.GPS_L2L, Signal.GPS_L5Q)
    observations = [Observation(Sat.gps(n), s, 2.2e7 + 1000.0 * n + k, 2.2e7 + 1000.0 * n + k + 0.5, 7, 0.19,
                                3000, 100.0, None)
                    for n in range(1, 21) for k, s in enumerate(signals)]

    class Bare:
        class timebase:
            @staticmethod
            def leap_at(t):
                return 18

    msms = Observations(Bare(), VrsConfig()).encode(GnssTime(2295, 87400.0), observations)
    counts = [len(m.fine_range) for m in msms]
    assert counts == [64, 16] and [bin(m.sat_mask).count("1") for m in msms] == [16, 4]
    assert msms[0].sat_mask >> (64 - 16) == (1 << 16) - 1                  # PRN 1 to 16 first
    for m in msms:
        rtcm3.decode(next(rtcm3.scan(rtcm3.encode(m))))                       # each encodes and decodes


def test_the_example_runs_as_a_script(capsys):
    import runpy
    from vrs_run import EXAMPLE
    runpy.run_path(str(EXAMPLE), run_name="__main__")
    out = capsys.readouterr().out
    assert "301 epochs" in out and "1077: 301 messages" in out and "1005: 31 messages" in out
