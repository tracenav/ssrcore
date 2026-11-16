"""The virtual base station behind one call per epoch (DESIGN 2.6), without
observations: a stand-in source returns fixed MSMs and records what the
station told it, so placement, moves, the station ID, what a move resets
and the message order and rates are tested through ``VirtualBase.epoch``,
and the whole epoch is encoded and decoded back with gnsswire.
"""

from __future__ import annotations

import math

import pytest
from gnsswire import rtcm3
from msm_samples import msm4

from ssrcore.model.types import GnssTime, System, VrsConfig
from ssrcore.toolkit.geodesy import enu_to_ecef, geodetic_to_ecef
from ssrcore.vrs.placement import Station
from ssrcore.vrs.station import VirtualBase

ORIGIN = (math.radians(-33.9), math.radians(18.4), 30.0)
SITE = geodetic_to_ecef(ORIGIN)


def east(metres: float):
    d = enu_to_ecef((metres, 0.0, 0.0), ORIGIN[0], ORIGIN[1])
    return (SITE[0] + d[0], SITE[1] + d[1], SITE[2] + d[2])


class Source:
    """Returns one MSM4 per system, in a scrambled order, and records every call."""

    def __init__(self, systems=(System.QZSS, System.GALILEO, System.GPS, System.BEIDOU, System.GLONASS)):
        self.systems = systems
        self.calls: list[tuple[GnssTime, Station]] = []

    def observe(self, t: GnssTime, station: Station):
        self.calls.append((t, station))
        return [msm4(s, epoch_ms=int(t.seconds * 1000) % (1 << 27)) for s in self.systems]


def epochs(start: float, count: int, week: int = 2350):
    return [GnssTime(week, start + k) for k in range(count)]


def numbers(messages):
    return [m.number for m in messages]


MSM_ORDER = [1074, 1084, 1094, 1124, 1114]


def test_nothing_before_the_first_rover_position():
    source = Source()
    station = VirtualBase(None, VrsConfig(), source)
    assert station.epoch(GnssTime(2350, 10.0)) == []
    assert station.epoch(GnssTime(2350, 11.0), None) == [] and source.calls == []
    assert station.station is None


def test_first_epoch_sends_station_messages_then_the_msms_in_order():
    source = Source()
    station = VirtualBase("engine", VrsConfig(station_id=42), source)
    out = station.epoch(GnssTime(2350, 100.0), east(0.0))
    assert numbers(out) == [1005, 1033, 1230] + MSM_ORDER
    assert all(m.station == 42 for m in out)
    assert [m.multiple for m in out[3:]] == [1, 1, 1, 1, 0]
    t, seen = source.calls[0]
    assert seen.new and seen.station_id == 42 and seen.position == station.station.position
    assert station.engine == "engine"


def test_rates_station_messages_every_ten_seconds_msms_every_epoch():
    source = Source()
    station = VirtualBase(None, VrsConfig(), source)
    with_station = []
    for t in epochs(200.0, 31):
        out = station.epoch(t, east(10.0))
        assert numbers(out)[-5:] == MSM_ORDER
        if out[0].number == 1005:
            with_station.append(t.seconds)
            assert numbers(out)[:3] == [1005, 1033, 1230]
        else:
            assert numbers(out) == MSM_ORDER
    assert with_station == [200.0, 210.0, 220.0, 230.0]
    assert [s.new for _, s in source.calls] == [True] + [False] * 30


def test_no_1230_and_indicators_follow_the_systems_sent():
    station = VirtualBase(None, VrsConfig(), Source((System.GALILEO, System.GPS)))
    out = station.epoch(GnssTime(2350, 0.0), east(0.0))
    assert numbers(out) == [1005, 1033, 1074, 1094]
    assert (out[0].gps, out[0].glonass, out[0].galileo) == (1, 0, 1)


def test_move_sends_station_messages_first_with_a_new_id_and_resets_the_source():
    """DESIGN 2.6.1: after 10 consecutive epochs beyond 5 km the station moves; that
    epoch sends 1005 and 1033 (and 1230) first with the new position and the ID
    incremented, the MSMs carry the new ID, and the source is told the placement
    is new, so it restarts every lock time (DESIGN 2.6.3, RTCM DF402)."""
    source = Source()
    config = VrsConfig(station_id=100)
    station = VirtualBase(None, config, source)
    times = epochs(1000.0, 25)
    first = station.epoch(times[0], east(0.0))
    for k, t in enumerate(times[1:4], 1):          # within 5 km
        assert numbers(station.epoch(t, east(4000.0))) == MSM_ORDER
    for t in times[4:13]:                           # 9 epochs beyond: not yet
        out = station.epoch(t, east(6000.0))
        assert numbers(out)[-5:] == MSM_ORDER and all(m.station == 100 for m in out)
        if t.seconds == 1010.0:                     # the regular 10 s station messages, old place
            assert numbers(out)[:3] == [1005, 1033, 1230] and out[0].x == first[0].x
        else:
            assert len(out) == 5
    moved = station.epoch(times[13], east(6000.0))  # the 10th: moves
    assert numbers(moved) == [1005, 1033, 1230] + MSM_ORDER
    assert all(m.station == 101 for m in moved)
    assert (moved[0].x, moved[0].y, moved[0].z) != (first[0].x, first[0].y, first[0].z)
    assert tuple(v * 0.0001 for v in (moved[0].x, moved[0].y, moved[0].z)) == station.station.position
    news = [(s.new, s.placement, s.station_id) for _, s in source.calls]
    assert news[0] == (True, 1, 100) and news[13] == (True, 2, 101)
    assert all(not new for new, _, _ in news[1:13])
    # the interval restarts at the move: next station messages 10 s later
    later = [numbers(station.epoch(t, east(6000.0)))[0] for t in times[14:24]]
    assert later.count(1005) == 1 and later.index(1005) == 9
    assert not source.calls[14][1].new


def test_pinned_station_ignores_the_rover():
    pinned = east(0.0)
    source = Source((System.GPS,))
    station = VirtualBase(None, VrsConfig(pinned_position=pinned, move_confirm=1), source)
    out = station.epoch(GnssTime(2350, 0.0))
    assert numbers(out) == [1005, 1033, 1074]
    for t in epochs(1.0, 12):
        assert numbers(station.epoch(t, east(50_000.0)))[-1] == 1074
    assert {s.placement for _, s in source.calls} == {1}


def test_every_message_of_an_epoch_decodes_back_to_what_was_given():
    """Encode the epoch with gnsswire, scan the bytes back and decode each frame."""
    station = VirtualBase(None, VrsConfig(station_id=4000), Source())
    messages = station.epoch(GnssTime(2350, 300.0), east(0.0))
    wire = b"".join(rtcm3.encode(m) for m in messages)
    decoded = [rtcm3.decode(frame, strict=True) for frame in rtcm3.scan(wire)]
    assert decoded == messages


def test_off_grid_or_repeated_epochs_are_refused():
    station = VirtualBase(None, VrsConfig(), Source())
    with pytest.raises(ValueError):
        station.epoch(GnssTime(2350, 0.5), east(0.0))
    station.epoch(GnssTime(2350, 1.0), east(0.0))
    with pytest.raises(ValueError):
        station.epoch(GnssTime(2350, 1.0), east(0.0))
