"""Messages and rates of the virtual base station (DESIGN 2.6.2), and the
station messages 1005, 1033 and 1230 (DESIGN 2.6, 2.6.4).

Every station message is encoded with gnsswire, framed, parsed and decoded
back, and compared with the values it was given; the 1005 is also checked
bit for bit against the worked example of RTCM 10403.3 section 4.2.
"""

from __future__ import annotations

import dataclasses

import pytest
from gnsswire import rtcm3
from msm_samples import msm4

from ssrcore import __version__
from ssrcore.model.types import GnssTime, System, VrsConfig
from ssrcore.vrs.placement import Station
from ssrcore.vrs.schedule import (Schedule, message_1005, message_1033, message_1230, msm_system, order_msms,
                                  station_messages)

ALL = (System.GPS, System.GLONASS, System.GALILEO, System.BEIDOU, System.QZSS)
STATION = Station((40_278_930_000, 3_070_450_000, 49_194_740_000), 1234, 1, True)

# RTCM 10403.3 section 4.2, the Hex-ASCII example of a 1005.
EXAMPLE_1005 = bytes.fromhex("D300133ED7D30202980EDEEF34B4BD62AC0941986F33360B98")


def round_trip(message):
    """Encode with gnsswire, parse the frame back and decode it."""
    wire = rtcm3.encode(message)
    frame = rtcm3.parse(wire, 0)
    assert isinstance(frame, rtcm3.Frame) and frame.length == len(wire)
    return rtcm3.decode(frame, strict=True)


# ---- the epoch grid ---------------------------------------------------------------------


@pytest.mark.parametrize("rate, on, off", [
    (1.0, [(2300, 0.0), (2300, 604799.0), (2301, 17.0)], [(2300, 0.5), (2300, 17.000002)]),
    (2.0, [(2300, 0.5), (2300, 10.0)], [(2300, 0.25)]),
    (0.5, [(2300, 2.0), (2300, 604798.0)], [(2300, 1.0), (2300, 604799.0)]),
    (1.5, [(2300, 2.0), (2300, 0.0) ], [(2300, 1.0)]),
])
def test_epochs_lie_on_whole_multiples_of_the_period_in_gps_time(rate, on, off):
    """DESIGN 2.6.2; RTCM 3.1.4 for the 1 µs. The grid runs on across weeks
    (604800 is a multiple of 2 s), and at 1.5 Hz every 2 s is on it."""
    schedule = Schedule(VrsConfig(rate=rate))
    assert all(schedule.is_epoch(GnssTime(w, s)) for w, s in on)
    assert not any(schedule.is_epoch(GnssTime(w, s)) for w, s in off)
    assert schedule.is_epoch(GnssTime(2300, 100.0, 0.4e-6) if rate != 1.5 else GnssTime(2300, 100.0))


def test_next_epoch():
    schedule = Schedule(VrsConfig(rate=2.0))
    assert schedule.next_epoch(GnssTime(2300, 10.2)) == GnssTime(2300, 10.5)
    assert schedule.next_epoch(GnssTime(2300, 10.5)) == GnssTime(2300, 10.5)
    assert schedule.next_epoch(GnssTime(2300, 604799.9)) == GnssTime(2301, 0.0)
    at_one_and_half = Schedule(VrsConfig(rate=1.5))
    assert at_one_and_half.is_epoch(at_one_and_half.next_epoch(GnssTime(2300, 100.1)))


def test_epochs_must_be_on_the_grid_and_increase():
    schedule = Schedule(VrsConfig())
    schedule.advance(GnssTime(2300, 10.0))
    with pytest.raises(ValueError):
        schedule.advance(GnssTime(2300, 10.0))
    with pytest.raises(ValueError):
        schedule.advance(GnssTime(2300, 9.0))
    with pytest.raises(ValueError):
        schedule.advance(GnssTime(2300, 11.5))
    with pytest.raises(ValueError):
        schedule.station_due(GnssTime(2300, 12.0), False)   # not advanced to


# ---- rates ------------------------------------------------------------------------------


def due_times(config: VrsConfig, start: float, seconds: float, new_at=()):
    schedule = Schedule(config)
    t = schedule.next_epoch(GnssTime(2300, start))
    out = []
    while t - GnssTime(2300, start) < seconds:
        schedule.advance(t)
        if schedule.station_due(t, (t - GnssTime(2300, 0.0)) in new_at):
            out.append(t - GnssTime(2300, 0.0))
        t = t + schedule.period
    return out


def test_station_messages_every_ten_seconds_at_one_hertz():
    """DESIGN 2.6.2: 1005, 1033 (and 1230) every station_interval, 10 s by default,
    from the first epoch."""
    assert due_times(VrsConfig(), 100.0, 60.0) == [100.0, 110.0, 120.0, 130.0, 140.0, 150.0]


def test_station_messages_at_two_hertz_and_on_every_move():
    assert due_times(VrsConfig(rate=2.0), 0.0, 25.0) == [0.0, 10.0, 20.0]
    # a move at 13 s sends them at once and restarts the interval
    assert due_times(VrsConfig(), 0.0, 30.0, new_at={13.0}) == [0.0, 10.0, 13.0, 23.0]


def test_no_gap_exceeds_the_interval_when_it_is_not_a_multiple_of_the_period():
    """VrsConfig allows 15 s at 0.5 Hz: the messages go out every 14 s, never 16."""
    times = due_times(VrsConfig(rate=0.5, station_interval=15.0), 0.0, 100.0)
    gaps = [b - a for a, b in zip(times, times[1:])]
    assert gaps and max(gaps) <= 15.0 and gaps == [14.0] * len(gaps)


# ---- order --------------------------------------------------------------------------------


def test_msms_are_ordered_gps_glonass_galileo_beidou_qzss_with_the_multiple_message_bit():
    """DESIGN 2.6.2 item 3; RTCM DF393: 1 on every MSM of the epoch but the last."""
    given = [msm4(System.QZSS), msm4(System.BEIDOU), msm4(System.GLONASS), msm4(System.GPS),
             msm4(System.GALILEO, multiple=0)]
    ordered = order_msms(given, 321)
    assert [msm_system(m) for m in ordered] == list(ALL)
    assert [m.multiple for m in ordered] == [1, 1, 1, 1, 0]
    assert all(m.station == 321 for m in ordered)
    single = order_msms([msm4(System.GALILEO)], 5)
    assert single[0].multiple == 0
    assert order_msms([], 5) == []


def test_msms_of_one_system_keep_their_order():
    first, second = msm4(System.GPS, satellite=3), msm4(System.GPS, satellite=9)
    ordered = order_msms([second, msm4(System.GLONASS), first], 1)
    assert [m.sat_mask for m in ordered[:2]] == [second.sat_mask, first.sat_mask]


def test_sbas_and_other_messages_are_refused():
    sbas = dataclasses.replace(msm4(System.GPS), number=1104)
    with pytest.raises(ValueError):
        order_msms([sbas], 1)
    with pytest.raises(ValueError):
        msm_system(dataclasses.replace(msm4(System.GPS), number=1070))


# ---- station messages -------------------------------------------------------------------


def test_station_messages_order_and_1230_only_with_glonass():
    """DESIGN 2.6.2 items 1 and 2; RTCM Table 3.1-2, 1230 "if GLONASS data is transmitted"."""
    assert [m.number for m in station_messages(STATION, ALL)] == [1005, 1033, 1230]
    assert [m.number for m in station_messages(STATION, {System.GPS, System.GALILEO})] == [1005, 1033]
    assert all(m.station == 1234 for m in station_messages(STATION, ALL))


@pytest.mark.parametrize("systems", [set(), {System.GPS}, {System.GLONASS, System.BEIDOU}, set(ALL)])
def test_1005_decodes_to_the_values_given(systems):
    """RTCM Table 3.5-6: ARP, ID, the systems sent, non-physical (DF141 = 1),
    single oscillator (DF142 = 1), DF021 and DF364 zero."""
    decoded = round_trip(message_1005(STATION, systems))
    assert decoded.station == STATION.station_id
    assert (decoded.x, decoded.y, decoded.z) == STATION.arp
    assert tuple(v * 0.0001 for v in (decoded.x, decoded.y, decoded.z)) == STATION.position
    assert (decoded.gps, decoded.glonass, decoded.galileo) == (
        int(System.GPS in systems), int(System.GLONASS in systems), int(System.GALILEO in systems))
    assert (decoded.reference_station, decoded.single_oscillator) == (1, 1)
    assert (decoded.itrf_year, decoded.quarter_cycle, decoded.reserved) == (0, 0, 0)


def test_1005_matches_the_worked_example_of_rtcm_4_2():
    """RTCM 10403.3 4.2: station 2003, GPS only, ARP (1114104.5999, −4850729.7108,
    3975521.4643) m. The example is a physical station with DF142 = 0; with
    those two indicators as in the example, the message encodes to its bytes."""
    station = Station((11141045999, -48507297108, 39755214643), 2003, 1, True)
    ours = message_1005(station, {System.GPS})
    as_example = dataclasses.replace(ours, reference_station=0, single_oscillator=0)
    assert rtcm3.encode(as_example) == EXAMPLE_1005
    assert rtcm3.decode(rtcm3.parse(EXAMPLE_1005, 0)) == as_example


def test_1033_decodes_to_the_values_given():
    """RTCM Table 3.5-36; 3.1.6 ADVNULLANTENNA; DF031 0; unknown serials have zero counters."""
    decoded = round_trip(message_1033(STATION))
    assert decoded.station == 1234
    assert (decoded.antenna_count, decoded.antenna) == (14, b"ADVNULLANTENNA")
    assert decoded.antenna_setup == 0
    assert (decoded.antenna_serial_count, decoded.antenna_serial) == (0, b"")
    assert decoded.receiver == b"ssrcore" and decoded.receiver_count == 7
    assert decoded.firmware == __version__.encode() and decoded.firmware_count == len(__version__)
    assert (decoded.receiver_serial_count, decoded.receiver_serial) == (0, b"")
    assert all(n <= 31 for n in (decoded.antenna_count, decoded.receiver_count, decoded.firmware_count))


def test_1230_decodes_to_the_values_given():
    """RTCM Table 3.5-109: aligned (DF421 = 1), L1 C/A and L2 C/A present and zero."""
    decoded = round_trip(message_1230(STATION))
    assert decoded.station == 1234
    assert decoded.aligned == 1 and decoded.mask == 0b1010
    assert (decoded.l1_ca, decoded.l1_p, decoded.l2_ca, decoded.l2_p) == (0, None, 0, None)
