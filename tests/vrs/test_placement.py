"""Placement, movement and identity of the virtual base station (DESIGN 2.6.1, 2.6.4).

Placement at the first
rover position, the move after ``move_distance`` with its confirmation, and
the station ID change. Rover positions are synthetic, built in the local
frame of a chosen site with the toolkit, and the GGA of the u-blox F9 HPG
1.51 2.7.5.1 example.
"""

from __future__ import annotations

import math

import pytest
from gnsswire import nmea

from ssrcore.adapters.nmea import rover_position
from ssrcore.model.types import RoverPosition, VrsConfig
from ssrcore.toolkit.geodesy import ecef_to_geodetic, enu_to_ecef, geodetic_to_ecef
from ssrcore.vrs.placement import ARP_UNIT, Placement, Station, arp_units, horizontal_distance

# A site at 52° N, 13° E, 80 m on the ellipsoid.
SITE_GEODETIC = (math.radians(52.0), math.radians(13.0), 80.0)
SITE = geodetic_to_ecef(SITE_GEODETIC)


def offset(east: float, north: float, up: float = 0.0, origin=SITE_GEODETIC):
    """ECEF of a point east, north, up metres of the site (straight line in the local frame)."""
    base = geodetic_to_ecef(origin)
    d = enu_to_ecef((east, north, up), origin[0], origin[1])
    return (base[0] + d[0], base[1] + d[1], base[2] + d[2])


def test_rtcm_4_2_example_coordinates_in_arp_units():
    """RTCM 10403.3 4.2: the 1005 example's ARP (1114104.5999, −4850729.7108,
    3975521.4643 m) is these DF025 to DF027 integers (decoded from the example's
    bytes with gnsswire)."""
    assert arp_units((1114104.5999, -4850729.7108, 3975521.4643)) == (11141045999, -48507297108, 39755214643)


def test_arp_units_round_half_away_from_zero_and_refuse_out_of_range():
    assert arp_units((0.00005, -0.00005, 0.00004)) == (1, -1, 0)
    assert arp_units((12.34, -12.36, 0.05), step=1000) == (123000, -124000, 1000)
    with pytest.raises(ValueError):
        arp_units((13_743_895.4, 0.0, 0.0))
    with pytest.raises(ValueError):
        arp_units((math.nan, 0.0, 0.0))


def test_first_rover_position_places_the_station_rounded_to_a_tenth_of_a_metre():
    rover = offset(123.456, -78.987, 1.234)
    placement = Placement(VrsConfig(station_id=17))
    assert placement.update(None) is None
    station = placement.update(rover)
    assert station.new and station.placement == 1 and station.station_id == 17
    assert all(u % 1000 == 0 for u in station.arp)                     # 0.1 m in 0.1 mm units
    assert all(abs(s - r) <= 0.05 + 1e-9 for s, r in zip(station.position, rover))
    # the rover's ellipsoidal height, to the rounding
    assert abs(ecef_to_geodetic(station.position)[2] - ecef_to_geodetic(rover)[2]) < 0.09
    # new only at the first epoch of the placement
    again = placement.update(rover)
    assert again.arp == station.arp and not again.new


def test_position_is_exactly_what_the_1005_carries():
    station = Placement(VrsConfig()).update(offset(1.0, 2.0))
    assert station.position == tuple(u * ARP_UNIT for u in station.arp)


def test_placement_from_the_gga_example():
    """u-blox F9 HPG 1.51 2.7.5.1 example: 47°17.11399′ N, 8°33.91590′ E, 499.6 m
    above mean sea level with a geoid separation of 48.0 m, so 547.6 m on the ellipsoid."""
    gga = nmea.decode(nmea.parse(b"$GPGGA,092725.00,4717.11399,N,00833.91590,E,1,08,1.01,499.6,M,48.0,M,,*5B\r\n", 0))
    rover = rover_position(gga)
    station = Placement(VrsConfig()).update(rover)
    lat, lon, h = ecef_to_geodetic(station.position)
    assert math.degrees(lat) == pytest.approx(47 + 17.11399 / 60, abs=1e-6)
    assert math.degrees(lon) == pytest.approx(8 + 33.91590 / 60, abs=1e-6)
    assert h == pytest.approx(547.6, abs=0.09)


def test_horizontal_distance_ignores_height():
    assert horizontal_distance(SITE, offset(3000.0, 4000.0)) == pytest.approx(5000.0, abs=1e-6)
    assert horizontal_distance(SITE, offset(0.0, 0.0, 2500.0)) == pytest.approx(0.0, abs=1e-6)


def test_station_stays_while_the_rover_is_within_move_distance():
    placement = Placement(VrsConfig())
    first = placement.update(SITE)
    for k in range(100):
        s = placement.update(offset(4999.0 * math.cos(k), 4999.0 * math.sin(k), 300.0))
        assert s.arp == first.arp and s.station_id == first.station_id and not s.new
    assert placement.pending == 0


def test_a_height_change_alone_never_moves_the_station():
    placement = Placement(VrsConfig(move_confirm=1))
    first = placement.update(SITE)
    for _ in range(5):
        assert placement.update(offset(0.0, 0.0, 6000.0)).arp == first.arp


def test_move_after_move_confirm_consecutive_epochs_beyond_move_distance():
    """DESIGN 2.6.1: the station moves when the rover has been beyond 5 km for 10
    consecutive epochs; the confirming epoch is the first at the new place, which
    is that epoch's rover position, and the station ID increments."""
    config = VrsConfig(station_id=7)
    placement = Placement(config)
    first = placement.update(SITE)
    assert placement.confirm_epochs == 10                    # 10 s at 1 Hz
    rovers = [offset(5001.0 + 10.0 * k, 0.0) for k in range(placement.confirm_epochs)]
    for k, rover in enumerate(rovers[:-1], 1):
        s = placement.update(rover)
        assert s.arp == first.arp and s.station_id == 7 and not s.new
        assert placement.pending == k
    moved = placement.update(rovers[-1])
    assert moved.new and moved.placement == 2 and moved.station_id == 8
    assert moved.arp == arp_units(rovers[-1], step=1000)
    assert placement.pending == 0
    after = placement.update(rovers[-1])
    assert after.arp == moved.arp and not after.new and after.station_id == 8


def test_move_confirm_is_a_time_counted_in_epochs_of_the_rate():
    """DESIGN 2.6.1: ``move_confirm`` seconds are ⌈move_confirm · rate⌉ epochs, at least one."""
    assert Placement(VrsConfig(rate=2.0)).confirm_epochs == 20
    assert Placement(VrsConfig(rate=0.5)).confirm_epochs == 5
    assert Placement(VrsConfig(rate=0.5, move_confirm=3.0)).confirm_epochs == 2
    assert Placement(VrsConfig(move_confirm=0.0)).confirm_epochs == 1
    placement = Placement(VrsConfig(rate=2.0, move_confirm=1.5))
    placement.update(SITE)
    for _ in range(2):
        assert not placement.update(offset(6000.0, 0.0)).new
    assert placement.update(offset(6000.0, 0.0)).new


def test_a_position_back_within_move_distance_restarts_the_count():
    placement = Placement(VrsConfig(move_confirm=3))
    first = placement.update(SITE)
    far, near = offset(6000.0, 0.0), offset(100.0, 0.0)
    placement.update(far)
    placement.update(far)
    placement.update(near)
    assert placement.pending == 0
    placement.update(far)
    assert placement.update(far).arp == first.arp
    assert placement.update(far).new


def test_epochs_without_a_position_neither_confirm_nor_break_the_count():
    placement = Placement(VrsConfig(move_confirm=3))
    first = placement.update(SITE)
    far = offset(0.0, -7000.0)
    placement.update(far)
    for _ in range(20):
        assert placement.update(None).arp == first.arp
    assert placement.pending == 1
    placement.update(far)
    moved = placement.update(far)
    assert moved.new and moved.arp == arp_units(far, step=1000)


def test_station_id_wraps_modulo_4096():
    """RTCM 10403.3 DF003: 0 to 4095."""
    placement = Placement(VrsConfig(station_id=4095, move_confirm=1))
    assert placement.update(SITE).station_id == 4095
    assert placement.update(offset(10000.0, 0.0)).station_id == 0
    assert placement.update(offset(20000.0, 0.0)).station_id == 1


def test_each_move_is_measured_from_the_new_station():
    placement = Placement(VrsConfig(move_confirm=1))
    placement.update(SITE)
    second = placement.update(offset(6000.0, 0.0))
    assert placement.update(offset(10000.0, 0.0)).arp == second.arp     # 4 km from the new station
    third = placement.update(offset(11500.0, 0.0))
    assert third.placement == 3 and third.station_id == 2


def test_pinned_station_never_moves():
    pinned = (3_800_000.12345, 880_000.98765, 5_000_000.55555)
    placement = Placement(VrsConfig(pinned_position=pinned, station_id=99, move_confirm=1))
    assert placement.station is not None
    first = placement.update(None)
    assert first.new and first.station_id == 99
    assert first.arp == arp_units(pinned)                                # the 0.1 mm of DF025
    for _ in range(3):
        s = placement.update(offset(50_000.0, 0.0))
        assert s.arp == first.arp and s.station_id == 99 and not s.new


def test_rover_position_type_is_checked():
    with pytest.raises(ValueError):
        Placement(VrsConfig()).update((1.0, 2.0))
    assert isinstance(Placement(VrsConfig()).update(RoverPosition(0.9, 0.2, 100.0, 0.0)), Station)
