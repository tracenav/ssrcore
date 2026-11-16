"""Where the station is, when it moves, and what identifies it.

DESIGN.md sections 2.6.1 and 2.6.4. The station is a point on a grid of
0.1 mm, the resolution of the ARP coordinates of RTCM 1005 (DF025 to DF027),
so the position the observations are computed at is exactly the position the
receiver decodes. A placement made from the rover is rounded to 0.1 m
(DESIGN 2.6.1), which lies on that grid.

The distance is horizontal, a move needs the rover beyond it for
``move_confirm`` seconds of consecutive epochs, and the first placement keeps
the configured station ID.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ssrcore.model.types import RoverPosition, VrsConfig
from ssrcore.toolkit import Vector
from ssrcore.toolkit.geodesy import ecef_to_enu, ecef_to_geodetic, geodetic_to_ecef

# [spec: RTCM 10403.3 Table 3.4-1, DF025 to DF027: resolution 0.0001 m] the ARP
# coordinates in units of 0.1 mm, and their range, ±13 743 895.3471 m (int38).
ARP_UNIT = 0.0001
_ARP_LIMIT = (1 << 37) - 1
# [design: DESIGN 2.6.1, 0.1 m rounding keeps the coordinates stable across GGA
# noise] one placement step, in ARP units.
_PLACEMENT_STEP = 1000
# [spec: RTCM 10403.3 Table 3.4-1, DF003 range 0 to 4095] the number of station IDs.
STATION_IDS = 4096


@dataclass(frozen=True, slots=True)
class Station:
    """The station of one placement: its ARP in DF025 to DF027 units (0.1 mm),
    its reference station ID (DF003), the count of placements so far (1 for
    the first) and whether this epoch is the first of the placement.

    [spec: RTCM 10403.3 3.5.3, the ARP is the point the coordinates refer to;
    Table 3.4-1 DF003, DF025 to DF027; design: DESIGN 2.6.1, the placement
    count and ``new`` tell the observations when every arc restarts]
    """

    arp: tuple[int, int, int]
    station_id: int
    placement: int
    new: bool

    @property
    def position(self) -> Vector:
        """[spec: RTCM 10403.3 Table 3.4-1, DF025 to DF027 in 0.0001 m] the ARP in ECEF metres,
        exactly the value a decoder of the 1005 obtains."""
        x, y, z = self.arp
        return (x * ARP_UNIT, y * ARP_UNIT, z * ARP_UNIT)


def arp_units(position: Vector, step: int = 1) -> tuple[int, int, int]:
    """ECEF metres to DF025 to DF027 units, each coordinate rounded to the nearest
    multiple of ``step`` units (half away from zero).

    [spec: RTCM 10403.3 Table 3.4-1, DF025 to DF027: 0.0001 m, int38; design:
    DESIGN 2.6.1, ``step`` 1000 gives the 0.1 m of a placement; rounding half
    away from zero treats both signs alike]
    """
    if len(position) != 3:
        raise ValueError("an ARP is an ECEF (x, y, z)")
    units = []
    for value in position:
        if not math.isfinite(value):
            raise ValueError(f"ARP coordinate {value} is not finite")
        scaled = abs(value) / (ARP_UNIT * step)
        n = int(math.floor(scaled + 0.5)) * step
        n = n if value >= 0 else -n
        if abs(n) > _ARP_LIMIT:
            raise ValueError(f"ARP coordinate {value} m outside the int38 range of DF025 to DF027")
        units.append(n)
    return (units[0], units[1], units[2])


def rover_ecef(rover: Vector | RoverPosition) -> Vector:
    """The rover as ECEF metres: an ECEF tuple as given, a GGA position
    (``adapters.nmea.rover_position``) through WGS-84.

    [spec: RTCM 10403.3 3.5.12.1 and 3.5.12.2, the rover sends its approximate
    position as NMEA GGA; TM-23 B.1.1 for the conversion]
    """
    if isinstance(rover, RoverPosition):
        return geodetic_to_ecef((rover.latitude, rover.longitude, rover.height))
    if len(rover) != 3:
        raise ValueError("a rover position is an ECEF (x, y, z) or a RoverPosition")
    return (float(rover[0]), float(rover[1]), float(rover[2]))


def horizontal_distance(station: Vector, rover: Vector) -> float:
    """Distance from the station to the rover in the station's local horizontal plane.

    [spec: TM-23 B.1.2, B.2.2 (local frame at the station); design: DESIGN 2.6.1
    measures the move horizontally, so a change of height alone, which the
    receiver's own troposphere model handles, never moves the station; at 5 km
    the plane departs from the ellipsoid by under 2 m vertically, and the
    horizontal distance by under 1 mm]
    """
    latitude, longitude, _ = ecef_to_geodetic(station)
    east, north, _ = ecef_to_enu((rover[0] - station[0], rover[1] - station[1], rover[2] - station[2]),
                                 latitude, longitude)
    return math.hypot(east, north)


class Placement:
    """The station's place and identity, updated once per epoch from the rover.

    [design: DESIGN 2.6.1 and 2.6.4]

    - The first rover position places the station there, rounded to 0.1 m,
      with the configured ``station_id``.
    - The station stays while the rover is within ``move_distance``
      horizontally; it moves to the rover once the rover has been beyond it
      for ``move_confirm`` seconds: ``confirm_epochs``, ⌈move_confirm · rate⌉
      and at least one, consecutive epochs that carried a position. The epoch
      that confirms is the first at the new place, and the station ID
      increments by one modulo 4096. [design: DESIGN 2.6.1, a time gives the
      same protection at every rate]
    - An epoch without a rover position neither confirms nor breaks the
      count. [design: an absent GGA says nothing about where the rover is,
      and resetting on it would delay a move indefinitely on a link that
      drops sentences]
    - With ``pinned_position`` the station is there from the start, at the
      0.1 mm of the ARP, and never moves.
    """

    def __init__(self, config: VrsConfig) -> None:
        """[design: DESIGN 2.6.1, a pinned station exists before any rover position]"""
        self._config = config
        self._station: Station | None = None
        self._beyond = 0
        # a hair below the product, so that 10 s at 1 Hz is 10 epochs despite rounding
        self.confirm_epochs = max(1, math.ceil(config.move_confirm * config.rate - 1e-9))
        if config.pinned_position is not None:
            self._station = Station(arp_units(config.pinned_position), config.station_id, 1, True)

    @property
    def station(self) -> Station | None:
        """[design: DESIGN 2.6.1] the current station; None until the first rover position."""
        return self._station

    @property
    def pending(self) -> int:
        """[design: DESIGN 2.6.1] consecutive epochs the rover has been beyond ``move_distance``."""
        return self._beyond

    def update(self, rover: Vector | RoverPosition | None) -> Station | None:
        """The station for this epoch, given the rover position of the epoch (or None).

        [spec: RTCM 10403.3 3.5.12.1, the station is generated at the rover's
        approximate position; Table 3.4-1 DF003, 0 to 4095; design: DESIGN
        2.6.1 and 2.6.4, see the class]
        """
        config = self._config
        station = self._station
        if config.pinned_position is not None:
            return self._emit(station)
        if rover is None:
            return self._emit(station)
        here = rover_ecef(rover)
        if station is None:
            self._beyond = 0
            return self._place(here, config.station_id, 1)
        if horizontal_distance(station.position, here) <= config.move_distance:
            self._beyond = 0
            return self._emit(station)
        self._beyond += 1
        if self._beyond < self.confirm_epochs:
            return self._emit(station)
        self._beyond = 0
        return self._place(here, (station.station_id + 1) % STATION_IDS, station.placement + 1)

    def _place(self, here: Vector, station_id: int, count: int) -> Station:
        return self._emit(Station(arp_units(here, _PLACEMENT_STEP), station_id, count, True))

    def _emit(self, station: Station | None) -> Station | None:
        """The station for this epoch, ``new`` at the first epoch of a placement only;
        what is kept for the next epoch is no longer new."""
        if station is not None:
            self._station = Station(station.arp, station.station_id, station.placement, False)
        return station
