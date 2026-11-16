"""The station: placement, schedule and observations behind one call per epoch.

DESIGN.md section 2.6::

    station = VirtualBase(engine, vrs_config)
    messages = station.epoch(t, rover_position)   # gnsswire.rtcm3 dataclasses, in send order
    wire = b"".join(rtcm3.encode(m) for m in messages)

The observations of an epoch (DESIGN 2.6.3) are produced by ``source``, an
:class:`ObservationSource`; by default ``vrs.observations.Observations`` on
the engine. The
station tells it, with each epoch, where the station is, its ID and whether
the epoch is the first of a placement, at which every arc and lock time
restarts (DESIGN 2.6.1); the station sets DF003 and DF393 of the MSMs it
returns and puts them in order (2.6.2).
"""

from __future__ import annotations

from typing import Any, Protocol, Sequence

from gnsswire import rtcm3

from ssrcore.model.types import GnssTime, RoverPosition, VrsConfig
from ssrcore.toolkit import Vector
from ssrcore.vrs.observations import Observations
from ssrcore.vrs.placement import Placement, Station
from ssrcore.vrs.schedule import Schedule, msm_system, order_msms, station_messages


class ObservationSource(Protocol):
    """What the station needs from the observations of DESIGN 2.6.3.

    ``observe(t, station)`` returns the MSMs of epoch ``t`` at
    ``station.position``, in any order; their DF003 and DF393 are the
    station's to set. When ``station.new`` is true, the source restarts every
    arc: lock times from zero, a new integer N_s per signal, wind-up unwrapping
    from scratch. [spec: RTCM 10403.3 Table 3.4-1 DF402, lock time restarts on
    any loss of continuity; design: DESIGN 2.6.1, a move resets every lock time]
    """

    def observe(self, t: GnssTime, station: Station) -> Sequence[rtcm3.Msm]:
        """[design: DESIGN 2.6.3] the MSMs of one epoch."""
        ...


class VirtualBase:
    """A virtual base station on the engine (DESIGN 2.6).

    ``engine`` is the engine the observations are drawn from; ``source``,
    when given, replaces the default ``Observations(engine, config)``. The
    station itself needs nothing of the engine. [design: DESIGN 2.6, the
    station uses only the public query, as any caller does; the source
    holds the engine queries, so placement and schedule are tested without
    one]
    """

    def __init__(self, engine: Any, config: VrsConfig, source: ObservationSource | None = None) -> None:
        """[design: DESIGN 2.6, configured by ``Config.vrs`` alone; the default source is
        ``vrs.observations.Observations``]"""
        self.engine = engine
        self.config = config
        self._source = Observations(engine, config) if source is None else source
        self._placement = Placement(config)
        self._schedule = Schedule(config)

    @property
    def station(self) -> Station | None:
        """[design: DESIGN 2.6.1] the current station, None before the first rover position."""
        return self._placement.station

    @property
    def source(self) -> ObservationSource:
        """[design: DESIGN 2.6] the observation source, for a caller that reads its state."""
        return self._source

    @property
    def schedule(self) -> Schedule:
        """[design: DESIGN 2.6.2] the epoch grid, for a caller choosing its epochs."""
        return self._schedule

    def epoch(self, t: GnssTime, rover: Vector | RoverPosition | None = None) -> list:
        """The messages of epoch ``t``, in send order; empty before the station is placed.

        ``rover`` is the rover's position for this epoch, ECEF metres or a
        GGA ``RoverPosition``, or None when no new position has come. ``t``
        must lie on the epoch grid (``schedule.next_epoch``) and follow the
        previous epoch.
        [spec: RTCM 10403.3 3.5.16.2, position messages before the
        observations of a station that has moved; Table 3.1-2, GNSS RTK
        service: MSM, 1005 or 1006, 1033, 1230 if GLONASS is sent; design:
        DESIGN 2.6.1, 2.6.2]
        """
        self._schedule.advance(t)
        station = self._placement.update(rover)
        if station is None:
            return []
        due = self._schedule.station_due(t, station.new)
        msms = order_msms(self._source.observe(t, station), station.station_id)
        out: list = []
        if due:
            out.extend(station_messages(station, {msm_system(m) for m in msms}))
        out.extend(msms)
        return out
