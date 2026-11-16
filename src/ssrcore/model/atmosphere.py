"""SPARTN area and atmosphere store (GAD, SM 2-0; HPAC, SM 1-0 to 1-4).

DESIGN.md sections 3.2 and 3.3, rules A1 to A8. The store keeps, per stream
(rule S1) and area ID (A3), the GAD definitions and the HPAC troposphere (per
area, A5) and ionosphere (per satellite, A6) in time order, each under the
issue it was sent with (a query at t uses the newest item not after t,
DESIGN 2.5): the SIOU as the stream state gives it (S3, S4) and the
AIOU with its generation (A2). A query chooses one area for a receiver
position (A4) and gives the zenith troposphere and the slant TEC of a
satellite in that area, polynomial and grid (A7), or the reason they are
absent. Mapping functions, the height reduction and the conversion of TEC to
metres per signal need the elevation, the day of year and the signal, and are
the engine's (toolkit ``troposphere`` and ``ionosphere``).

Times are GPS times (``GnssTime``): an input's time is its SPARTN time tag
resolved and converted by the engine (rules T1 to T3); the time of a query is
the query's ``t_rx``. The store reads no clock.

**Stream state.** The issue rules S3, S4 and A2 are the stream store's
(``model/stream.py``). The atmosphere store shares one ``Streams``
with the corrections store and notes each GAD and HPAC message there
(``Streams.observe``), as the corrections store notes each OCB message; the
``Issues`` it gets back key what it stores. A query is given the SIOU
``Issue`` the selection uses (DESIGN 3.4 step 2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ssrcore.model.types import (
    AreaDefinition,
    FeedOutcome,
    FeedReason,
    FeedResult,
    GadInput,
    GnssTime,
    HpacInput,
    IonosphereSatellite,
    MaxAge,
    Missing,
    Sat,
    StreamId,
    TroposphereBlock,
)
from ssrcore.model.stream import AIOU_LIFETIME, Issue, Issues, Streams
from ssrcore.toolkit import area as _area


@dataclass(frozen=True, slots=True)
class Absent:
    """No value, and why. [design: DESIGN 2.4, absent with a reason, never silently zero]"""

    reason: Missing


@dataclass(frozen=True, slots=True)
class AreaChoice:
    """The area a query uses (rule A4), with the issue whose data it combines (A1).

    ``latitude`` and ``longitude`` are the receiver's, in degrees.
    """

    stream: StreamId
    area_id: int
    definition: AreaDefinition
    issue: Issues
    latitude: float
    longitude: float


@dataclass(frozen=True, slots=True)
class TroposphereZenith:
    """Zenith troposphere of an area at the receiver, metres (SPARTN 8.12).

    ``hydrostatic`` is the area's average hydrostatic delay at 0 m on the
    ellipsoid (SF043 with its offset); ``polynomial`` the wet residual T_P;
    ``grid`` the interpolated wet residual T_G, 0.0 without a grid block and
    None when absent (``grid_missing`` says why; rule A7). ``quality`` is the
    SF042 index; ``age`` seconds from the HPAC time to the query time.
    """

    area_id: int
    hydrostatic: float
    polynomial: float
    grid: float | None
    grid_missing: Missing | None
    quality: int
    time: GnssTime
    age: float

    @property
    def wet(self) -> float:
        """[spec: SPARTN v2.0.3 8.12, T_P + T_G; A7, the polynomial alone without a usable grid]"""
        return self.polynomial + (self.grid or 0.0)


@dataclass(frozen=True, slots=True)
class IonosphereTec:
    """Slant TEC of a satellite in an area at the receiver, TECU (SPARTN 8.11).

    ``polynomial`` is I_P; ``grid`` I_G, 0.0 without a grid block and None
    when absent (``grid_missing`` says why; rule A7). ``quality`` is the
    SF055 index; ``age`` seconds from the HPAC time to the query time.
    """

    sat: Sat
    area_id: int
    polynomial: float
    grid: float | None
    grid_missing: Missing | None
    quality: int
    time: GnssTime
    age: float

    @property
    def slant(self) -> float:
        """[spec: SPARTN v2.0.3 8.11, I_P + I_G; A7, the polynomial alone without a usable grid]"""
        return self.polynomial + (self.grid or 0.0)


def area_of(d: AreaDefinition) -> _area.Area:
    """[spec: SPARTN v2.0.3 Table 6.24, SF032 to SF037] the toolkit's ``Area`` of a GAD area."""
    return _area.Area(d.latitude, d.longitude, d.latitude_nodes, d.longitude_nodes,
                      d.latitude_spacing, d.longitude_spacing)


@dataclass(slots=True)
class _Held:
    """A stored item and the time of the message it came in."""

    value: object
    time: GnssTime
    equation_type: int = 0


class AtmosphereStore:
    """GAD areas and HPAC atmosphere by stream and area ID (DESIGN 3.2, rules A1 to A8)."""

    def __init__(self, streams: Streams, max_age: float = MaxAge().atmosphere) -> None:
        """[design: DESIGN 2.5, ``config.max_age.atmosphere`` limits the age of HPAC data;
        DESIGN 3.2, one stream store serves the corrections and the atmosphere stores]"""
        if not 0.0 < max_age < AIOU_LIFETIME:
            raise ValueError(f"max_age {max_age} s outside (0, {AIOU_LIFETIME}) (SF068)")
        self.streams = streams
        self._max_age = float(max_age)
        # issue -> latest time a GAD or HPAC message of that issue was received
        self._seen: dict[tuple[StreamId, Issues], GnssTime] = {}
        self._gad: dict[tuple[StreamId, int], dict[Issues, _Held]] = {}
        # HPAC items per issue, in time order (``_add``)
        self._tropo: dict[tuple[StreamId, int], dict[Issues, list[_Held]]] = {}
        self._iono: dict[tuple[StreamId, int, Sat], dict[Issues, list[_Held]]] = {}

    # ---- input ---------------------------------------------------------------------------

    def _issue(self, stream: StreamId, siou: int, aiou: int, time: GnssTime) -> Issues | None:
        """The issues of a message, noted in the stream store; None when it is stale there.

        [spec: SPARTN v2.0.3 6.2, Table 6.2 SF005, SF068 (rules S3, S4, A2), through
        ``Streams.observe``]
        """
        issues = self.streams.observe(stream, time, siou, aiou)
        if issues is None:
            return None
        seen = self._seen.get((stream, issues))
        if seen is None or time > seen:
            self._seen[(stream, issues)] = time
        return issues

    def _last(self, stream: StreamId) -> GnssTime | None:
        """[spec: SPARTN Table 6.2 SF068] the time of the stream's latest GAD or HPAC message."""
        s = self.streams.state(stream)
        return None if s is None else s.aiou.last

    @staticmethod
    def _put(slot: dict[Issues, _Held], issue: Issues, held: _Held) -> bool:
        """[spec: SPARTN v2.0.3 6.5.3 (rule A3): a GAD area is replaced by the latest
        definition of its issue; design: a message of an older epoch does not replace a
        newer one]"""
        old = slot.get(issue)
        if old is not None and held.time < old.time:
            return False
        slot[issue] = held
        return True

    @staticmethod
    def _add(slot: dict[Issues, list[_Held]], issue: Issues, held: _Held) -> bool:
        """Keep an HPAC item in the time order of its issue.

        [spec: SPARTN v2.0.3 8.4 item 2 (rule A5): the latest received model of an epoch
        is used, so an item of the same time tag replaces the stored one; Table 7.1
        TF009, the tag is the time of the data; design: DESIGN 2.5, an item applies from
        its time tag until the next, so the earlier items of the issue are kept for a
        query before the newest tag, as the corrections store keeps its blocks; an item
        older than the newest is refused, as there (C3)]
        """
        items = slot.setdefault(issue, [])
        if items and held.time < items[-1].time:
            return False
        if items and held.time == items[-1].time:
            items[-1] = held
        else:
            items.append(held)
        return True

    @staticmethod
    def _at(items: list[_Held], t: GnssTime) -> _Held | None:
        """[spec: SPARTN Table 7.1 TF009; design: DESIGN 2.5, an item applies from its
        time tag until the next] the newest item not after ``t``."""
        for h in reversed(items):
            if h.time <= t:
                return h
        return None

    def _prune(self, stream: StreamId) -> None:
        """Forget what can no longer be used: HPAC data older than the stream's latest
        atmosphere message by the AIOU outage, and GAD definitions of an issue not seen
        for as long.

        [spec: SPARTN v2.0.3 Table 6.2 SF068 (rule A2); design: such items fail the age
        limit of DESIGN 2.5 at every later query, so keeping them only costs memory]
        """
        last = self._last(stream)
        assert last is not None   # the stream has just been observed with an AIOU
        for table in (self._tropo, self._iono):
            for key in [k for k in table if k[0] == stream]:
                slot = table[key]
                for issue in list(slot):
                    slot[issue] = [h for h in slot[issue] if last - h.time < AIOU_LIFETIME]
                    if not slot[issue]:
                        del slot[issue]
                if not slot:
                    del table[key]
        old = {i for (s, i), t in self._seen.items() if s == stream and last - t >= AIOU_LIFETIME}
        for key in [k for k in self._gad if k[0] == stream]:
            slot = self._gad[key]
            for issue in [i for i in slot if i in old]:
                del slot[issue]
            if not slot:
                del self._gad[key]
        for issue in old:
            del self._seen[(stream, issue)]

    @staticmethod
    def _result(stored: list[int], sats: list[Sat], older: int) -> FeedResult:
        """[design: DESIGN 2.1, the result says what the message changed]"""
        if not stored and older:
            return FeedResult(FeedOutcome.SUPERSEDED, FeedReason.OLDER,
                              detail="every block is older than the one stored")
        return FeedResult(FeedOutcome.STORED, sats=tuple(sats), areas=tuple(dict.fromkeys(stored)),
                          detail=f"{older} older block(s) not stored" if older else "")

    def feed_gad(self, gad: GadInput, time: GnssTime) -> FeedResult:
        """Store the area definitions of a GAD message whose time tag is ``time``.

        [spec: SPARTN v2.0.3 6.5.3, Table 6.24; Table 6.2 SF031 (rule A3: areas
        keyed by stream and area ID), SF005, SF068 (rules A1, A2)]
        """
        issue = self._issue(gad.stream, gad.siou, gad.aiou, time)
        if issue is None:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.STALE,
                              detail="older than the stream's latest by the issue lifetime (SF005, SF068)")
        stored, older = [], 0
        for d in gad.areas:
            slot = self._gad.setdefault((gad.stream, d.area_id), {})
            if self._put(slot, issue, _Held(d, time)):
                stored.append(d.area_id)
            else:
                older += 1
        self._prune(gad.stream)
        return self._result(stored, [], older)

    def feed_hpac(self, hpac: HpacInput, time: GnssTime) -> FeedResult:
        """Store the troposphere and ionosphere of an HPAC message whose time tag is ``time``.

        The troposphere is kept per area, whichever constellation's message
        carried it (rule A5); the ionosphere per satellite of the message's
        constellation (A6). An area without a troposphere or ionosphere block
        leaves the stored one as it is. [spec: SPARTN v2.0.3 6.5.2, Tables
        6.13 to 6.22; 8.4 item 2 (A5); Table 6.2 SF031, SF005, SF068 (A1 to
        A3); design: a missing block is not a withdrawal, since 8.4 item 2 lets
        providers send the troposphere in one constellation's message only]
        """
        issue = self._issue(hpac.stream, hpac.siou, hpac.aiou, time)
        if issue is None:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.STALE,
                              detail="older than the stream's latest by the issue lifetime (SF005, SF068)")
        stored: list[int] = []
        sats: list[Sat] = []
        older = 0
        for a in hpac.areas:
            if a.troposphere is not None:
                slot = self._tropo.setdefault((hpac.stream, a.area_id), {})
                if self._add(slot, issue, _Held(a.troposphere, time, a.troposphere.equation_type)):
                    stored.append(a.area_id)
                else:
                    older += 1
            if a.ionosphere is not None:
                for s in a.ionosphere.satellites:
                    slot = self._iono.setdefault((hpac.stream, a.area_id, s.sat), {})
                    if self._add(slot, issue, _Held(s, time, a.ionosphere.equation_type)):
                        stored.append(a.area_id)
                        sats.append(s.sat)
                    else:
                        older += 1
        self._prune(hpac.stream)
        return self._result(stored, list(dict.fromkeys(sats)), older)

    # ---- queries -------------------------------------------------------------------------

    def _usable(self, held: _Held, t: GnssTime) -> bool:
        """[spec: SPARTN v2.0.3 Table 6.2 SF068 (450 s); design: DESIGN 2.5, an item applies
        from its time tag and for at most ``max_age``, which is below 450 s]"""
        return 0.0 <= t - held.time <= self._max_age

    def _data(self, stream: StreamId, area_id: int, t: GnssTime) -> list[tuple[Issues, _Held]]:
        """The HPAC items of an area that apply at ``t`` (``_at``), troposphere and
        ionosphere, each with its issue; an issue whose items all come after ``t`` has none.
        [design: DESIGN 2.5, an item applies from its time tag]"""
        slots = [self._tropo.get((stream, area_id), {})]
        slots += [slot for key, slot in self._iono.items() if key[0] == stream and key[1] == area_id]
        out = []
        for slot in slots:
            for issue, items in slot.items():
                held = self._at(items, t)
                if held is not None:
                    out.append((issue, held))
        return out

    def _held(self, table: dict, key: tuple, issue: Issues, t: GnssTime) -> _Held | None:
        """[design: ``_at`` of one stored slot and issue]"""
        items = table.get(key, {}).get(issue)
        return None if items is None else self._at(items, t)

    def area(self, stream: StreamId, siou: Issue, latitude: float, longitude: float,
             t: GnssTime) -> AreaChoice | Absent:
        """The area whose atmosphere serves a receiver at (``latitude``, ``longitude``),
        degrees, at time ``t``, for the SIOU issue ``siou`` the selection uses.

        The receiver must be inside the area and within ±85° latitude; the GAD and
        the HPAC data must have that SIOU and the same AIOU, with no outage of 450 s
        or more up to ``t``. Of the areas that qualify and hold data usable at
        ``t``, the one whose centroid is nearest is chosen (ties: the lowest area
        ID); within an area, the issue with the newest data. The same choice serves
        the troposphere and the ionosphere of the query.

        [spec: SPARTN v2.0.3 8.10 (inside the area, ±85°), 6.2 and Table 6.2
        SF005, SF068 (rules A1, A2), SF031 (A3); design (A4): the nearest centroid
        gives the most interior polynomial expansion, and one area per query keeps
        the terms consistent; distance on the local plane, δφ² + (δλ cos φ)², which
        is enough for areas at most 22.4° wide; an area without usable data cannot
        serve, so it is passed over]
        """
        last = self._last(stream)
        if last is None:
            return Absent(Missing.NO_ATMOSPHERE)
        if t - last >= AIOU_LIFETIME:
            return Absent(Missing.STALE)
        gads = [(key[1], issue, held.value) for key, slot in self._gad.items() if key[0] == stream
                for issue, held in slot.items()]
        if not gads:
            return Absent(Missing.NO_ATMOSPHERE)
        gads = [g for g in gads if g[1].siou == siou]
        if not gads:
            return Absent(Missing.SIOU_MISMATCH)
        gads = [g for g in gads if _area.contains(area_of(g[2]), latitude, longitude)]
        if not gads:
            return Absent(Missing.LATITUDE_LIMIT if abs(latitude) > 85.0
                          else Missing.OUTSIDE_AREA)
        candidates = []
        reasons = set()
        for area_id, issue, definition in gads:
            data = self._data(stream, area_id, t)
            times = [h.time for i, h in data if i == issue and self._usable(h, t)]
            if times:
                candidates.append((area_id, issue, definition, max(times)))
            elif any(i == issue for i, _ in data):
                reasons.add(Missing.STALE)
            elif any(i.siou == siou for i, _ in data):
                reasons.add(Missing.AIOU_MISMATCH)
            elif data:
                reasons.add(Missing.SIOU_MISMATCH)
            else:
                reasons.add(Missing.NO_ATMOSPHERE)
        if not candidates:
            for reason in (Missing.STALE, Missing.AIOU_MISMATCH,
                           Missing.SIOU_MISMATCH, Missing.NO_ATMOSPHERE):
                if reason in reasons:
                    return Absent(reason)
        newest: dict[int, tuple] = {}
        for c in candidates:
            if c[0] not in newest or c[3] > newest[c[0]][3]:
                newest[c[0]] = c
        cos_lat = math.cos(math.radians(latitude))

        def distance(c: tuple) -> tuple[float, int]:
            dlat, dlon = _area.offsets(area_of(c[2]), latitude, longitude)
            return dlat * dlat + (dlon * cos_lat) ** 2, c[0]

        area_id, issue, definition, _ = min(newest.values(), key=distance)
        return AreaChoice(stream, area_id, definition, issue, latitude, longitude)

    @staticmethod
    def _grid(area: _area.Area, residuals: tuple[float | None, ...] | None, latitude: float,
              longitude: float) -> tuple[tuple[float, float, float, float] | None,
                                         tuple[float, float, float, float] | None,
                                         Missing | None]:
        """(weights, residuals of the four nodes, reason) for the grid, or zero weights
        without a grid block.

        [spec: SPARTN v2.0.3 8.11, 8.12 ("in the case where no grid information is
        present, I_G is assumed to be zero"); node order of the 8.10 figure and the
        weights of 8.13 through ``toolkit.area.grid``; design (rule A7): a residual
        count other than the GAD's node count cannot be placed on the nodes, and a
        node at the invalid value leaves the grid residual absent, not extrapolated]
        """
        if residuals is None:
            return None, None, None
        if len(residuals) != area.latitude_nodes * area.longitude_nodes:
            return None, None, Missing.GRID_MISMATCH
        nodes, weights = _area.grid(area, latitude, longitude)
        values = tuple(residuals[n] for n in nodes)
        if any(v is None for v in values):
            return None, None, Missing.INVALID_GRID_NODE
        return weights, values, None  # type: ignore[return-value]

    def troposphere(self, choice: AreaChoice, t: GnssTime) -> TroposphereZenith | Absent:
        """The zenith troposphere of the chosen area at the receiver, at time ``t``.

        [spec: SPARTN v2.0.3 8.12 (T_h, T_P from the polynomial of SF041 about the
        area centroid, T_G bilinear from the grid); rules A1, A5, A7]
        """
        held = self._held(self._tropo, (choice.stream, choice.area_id), choice.issue, t)
        if held is None:
            return Absent(Missing.NO_ATMOSPHERE)
        if not self._usable(held, t):
            return Absent(Missing.STALE)
        block: TroposphereBlock = held.value  # type: ignore[assignment]
        area = area_of(choice.definition)
        dlat, dlon = _area.offsets(area, choice.latitude, choice.longitude)
        polynomial = _area.polynomial(block.equation_type, block.coefficients, dlat, dlon)
        weights, values, why = self._grid(area, block.residuals, choice.latitude, choice.longitude)
        if why is not None:
            grid = None
        elif weights is None:
            grid = 0.0
        else:
            grid = _area.interpolate(weights, values)
        return TroposphereZenith(choice.area_id, block.hydrostatic, polynomial, grid, why,
                                 block.quality, held.time, t - held.time)

    def ionosphere(self, choice: AreaChoice, sat: Sat, t: GnssTime) -> IonosphereTec | Absent:
        """The slant TEC of ``sat`` in the chosen area at the receiver, at time ``t``.

        [spec: SPARTN v2.0.3 8.11 (I_P from the polynomial of SF054 about the area
        centroid, I_G bilinear from the grid); 6.5.2 Table 6.19 (rule A6, the
        satellite's block in its constellation's HPAC); rules A1, A7; design: each
        satellite block is complete in itself, so a satellite keeps its latest
        block until the age limit even when a later message omits it]
        """
        held = self._held(self._iono, (choice.stream, choice.area_id, sat), choice.issue, t)
        if held is None:
            return Absent(Missing.NO_ATMOSPHERE)
        if not self._usable(held, t):
            return Absent(Missing.STALE)
        block: IonosphereSatellite = held.value  # type: ignore[assignment]
        area = area_of(choice.definition)
        dlat, dlon = _area.offsets(area, choice.latitude, choice.longitude)
        polynomial = _area.polynomial(held.equation_type, block.coefficients, dlat, dlon)
        weights, values, why = self._grid(area, block.residuals, choice.latitude, choice.longitude)
        if why is not None:
            grid = None
        elif weights is None:
            grid = 0.0
        else:
            grid = _area.interpolate(weights, values)
        return IonosphereTec(sat, choice.area_id, polynomial, grid, why, block.quality,
                             held.time, t - held.time)
