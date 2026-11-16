"""Pseudorange, phase range, rate and lock time of each signal at the station.

DESIGN.md section 2.6.3. ``Observations`` is the observation source of
``VirtualBase`` (``vrs/station.py``): at each epoch it queries the engine at
the station's ARP for every satellite with a broadcast ephemeris, forms

    P_s = pseudorange_model(s)
    Φ_s = phase_model(s) + λs·(k + N_s)

from the line-of-sight record (DESIGN 2.4; receiver clock and receiver terms
zero), keeps the arc of every signal (its integer N_s, its start for the lock
time, what ends it) and the unwrapped wind-up of every satellite, and encodes
the epoch as MSMs (``vrs/msm.py``).

The lock time follows the engine's discontinuity counts and the arc rules
of DESIGN 2.6.3 rather than a cycle-slip detector on the synthesised phases
(which cannot slip), the integer is chosen at the start of an arc, and no
CNR is invented.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from gnsswire import rtcm3

from ssrcore.model.providers import Kind
from ssrcore.model.types import SPEED_OF_LIGHT, GnssTime, Sat, Signal, System, Term, VrsConfig
from ssrcore.record import LineOfSight, MissingTermError
from ssrcore.toolkit import windup as _windup
from ssrcore.vrs import msm
from ssrcore.vrs.placement import Station

# [spec: RTCM 10403.3 Table 3.5-78 notes, "X≤64": the cells of one MSM]
MAX_CELLS = 64


@dataclass(frozen=True, slots=True)
class Observation:
    """One signal of one satellite at the station, before encoding.

    ``pseudorange`` P_s and ``phase`` Φ_s in metres (``phase`` None when
    withheld, DESIGN 2.6.3), ``cycles`` the whole cycles k + N_s added to the
    record's phase model, ``wavelength`` λs (m), ``lock`` the milliseconds
    since the arc began (0 without a phase), ``rate`` the phase range rate
    (m/s), ``channel`` the GLONASS frequency channel.
    [design: DESIGN 2.6.3, the observation as the record and the arc give it]
    """

    sat: Sat
    signal: Signal
    pseudorange: float
    phase: float | None
    cycles: int | None
    wavelength: float
    lock: int
    rate: float
    channel: int | None


@dataclass(frozen=True, slots=True)
class _Arc:
    """A signal's phase arc: when it began, its integer N_s and what ends it (DESIGN 2.6.3),
    ``origins`` the origins of the record's troposphere and ionosphere (item 7)."""

    start: GnssTime
    ambiguity: int
    discontinuities: int
    fixed: bool
    siou: int | None
    origins: tuple = (None, None)


# [design: DESIGN 2.6.3 *Satellites*, the reasons a satellite whose terms are present
# is left out because a term is not one the station may use (rule A11)]
TROPOSPHERE_NOT_A_CORRECTION = "troposphere_not_a_correction"
IONOSPHERE_NOT_WITH_BIASES = "ionosphere_not_with_biases"


@dataclass(frozen=True, slots=True)
class _Satellite:
    """One satellite's cells in an MSM: rough range, rate, channel and per signal ID the cell."""

    sat_id: int
    rough: tuple[int, int]
    rate: tuple[int, int]
    channel: int | None
    cells: dict          # signal ID -> (fine pseudorange, fine phase range, lock indicator)


class Observations:
    """The observation source of the virtual base station (DESIGN 2.6.3).

    [design: DESIGN 2.6, it uses only the engine's public query; the state
    kept from one epoch to the next is the arcs and the wind-up, nothing
    else, so the same inputs and epochs give the same messages]
    """

    def __init__(self, engine: Any, config: VrsConfig) -> None:
        """[design: DESIGN 2.6, the engine to query and the station's configuration]"""
        self.engine = engine
        self.config = config
        self._arcs: dict[tuple[Sat, Signal], _Arc] = {}
        self._windup: dict[Sat, float] = {}
        # [design: the observations of the last epoch, for a caller comparing them
        # with the engine (the tests)]
        self.last: list[Observation] = []
        # [design: DESIGN 2.6.3 *Satellites*, why each satellite with corrections above
        # the mask was left out at the last epoch for its atmosphere or tide]
        self.left: dict[Sat, str] = {}

    # ---- the observations ---------------------------------------------------------

    def observations(self, t: GnssTime, station: Station) -> list[Observation]:
        """The observations of epoch ``t`` at ``station.position``, by satellite and signal.

        Advances the arcs: call once per epoch, in time order. [spec: RTCM 10403.3
        DF401 (the integer kept for the arc), DF402 (lock time); design: DESIGN
        2.6.3, every rule of the satellites, signals, phase and arcs]
        """
        if station.new:
            # [design: DESIGN 2.6.1, a move restarts every arc]
            self._arcs.clear()
            self._windup.clear()
        out: list[Observation] = []
        arcs: dict[tuple[Sat, Signal], _Arc] = {}
        windup: dict[Sat, float] = {}
        left: dict[Sat, str] = {}
        for sat in self.engine.ephemerides.satellites():
            if msm.satellite_id(sat) is None:
                continue
            los = self.engine.query(sat, t, station.position)
            out.extend(self._satellite(t, los, arcs, windup, left))
        self._arcs, self._windup = arcs, windup
        self.last, self.left = out, left
        return out

    @staticmethod
    def _unusable(los: LineOfSight) -> str | None:
        """Why the station may not use the record's troposphere, ionosphere or tide, or None.

        [design: DESIGN 2.6.3 *Satellites*: a term missing leaves the satellite out with
        the record's reason; a troposphere that is not a correction, or an ionosphere
        that is not a correction consistent with the code biases (rule A11), leaves it
        out as well, so that every observation the station sends is formed from
        corrections of one consistent set: a model is not a correction of this
        atmosphere (UNB3m's zenith error is 4.9 cm, 1σ), and an ionosphere from
        another source than the code biases differs from one consistent with them
        by a constant per satellite]
        """
        tide = None if los.modelled is None else los.modelled.solid_tide
        for term, value in ((Term.TROPOSPHERE, los.troposphere), (Term.IONOSPHERE, los.ionosphere),
                            (Term.SOLID_TIDE, tide)):
            if value is None:
                return los.why(term) or term     # a record built without its reason: the term
        tropo, iono = los.troposphere, los.ionosphere
        if tropo.origin is None or tropo.origin.kind != Kind.CORRECTION:
            return TROPOSPHERE_NOT_A_CORRECTION
        if iono.origin is None or iono.origin.kind != Kind.CORRECTION or iono.with_biases is False:
            return IONOSPHERE_NOT_WITH_BIASES
        return None

    def _satellite(self, t: GnssTime, los: LineOfSight, arcs: dict, windup: dict,
                   left: dict) -> list[Observation]:
        """[design: DESIGN 2.6.3 *Satellites*, *Signals*, *Phase*, *Integer*, *Arcs*]"""
        geometry = los.geometry
        if not los.corrected or geometry is None or geometry.elevation < self.config.elevation_mask:
            return []
        sat = los.sat
        reason = self._unusable(los)
        if reason is not None:
            left[sat] = reason
            return []
        origins = (los.troposphere.origin, los.ionosphere.origin)
        channel = los.satellite.frequency_channel
        # [spec: RTCM 10403.3 DF399, DF404, the sign of the derivative of the phase range;
        # design: DESIGN 2.6.3, ρ̇ − c·ḋt_sat, the rate of the satellite's modelled range]
        rate = geometry.range_rate - SPEED_OF_LIGHT * los.satellite.clock_drift
        siou = None if los.quality.stream is None else los.quality.stream.siou
        w = los.modelled.windup
        # [spec: TM-23 5.5 through toolkit.windup.unwrap; design: DESIGN 2.6.3, per
        # satellite, from scratch when the satellite had no phase at the last epoch]
        unwrapped = None if w is None else (_windup.unwrap(self._windup[sat], w) if sat in self._windup else w)
        turns = None if w is None else round(unwrapped - w)
        signals = sorted((s for s in los.signals if msm.signal_id(s) is not None), key=msm.signal_id)
        out = []
        rough = None
        for signal in signals:
            try:
                pseudorange = los.pseudorange_model(signal)
            except MissingTermError:
                continue       # [design: DESIGN 2.6.3, no pseudorange: the signal is not sent]
            if rough is None:
                rough = msm.rough_range(pseudorange)
                if rough is None:
                    return []
            rough_m = msm.rough_metres(*rough)
            if msm.fine_pseudorange(pseudorange, rough_m, self.config.msm) is None:
                continue       # [spec: RTCM 10403.3 3.5.16.3.4, a pseudorange that does not fit]
            wavelength = SPEED_OF_LIGHT / los.frequency[signal]
            phase = cycles = None
            lock = 0
            bias = los.phase_bias.get(signal)
            if bias is not None and turns is not None and (bias.fixed or self.config.float_phase):
                model = los.phase_model(signal) + wavelength * turns
                arc = self._continued(sat, signal, bias, siou, origins)
                if arc is None or msm.fine_phaserange(model + wavelength * arc.ambiguity, rough_m,
                                                      self.config.msm) is None:
                    # [spec: RTCM 10403.3 DF401, the integer that matches the phase range to
                    # the pseudorange at the start of the arc]
                    arc = _Arc(t, msm.round_half_away((pseudorange - model) / wavelength),
                               bias.discontinuities, bias.fixed, siou, origins)
                cycles = turns + arc.ambiguity
                phase = model + wavelength * arc.ambiguity
                lock = msm.round_half_away((t - arc.start) * 1000.0)
                arcs[(sat, signal)] = arc
                windup[sat] = unwrapped
            out.append(Observation(sat, signal, pseudorange, phase, cycles, wavelength, lock, rate, channel))
        return out

    def _continued(self, sat: Sat, signal: Signal, bias, siou: int | None, origins: tuple) -> _Arc | None:
        """The arc of the previous epoch when it goes on, None when it ends.

        [spec: RTCM 10403.3 DF402, 3.5.16.3.3; SPARTN Table 6.2 SF015 (rule C7);
        design: DESIGN 2.6.3 *Arcs* items 2 to 5 and 7 (an origin change of the
        troposphere or ionosphere, rule A9); item 1 is the clearing on a move,
        item 6 the field check of the caller]
        """
        arc = self._arcs.get((sat, signal))
        if arc is None or arc.discontinuities != bias.discontinuities or arc.fixed != bias.fixed \
                or arc.siou != siou or arc.origins != origins:
            return None
        return arc

    # ---- the messages -------------------------------------------------------------

    def observe(self, t: GnssTime, station: Station) -> list[rtcm3.Msm]:
        """The MSMs of epoch ``t`` at the station: the ``ObservationSource`` of ``vrs/station.py``.

        DF003 and DF393 are the station's to set (0 and 1 here).
        [spec: RTCM 10403.3 3.5.16.3, Tables 3.5-78 to 3.5-88; design: DESIGN 2.6.3]
        """
        return self.encode(t, self.observations(t, station))

    def encode(self, t: GnssTime, observations: list[Observation]) -> list[rtcm3.Msm]:
        """The MSMs of the configured level holding ``observations`` of epoch ``t``: per
        system in ascending satellite ID, split over several MSMs where the cells
        exceed 64; systems in the order of DESIGN 2.6.2.

        [spec: RTCM 10403.3 3.5.16.3.4, splitting satellite by satellite; Table
        3.5-78, X ≤ 64; design: DESIGN 2.6.3 *Fields of a cell*]
        """
        level = int(self.config.msm)
        by_system: dict[System, list[_Satellite]] = {}
        for sat in sorted({o.sat for o in observations}):
            mine = sorted((o for o in observations if o.sat == sat), key=lambda o: msm.signal_id(o.signal))
            by_system.setdefault(sat.system, []).append(self._cells(mine, level))
        leap = self.engine.timebase.leap_at(t) if System.GLONASS in by_system else None
        out = []
        for system in (System.GPS, System.GLONASS, System.GALILEO, System.BEIDOU, System.QZSS):
            satellites = by_system.get(system)
            if not satellites:
                continue
            fields = msm.epoch_fields(system, t, leap)
            if fields is None:
                continue       # [spec: rule T4; GLONASS epoch time needs the leap seconds]
            for group in _groups(satellites):
                out.append(_message(system, level, fields, group))
        return out

    @staticmethod
    def _cells(observations: list[Observation], level: int) -> _Satellite:
        """[spec: RTCM 10403.3 DF397 to DF408, DF419, DF420; design: DESIGN 2.6.3 *Fields of a cell*]"""
        first = observations[0]
        rough = msm.rough_range(first.pseudorange)
        rough_m = msm.rough_metres(*rough)
        cells = {}
        for o in observations:
            fine_p = msm.fine_pseudorange(o.pseudorange, rough_m, level)
            if o.phase is None:
                fine_l, indicator = msm.invalid_phaserange(level), 0
            else:
                fine_l, indicator = msm.fine_phaserange(o.phase, rough_m, level), msm.lock(o.lock, level)
                if fine_l is None:     # cannot happen for an arc checked in ``_satellite``; kept invalid
                    fine_l, indicator = msm.invalid_phaserange(level), 0
            cells[msm.signal_id(o.signal)] = (fine_p, fine_l, indicator)
        return _Satellite(msm.satellite_id(first.sat), rough, msm.rates(first.rate), first.channel, cells)


def _groups(satellites: list[_Satellite]) -> list[list[_Satellite]]:
    """Satellites in ID order cut into consecutive groups whose satellites × signals fit 64 cells.

    [spec: RTCM 10403.3 3.5.16.3.4, "splitting on satellite-by-satellite basis";
    Table 3.5-78, X ≤ 64]
    """
    groups: list[list[_Satellite]] = []
    current: list[_Satellite] = []
    signals: set[int] = set()
    for s in satellites:
        joined = signals | set(s.cells)
        if current and (len(current) + 1) * len(joined) > MAX_CELLS:
            groups.append(current)
            current, joined = [], set(s.cells)
        current.append(s)
        signals = joined
    if current:
        groups.append(current)
    return groups


def _message(system: System, level: int, epoch: dict, satellites: list[_Satellite]) -> rtcm3.Msm:
    """One MSM of ``satellites``.

    [spec: RTCM 10403.3 Table 3.5-78 (DF394 MSb is ID 1, DF395 MSb is ID 1, DF396 by
    satellite then signal), Tables 3.5-79 to 3.5-88 (internal looping); DF409 0,
    DF411 1, DF412 0, DF417 0, DF418 0 (DESIGN 2.6.2); DF420 0, DF403/DF408 0
    (DESIGN 2.6.3)]
    """
    sig_ids = sorted({i for s in satellites for i in s.cells})
    sat_mask = sum(1 << (64 - s.sat_id) for s in satellites)
    sig_mask = sum(1 << (32 - i) for i in sig_ids)
    cell_mask = 0
    fine_range, fine_phase, locks = [], [], []
    for s in satellites:
        for i in sig_ids:
            cell_mask <<= 1
            cell = s.cells.get(i)
            if cell is None:
                continue
            cell_mask |= 1
            fine_range.append(cell[0])
            fine_phase.append(cell[1])
            locks.append(cell[2])
    n_cells = len(fine_range)
    rates = msm.has_rates(level)
    return rtcm3.Msm(
        number=msm.PREFIX[system] * 10 + level, station=0, **epoch, multiple=1, iods=0, reserved=0,
        clock_steering=1, external_clock=0, smoothing=0, smoothing_interval=0,
        sat_mask=sat_mask, sig_mask=sig_mask, cell_mask=cell_mask,
        rough_ms=tuple(s.rough[0] for s in satellites),
        extended_info=tuple(msm.frequency_channel(s.channel) if system == System.GLONASS else 0
                            for s in satellites) if rates else None,
        rough_mod=tuple(s.rough[1] for s in satellites),
        rough_rate=tuple(s.rate[0] for s in satellites) if rates else None,
        fine_range=tuple(fine_range), fine_phase=tuple(fine_phase), lock=tuple(locks),
        half_cycle=(0,) * n_cells, cnr=(0,) * n_cells,
        fine_rate=tuple(s.rate[1] for s in satellites for i in sig_ids if i in s.cells) if rates else None)
