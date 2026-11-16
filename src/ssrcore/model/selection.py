"""Selection of one consistent set of data for a query, and the transmit time.

DESIGN.md section 3.4. For a satellite and a reception time, ``choose`` takes
steps 1 to 3: the satellite stream (S10: the named one, or the active one of
S2), the newest usable orbit and clock of one
issue with the biases of the same issue (``Corrections.select``, S3, C3 to
C6), and the ephemeris that orbit block names (E1, E2, E5, validity). When the
corrections cannot be applied, the satellite still has its broadcast state, of
the nearest valid issue (DESIGN 2.4, ``satellite.position``). A ``Choice``
carries the interval of reception times over which it cannot change while
nothing is fed (DESIGN 2.8), so that the engine need not choose again.

``transmit`` takes step 4: the light-time iteration of DESIGN 2.3, with the
broadcast orbit corrected by C1 (SPARTN 8.8), and an orbit correction with
rates evaluated at each step (C11).

Everything here returns model types or a ``Missing`` reason; the record is
the engine's (DESIGN 3.1).
"""

from __future__ import annotations

from dataclasses import dataclass

from ssrcore.model.corrections import Corrections, Selection
from ssrcore.model.ephemeris import (
    EphemerisStore,
    StoredEphemeris,
    glonass_orbit,
    healthy,
    kepler_clock,
    kepler_orbit,
    own_scale,
    reference_time,
    validity,
)
from ssrcore.model.stream import Streams
from ssrcore.model.types import (
    SPEED_OF_LIGHT,
    GlonassEphemeris,
    GnssTime,
    MaxAge,
    Missing,
    OrbitCorrection,
    Sat,
    StreamId,
)
from ssrcore.toolkit import KernelError, Vector, geodesy, orbits, relativity, ssr

# [design: DESIGN 2.3] the first guess of the light time, s: GPS, Galileo and
# BeiDou MEO ranges of 20 000 to 26 000 km are 0.067 to 0.087 s, GEO and IGSO
# up to 0.13 s, and the iteration contracts by |ρ̇|/c (under 1e-5) per step.
FIRST_GUESS = 0.075
# [design: DESIGN 2.3] 1 ps, 0.3 mm of range, on the change of the light time.
TOLERANCE = 1e-12
# [design: DESIGN 2.3] from the first guess the change falls below 1 ps at the
# fourth step for every Earth orbit; one more is margin.
MAX_ITERATIONS = 5


@dataclass(frozen=True, slots=True)
class Broadcast:
    """A stored ephemeris with the orbit and clock arguments of its kernel prepared once.

    [design: DESIGN 2.8, the conversion of the elements (semicircles, own
    time scale) is done once per choice, not once per light-time step]
    """

    stored: StoredEphemeris
    orbit: orbits.KeplerOrbit | orbits.GlonassOrbit
    clock: orbits.KeplerClock | None

    @property
    def glonass(self) -> bool:
        """[spec: GLONASS ICD 4.4] whether the ephemeris is GLONASS's (state vectors, no Kepler clock)."""
        return self.clock is None

    def state(self, t: GnssTime) -> tuple[orbits.OrbitState, float]:
        """[spec: the orbit kernels, ``toolkit.orbits`` (DESIGN 2.7)] broadcast position,
        velocity and clock (s, without dt_rel and group delays) at GPS time ``t``;
        raises ``KernelError`` when the kernel refuses (as ``ephemeris.broadcast_state``)."""
        if self.clock is None:
            return (orbits.glonass_state(self.orbit, t.week, t.seconds),
                    orbits.glonass_clock(self.orbit, t.week, t.seconds))
        week, sow = own_scale(self.stored.ephemeris.sat.system, t)
        return orbits.kepler_state(self.orbit, week, sow), orbits.kepler_clock(self.clock, week, sow)

    def drift(self, t: GnssTime) -> float:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; GLONASS ICD 3.3.3, 4.4]
        broadcast clock drift (s/s) at ``t``."""
        if self.clock is None:
            return orbits.glonass_clock_drift(self.orbit)
        week, sow = own_scale(self.stored.ephemeris.sat.system, t)
        return orbits.kepler_clock_drift(self.clock, week, sow)


def prepare(stored: StoredEphemeris) -> Broadcast:
    """[spec: DESIGN 2.7, ``ephemeris.kepler_orbit``, ``kepler_clock``, ``glonass_orbit``]
    The kernel arguments of a stored ephemeris."""
    eph = stored.ephemeris
    if isinstance(eph, GlonassEphemeris):
        return Broadcast(stored, glonass_orbit(eph), None)
    return Broadcast(stored, kepler_orbit(eph), kepler_clock(eph))


@dataclass(frozen=True, slots=True)
class Choice:
    """Steps 1 to 3 of DESIGN 3.4 for one satellite, made at time ``t``.

    ``selection`` is the corrections' choice (orbit, clock and biases of one
    issue) or None; ``corrections_missing`` says why the orbit and clock
    corrections are not applied, None when they are. ``broadcast`` is the
    ephemeris the state is computed from: the one the orbit block names when
    the corrections apply, otherwise the nearest valid healthy issue, or
    None with ``broadcast_missing``. The choice holds for reception times from
    ``t`` up to, not including, ``until`` and up to and including ``through``
    (None: no bound). [design: DESIGN 2.8]
    """

    sat: Sat
    t: GnssTime
    stream: StreamId | None
    selection: Selection | None
    corrections_missing: Missing | None
    broadcast: Broadcast | None
    broadcast_missing: Missing | None
    until: GnssTime | None = None
    through: GnssTime | None = None

    def holds(self, t: GnssTime) -> bool:
        """[design: DESIGN 2.8] whether ``choose`` at ``t`` gives this same choice, nothing having been fed."""
        return (self.t <= t and (self.until is None or t < self.until)
                and (self.through is None or t <= self.through))


def _earliest(*times: GnssTime | None) -> GnssTime | None:
    present = [t for t in times if t is not None]
    return min(present) if present else None


def _no_broadcast(ephemerides: EphemerisStore, sat: Sat) -> Missing:
    """[design: DESIGN 2.4, the most specific reason: nothing stored, only
    unhealthy issues (E5), or none valid at the time (2.5)]"""
    issues = ephemerides.issues(sat)
    if not issues:
        return Missing.NO_EPHEMERIS
    if not any(healthy(s.ephemeris) for s in issues):
        return Missing.UNHEALTHY
    return Missing.OUTSIDE_VALIDITY


def choose(streams: Streams, corrections: Corrections, ephemerides: EphemerisStore,
           sat: Sat, t: GnssTime, max_age: MaxAge, stream: StreamId | None = None) -> Choice:
    """DESIGN 3.4 steps 1 to 3 for ``sat`` at reception time ``t``.

    1. The satellite stream (S10): ``stream`` when given, the active stream
       (S2) otherwise; without one, no corrections.
    2. ``Corrections.select``: the newest issue with usable orbit and clock
       (S3, C3, C4, C6), and its biases (C5).
    3. ``EphemerisStore.find`` with the orbit block's IODE and the OCB
       header's ephemeris type (E1, E2, E5, validity of 2.5).

    Without step 2 or 3 the state is broadcast only, from
    ``EphemerisStore.broadcast``. [spec: SPARTN 8.3, 6.2; design: DESIGN 3.4,
    2.4 for the broadcast state; 2.8 for the interval the choice holds over:
    the next stored block of the satellite, the end of each chosen block's age
    limit and of the ephemeris validity; a broadcast-only choice holds at
    ``t`` alone, since the nearest issue changes with time; S10, a named stream
    is not replaced by another]
    """
    if stream is None:
        stream = streams.active
    selection: Selection | None = None
    missing: Missing | None = Missing.NO_CORRECTIONS
    until = through = None
    if stream is not None:
        found = corrections.select(stream, sat, t, max_age)
        until = corrections.next_change(stream, sat, t)
        if isinstance(found, Selection):
            selection, missing = found, None
        else:
            missing = found
    if selection is not None:
        orbit = selection.orbit
        lookup = ephemerides.find(sat, orbit.value.iode, orbit.value.ephemeris_type, t)
        if lookup.stored is not None:
            biases = selection.biases
            through = _earliest(orbit.time + max_age.orbit, selection.clock.time + max_age.clock,
                                None if biases is None else biases.time + max_age.bias,
                                reference_time(lookup.stored.ephemeris) + validity(lookup.stored.ephemeris))
            return Choice(sat, t, stream, selection, None, prepare(lookup.stored), None, until, through)
        missing = lookup.miss
    stored = ephemerides.broadcast(sat, t)
    if stored is None:
        return Choice(sat, t, stream, selection, missing, None, _no_broadcast(ephemerides, sat), until, t)
    return Choice(sat, t, stream, selection, missing, prepare(stored), None, until, t)


@dataclass(frozen=True, slots=True)
class Transmit:
    """The satellite at transmission (DESIGN 3.4 step 4).

    ``t_tx`` is the GPS time of transmission; ``state`` and ``clock`` the
    broadcast position, velocity and clock there; ``correction`` δX (ECEF, m,
    None without an orbit correction) and ``position`` X = X_b − δX; ``range``
    ρ = |X − r_r| and its rate; ``sagnac`` E, so that (ρ + E)/c is the light
    time to within ``TOLERANCE``; ``iterations`` the steps taken; ``orbit``
    the δO (radial, along, cross, m) that gave δX: as sent, or at ``t_tx``
    for a correction with rates (rule C11).
    """

    t_tx: GnssTime
    state: orbits.OrbitState
    clock: float
    correction: Vector | None
    position: Vector
    range: float
    range_rate: float
    sagnac: float
    iterations: int
    orbit: Vector | None = None


def transmit(broadcast: Broadcast, orbit: OrbitCorrection | None, receiver: Vector,
             t_rx: GnssTime) -> Transmit | Missing:
    """The transmit time and the satellite there, by the light-time iteration.

    From τ = ``FIRST_GUESS``, each step evaluates the broadcast orbit at
    t_tx = t_rx − τ (a ``GnssTime``, so that t_rx − t_tx is exact, DESIGN 4),
    corrects it (C1), and sets τ' = (ρ + E)/c with ρ = |X − r_r| and E the
    Sagnac term of the same positions; it stops when |τ' − τ| < ``TOLERANCE``
    and returns the state at that t_tx. A kernel that refuses its input, or no
    convergence within ``MAX_ITERATIONS``, gives ``not_converged``.
    [spec: Ashby 2003 section 2, the light path is straight in the inertial
    frame, c τ = |R_z(ω_E τ) X(t_rx − τ) − r_r|; SPARTN 8.8 for X; RTCM
    10403.3 3.5.13.5 eq. (3.12-6) for an orbit correction with rates, evaluated
    at each step's t_tx (C11); design: DESIGN 2.3 for the first guess, the
    tolerance and the count]
    """
    sent = d_o = None if orbit is None else (orbit.radial, orbit.along, orbit.cross)
    rate = None if orbit is None else orbit.rate
    tau = FIRST_GUESS
    for step in range(1, MAX_ITERATIONS + 1):
        t_tx = t_rx - tau
        try:
            state, clock = broadcast.state(t_tx)
            if rate is not None:
                # [spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6); rule C11, δO at this step's t_tx]
                d_o = ssr.orbit_at(sent, rate, t_tx - orbit.reference)
            position = state.position if d_o is None else ssr.corrected_position(state.position, state.velocity, d_o)
            rho, range_rate = geodesy.range_rate(position, state.velocity, receiver)
            sagnac = relativity.sagnac(position, receiver)
        except KernelError:
            return Missing.NOT_CONVERGED
        light = (rho + sagnac) / SPEED_OF_LIGHT
        if abs(light - tau) < TOLERANCE:
            correction = None if d_o is None else ssr.orbit_correction(state.position, state.velocity, d_o)
            return Transmit(t_tx, state, clock, correction, position, rho, range_rate, sagnac, step, d_o)
        tau = light
    return Missing.NOT_CONVERGED
