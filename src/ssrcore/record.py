"""The line-of-sight record and its parts, every term separate.

DESIGN.md section 2.4. Every term is in SI units and carries the sign it has
in the observation equations the record is built for (metres; ``s`` a signal
of wavelength λs):

    P_s = ρ + E + S − c·(dt_sat + dt_rel) + T + I_s + b_s       + R_tide
    Φ_s = ρ + E + S − c·(dt_sat + dt_rel) + T − I_s + B_s + λs·w + R_tide + λs·N

(receiver terms omitted). A term that is not valid is absent (None, or not in
its mapping) and ``missing`` says why; it is never silently zero (DESIGN 2.3).
The engine fills the record; this module only holds it and forms
the two sums.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

# [design: DESIGN 3.1, Term and Missing are defined once, in the model's types,
# so the stores can give the record's reasons; re-exported here]
from ssrcore.model.providers import Origin
from ssrcore.model.types import SPEED_OF_LIGHT, GnssTime, Missing, Sat, Signal, Term

Vector = tuple[float, float, float]


def _frozen_mapping(value: Mapping) -> Mapping:
    return MappingProxyType(dict(value))


class MissingTermError(LookupError):
    """A term a sum needs is absent. [design: DESIGN 2.4, the sums fail rather than skip a term]"""

    def __init__(self, term: Term, signal: Signal | None, reason: Missing | None):
        self.term, self.signal, self.reason = term, signal, reason
        where = term.value if signal is None else f"{term.value}[{signal.code}]"
        super().__init__(f"{where} is missing" + (f": {reason.value}" if reason else ""))


# ---- parts ----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Geometry:
    """Receiver to satellite: ``range`` ρ in metres, from the receiver to the
    corrected antenna phase centre at ``t_tx``, each in ECEF of its own epoch;
    ``range_rate`` ρ̇ = u·ṙ_s (m/s) for a receiver at rest in ECEF; ``unit``
    the unit vector towards the satellite in the ECEF frame of reception,
    ``azimuth`` and ``elevation`` in radians in the WGS-84 local frame (TM-23
    3.1.4; Ashby 2003 section 2; DESIGN 2.4)."""

    range: float
    range_rate: float
    unit: Vector
    azimuth: float
    elevation: float


@dataclass(frozen=True, slots=True)
class Modelled:
    """The terms ssrcore models, each None when absent.

    ``earth_rotation`` E (m) [spec: Ashby 2003 section 2]; ``shapiro`` S (m)
    [spec: Ashby 2003 section 10; SPARTN 8.16]; ``relativity_clock`` dt_rel
    (s), −2 r·ṙ/c², zero for GLONASS [spec: SPARTN 8.7]; ``solid_tide``
    R_tide = −u·Δr (m) and ``tide_displacement`` Δr (ECEF, m) [spec: TN36
    7.1.1]; ``windup`` w in cycles, in (−½, ½] [spec: TM-23 5.5; SPARTN 8.14].
    """

    earth_rotation: float | None = None
    shapiro: float | None = None
    relativity_clock: float | None = None
    solid_tide: float | None = None
    tide_displacement: Vector | None = None
    windup: float | None = None


@dataclass(frozen=True, slots=True)
class Broadcast:
    """The broadcast ephemeris used: position (ECEF, m) and clock (s) at
    ``t_tx`` before correction, its issue of data, the RTCM message it came
    from (its type, rule E2), its toe (tb for GLONASS) and age at ``t_rx`` (s)."""

    position: Vector
    clock: float
    iod: int
    message: int
    toe: GnssTime
    age: float


@dataclass(frozen=True, slots=True)
class SatelliteState:
    """The satellite at ``t_tx``.

    ``position`` X = X_b − δX (ECEF, m) [spec: SPARTN 8.8]; ``velocity`` from
    the broadcast orbit (m/s); ``clock`` dt_sat = t_b − δC/c (s) [spec: SPARTN
    8.7]; ``clock_drift`` the broadcast drift (s/s), since SPARTN carries no
    rate, less (C1 + 2·C2·(t_tx − t0))/c for a clock polynomial (rule C11),
    and without the rate of dt_rel (DESIGN 2.4); ``frequency_channel``
    the GLONASS channel, None elsewhere.
    """

    position: Vector
    velocity: Vector
    clock: float
    clock_drift: float
    broadcast: Broadcast
    frequency_channel: int | None = None


@dataclass(frozen=True, slots=True)
class OrbitTerm:
    """The orbit correction used: δO radial, along, cross as sent (or, with
    rates, at ``t_tx``, rule C11) and δX in ECEF (m), the IODE it applies to,
    yaw in degrees (None when absent, rule C8), the time tag (GPS time) and
    the age at ``t_rx`` (s); ``rate`` δȮ (m/s) and ``reference`` t0 as sent,
    None for SPARTN."""

    radial: float
    along: float
    cross: float
    ecef: Vector
    iode: int
    yaw: float | None
    time: GnssTime
    age: float
    rate: Vector | None = None
    reference: GnssTime | None = None


@dataclass(frozen=True, slots=True)
class ClockTerm:
    """The clock correction used: δC as sent (or, for a polynomial, at
    ``t_tx``, rule C11) (m), IODE continuity SF022 (s), URE SF024 index and its
    value in metres (None for the index 0, unknown, and 7, above 1 m), time
    tag and age (s); ``rate`` C1 (m/s), ``acceleration`` C2 (m/s²) and
    ``reference`` t0 as sent, None for SPARTN."""

    correction: float
    iode_continuity: float
    ure_index: int
    ure: float | None
    time: GnssTime
    age: float
    rate: float | None = None
    acceleration: float | None = None
    reference: GnssTime | None = None


@dataclass(frozen=True, slots=True)
class CodeBiasTerm:
    """b_s: the code bias the raw pseudorange contains (m), which SPARTN says
    to subtract [spec: SPARTN 8.9]; time tag and age (s)."""

    value: float
    time: GnssTime
    age: float


@dataclass(frozen=True, slots=True)
class PhaseBiasTerm:
    """B_s: the phase bias the raw phase contains (m) [spec: SPARTN 8.9], its
    fix flag SF023, continuity SF015 (s), the count of discontinuities of its
    arc (rule C7), time tag and age (s)."""

    value: float
    fixed: bool
    continuity: float
    discontinuities: int
    time: GnssTime
    age: float


@dataclass(frozen=True, slots=True)
class TroposphereTerm:
    """T and its parts [spec: SPARTN 8.12]: ``slant`` (m); the hydrostatic
    zenith delay at the receiver height (m); hydrostatic and wet mapping
    factors; the wet zenith polynomial and grid residual (m; ``grid`` None
    when absent, rule A7); quality SF042 index; area ID; age (s).

    ``wet_zenith`` is the wet zenith delay used (SPARTN T_P + T_G, or a
    model's), ``origin`` the provider that answered (rule A9) and ``sigma``
    its 1σ of the total zenith delay (m, None when it states none). For
    another provider than SPARTN, ``polynomial``, ``grid``, ``quality`` and
    ``area_id`` are None, and ``age`` is None for a model (DESIGN 2.4)."""

    slant: float
    hydrostatic_zenith: float
    hydrostatic_mapping: float
    wet_mapping: float
    polynomial: float | None
    grid: float | None
    quality: int | None
    area_id: int | None
    age: float | None
    wet_zenith: float | None = None
    origin: Origin | None = None
    sigma: float | None = None


@dataclass(frozen=True, slots=True)
class IonosphereTerm:
    """I_s [spec: SPARTN 8.11]: slant TEC in TECU from the polynomial and the
    grid residual (``grid`` None when absent, rule A7), and per requested
    signal the code delay +40.3·10¹⁶·TEC/f² in metres; the phase delay is its
    negative. Quality SF055 index, area ID, age (s).

    ``origin`` and ``sigma`` (TECU) as for ``TroposphereTerm``;
    ``with_biases`` True when the ionosphere came from the stream of the
    record's code biases, False when not (rule A11), None without code
    biases. For another provider than SPARTN, ``polynomial``, ``grid``,
    ``quality`` and ``area_id`` are None (DESIGN 2.4)."""

    slant_tec: float
    polynomial: float | None
    grid: float | None
    code_delay: Mapping[Signal, float]
    quality: int | None
    area_id: int | None
    age: float | None
    origin: Origin | None = None
    sigma: float | None = None
    with_biases: bool | None = None

    def __post_init__(self) -> None:
        """[design: the record is immutable throughout]"""
        object.__setattr__(self, "code_delay", _frozen_mapping(self.code_delay))

    def phase_delay(self, signal: Signal) -> float:
        """[spec: SPARTN 8.11, δI_L = −40.3·10¹⁶/f²·TEC, the code delay with opposite sign] metres."""
        return -self.code_delay[signal]


@dataclass(frozen=True, slots=True)
class StreamQuality:
    """Stream indicators (DESIGN 2.5): solution ID TF010, processor ID TF011,
    SIOU SF005, AIOU SF068 (None without atmosphere), datum SF009, whether
    the last OCB set ended with SF010 and whether a frame of it failed
    (DESIGN 3.5), the active-stream switches (S2) and datum changes (S6)
    so far, failed frames counted by (protocol, kind, message type,
    subtype) (S7), and notices likewise (S9)."""

    solution_id: int
    processor_id: int
    siou: int
    aiou: int | None
    datum: int
    end_of_set: bool
    failures: Mapping[tuple[str, str, int | None, int | None], int] = field(default_factory=dict)
    incomplete: bool = False
    switches: int = 0
    datum_changes: int = 0
    notices: Mapping[tuple[str, str, int | None, int | None], int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """[design: the record is immutable throughout]"""
        object.__setattr__(self, "failures", _frozen_mapping(self.failures))
        object.__setattr__(self, "notices", _frozen_mapping(self.notices))


@dataclass(frozen=True, slots=True)
class SatelliteQuality:
    """Satellite indicators (DESIGN 2.5): URE SF024 index, fix flag SF023,
    continuity SF015 (s), IODE continuity SF022 (s), yaw present, and the
    count of satellite discontinuities (rule C7)."""

    ure_index: int
    fixed: bool
    continuity: float
    iode_continuity: float
    yaw_present: bool
    discontinuities: int


@dataclass(frozen=True, slots=True)
class Quality:
    """The quality indicators of a record; atmosphere quality is in its terms."""

    stream: StreamQuality | None = None
    satellite: SatelliteQuality | None = None


# ---- the record -----------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class LineOfSight:
    """One satellite seen from one receiver at one time: every term separate.

    ``sat``, ``t_rx`` (GPS time of reception), ``receiver`` (ECEF, m) and
    ``signals`` are the query; ``t_tx`` the GPS time of transmission.
    ``frequency`` gives the carrier frequency (Hz) of each requested signal
    the engine could name one for. ``missing`` maps (term, signal or None) to
    the reason the term is absent.
    """

    sat: Sat
    t_rx: GnssTime
    receiver: Vector
    signals: tuple[Signal, ...] = ()
    t_tx: GnssTime | None = None
    geometry: Geometry | None = None
    satellite: SatelliteState | None = None
    orbit: OrbitTerm | None = None
    clock_correction: ClockTerm | None = None
    modelled: Modelled = field(default_factory=Modelled)
    code_bias: Mapping[Signal, CodeBiasTerm] = field(default_factory=dict)
    phase_bias: Mapping[Signal, PhaseBiasTerm] = field(default_factory=dict)
    troposphere: TroposphereTerm | None = None
    ionosphere: IonosphereTerm | None = None
    frequency: Mapping[Signal, float] = field(default_factory=dict)
    quality: Quality = field(default_factory=Quality)
    missing: Mapping[tuple[Term, Signal | None], Missing] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """[design: the record is immutable throughout, and names its terms and reasons by enumeration]"""
        object.__setattr__(self, "receiver", tuple(float(v) for v in self.receiver))
        object.__setattr__(self, "signals", tuple(self.signals))
        for name in ("code_bias", "phase_bias", "frequency"):
            object.__setattr__(self, name, _frozen_mapping(getattr(self, name)))
        object.__setattr__(self, "missing", _frozen_mapping(
            {(Term(t), s): Missing(r) for (t, s), r in self.missing.items()}))

    @property
    def available(self) -> bool:
        """[design: DESIGN 2.4] true when the satellite position and clock exist."""
        return self.satellite is not None

    @property
    def corrected(self) -> bool:
        """[spec: SPARTN 8.7, 8.8; design: DESIGN 2.4] true when the orbit and clock
        corrections are applied to the satellite state."""
        return self.satellite is not None and self.orbit is not None and self.clock_correction is not None

    def why(self, term: Term, signal: Signal | None = None) -> Missing | None:
        """[design: DESIGN 2.4] the reason a term (of a signal) is absent, or None."""
        return self.missing.get((Term(term), signal))

    def _need(self, value, term: Term, signal: Signal | None = None):
        if value is None:
            raise MissingTermError(term, signal, self.why(term, signal))
        return value

    def _common(self) -> float:
        """ρ + E + S − c·(dt_sat + dt_rel) + T + R_tide, the terms shared by both sums."""
        geometry = self._need(self.geometry, Term.GEOMETRY)
        satellite = self._need(self.satellite, Term.SATELLITE)
        m = self.modelled
        clock = satellite.clock + self._need(m.relativity_clock, Term.RELATIVITY_CLOCK)
        troposphere = self._need(self.troposphere, Term.TROPOSPHERE)
        return (geometry.range
                + self._need(m.earth_rotation, Term.EARTH_ROTATION)
                + self._need(m.shapiro, Term.SHAPIRO)
                - SPEED_OF_LIGHT * clock
                + troposphere.slant
                + self._need(m.solid_tide, Term.SOLID_TIDE))

    def _ionosphere(self, signal: Signal) -> float:
        ionosphere = self._need(self.ionosphere, Term.IONOSPHERE)
        return self._need(ionosphere.code_delay.get(signal), Term.IONOSPHERE, signal)

    def pseudorange_model(self, signal: Signal) -> float:
        """ρ + E + S − c·(dt_sat + dt_rel) + T + I_s + b_s + R_tide in metres, without
        receiver terms; raises ``MissingTermError`` when a term is absent.

        [spec: SPARTN 8.7 to 8.12, 8.16 for the terms and their signs; design:
        DESIGN 2.4, a convenience over the separate terms]
        """
        bias = self._need(self.code_bias.get(signal), Term.CODE_BIAS, signal)
        return self._common() + self._ionosphere(signal) + bias.value

    def phase_model(self, signal: Signal) -> float:
        """ρ + E + S − c·(dt_sat + dt_rel) + T − I_s + B_s + λs·w + R_tide in metres,
        without receiver terms and without the ambiguity λs·N; raises
        ``MissingTermError`` when a term is absent.

        [spec: SPARTN 8.7 to 8.12, 8.14, 8.16 for the terms and their signs;
        design: DESIGN 2.4, a convenience over the separate terms]
        """
        bias = self._need(self.phase_bias.get(signal), Term.PHASE_BIAS, signal)
        windup = self._need(self.modelled.windup, Term.WINDUP)
        frequency = self._need(self.frequency.get(signal), Term.FREQUENCY, signal)
        wavelength = SPEED_OF_LIGHT / frequency
        return self._common() - self._ionosphere(signal) + bias.value + wavelength * windup
