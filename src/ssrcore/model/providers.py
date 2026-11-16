"""Atmosphere providers: the seam between the engine and whatever gives it atmosphere.

DESIGN 3.6 is the design. The engine needs two atmosphere terms of
a query: the zenith troposphere at the receiver and the slant ionosphere of
a satellite. Each comes from an ordered list of providers; per term, the
first provider that answers wins (``walk_troposphere``,
``walk_ionosphere``). A provider is any object with the members of
``AtmosphereProvider``:

- ``name`` (unique in the list), ``kind`` (``Kind.CORRECTION`` or
  ``Kind.MODEL``) and ``terms``, the terms it can ever give. A term not in
  ``terms`` is one the provider will *never* have; the walk does not ask it.
- ``troposphere(where)`` and ``ionosphere(where, sat, look)``, each returning a
  value or ``Absent(reason)``: a term the provider has in general but not
  for this query, *not yet* or *not here*.

Three kinds are here: ``SpartnAtmosphere`` presents the SPARTN atmosphere
store through the interface without changing it (its area choice A4 and
its issue-of-update binding A1 stay inside it); ``Unb3mTroposphere`` is the
a-priori model of ``toolkit.unb3m``, labelled a model; anything else the
caller writes against the interface (an outside provider).

Everything here returns model types; the record is the engine's (DESIGN 3.1).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from functools import lru_cache
from typing import Protocol, runtime_checkable

from ssrcore.model.atmosphere import Absent, AreaChoice, AtmosphereStore, IonosphereTec, TroposphereZenith
from ssrcore.model.stream import Issue
from ssrcore.model.types import (
    CallerSource,
    Config,
    GnssTime,
    Missing,
    ModelSource,
    Sat,
    SpartnSource,
    StreamId,
    Term,
)
from ssrcore.toolkit import KernelError, ssr, troposphere, unb3m

# [spec: SPARTN Table 6.2 SF042, SF055] the quality fields, 1 sigma at zenith (m) and slant (TECU).
_SF_TROPOSPHERE_QUALITY = 42
_SF_IONOSPHERE_QUALITY = 55

ATMOSPHERE_TERMS = frozenset({Term.TROPOSPHERE, Term.IONOSPHERE})


class Kind(StrEnum):
    """What a value is. [design: DESIGN 3.6, a correction (a measurement of this
    atmosphere) and a model (a climatology) differ in what their error is, so the
    record says which it gave]"""

    CORRECTION = "correction"
    MODEL = "model"


@dataclass(frozen=True, slots=True)
class Origin:
    """Which provider gave a term, what kind of value it is, and from which SPARTN
    stream when it came from one (None otherwise).
    [design: DESIGN 3.6, interface]"""

    provider: str
    kind: Kind
    stream: StreamId | None = None


@dataclass(frozen=True, slots=True)
class Binding:
    """The stream and SIOU issue of the satellite corrections a query uses (DESIGN
    3.4 step 2), which rule A1 binds the same stream's atmosphere to.
    [spec: SPARTN 6.2, SF005 combines parts sent in separate messages]"""

    stream: StreamId
    issue: Issue


@dataclass(slots=True)
class Where:
    """What a provider is asked about: the query time ``t`` (GPS time), the
    receiver's geodetic position (``latitude``, ``longitude`` in radians and
    in degrees, ``height`` ellipsoidal in metres), the day of year counted from
    January 0.0 (DESIGN 2.7), and the satellite corrections' ``binding`` (None,
    with the selection's reason in ``unbound``, when the query has none).

    ``cache`` is the engine's per-epoch store of values that are functions of
    the stores and of these arguments alone (DESIGN 2.8); a provider may
    keep such values there under keys of its own.
    """

    t: GnssTime
    latitude: float
    longitude: float
    latitude_deg: float
    longitude_deg: float
    height: float
    doy: float
    binding: Binding | None = None
    unbound: Missing | None = None
    cache: dict = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Look:
    """The direction from the receiver to a satellite, as the record's geometry has it:
    ``elevation`` and ``azimuth`` in radians, azimuth clockwise from north (DESIGN 2.4).
    A provider whose slant value depends on the path (a global map evaluated at a
    pierce point) needs it; SPARTN's per-satellite polynomial does not."""

    elevation: float
    azimuth: float


@dataclass(frozen=True, slots=True)
class Zenith:
    """A zenith troposphere at the receiver.

    ``hydrostatic`` and ``wet`` are zenith delays in metres at the receiver's
    height, mapped by the engine with Niell m_h and m_w (SPARTN 8.12, UNB3m
    eq. (8)). ``time`` is the time the value refers to and ``age`` the query
    time less it (None for a model, which has neither). ``sigma`` is the
    provider's 1-sigma of the total zenith delay in metres, None when it
    states none or only a bound below. ``spartn`` is the SPARTN store's own
    answer when the provider is SPARTN's, so the engine records the
    polynomial, grid, quality and area exactly as before.
    """

    hydrostatic: float
    wet: float
    origin: Origin
    time: GnssTime | None = None
    age: float | None = None
    sigma: float | None = None
    spartn: TroposphereZenith | None = None


@dataclass(frozen=True, slots=True)
class SlantTec:
    """A slant ionosphere of one satellite at the receiver, TECU (the engine turns it
    into metres per signal, SPARTN 8.11). Fields as ``Zenith``; ``sigma`` in TECU."""

    sat: Sat
    tec: float
    origin: Origin
    time: GnssTime | None = None
    age: float | None = None
    sigma: float | None = None
    spartn: IonosphereTec | None = None


@runtime_checkable
class AtmosphereProvider(Protocol):
    """The interface every atmosphere provider presents. [design: DESIGN 3.6]"""

    name: str
    kind: Kind
    terms: frozenset[Term]

    def troposphere(self, where: Where) -> Zenith | Absent:
        """[design: DESIGN 3.6] the zenith troposphere at ``where``, or why not."""
        ...

    def ionosphere(self, where: Where, sat: Sat, look: Look) -> SlantTec | Absent:
        """[design: DESIGN 3.6] the slant TEC of ``sat``, seen along ``look``, at
        ``where``, or why not."""
        ...


def check(providers: Iterable[AtmosphereProvider]) -> tuple[AtmosphereProvider, ...]:
    """The providers as a tuple, each with the interface, a name used once, a ``Kind``
    and atmosphere terms only; raises TypeError or ValueError otherwise.
    [design: DESIGN 3.6, a list the engine can walk; a name says in the
    record which provider answered, so two providers may not share one]"""
    out = tuple(providers)
    names = set()
    for p in out:
        if not isinstance(p, AtmosphereProvider):
            raise TypeError(f"not an atmosphere provider: {type(p).__name__}")
        if not isinstance(p.kind, Kind):
            raise TypeError(f"provider {p.name!r}: kind must be a Kind, not {p.kind!r}")
        if not p.terms <= ATMOSPHERE_TERMS:
            raise ValueError(f"provider {p.name!r}: terms {sorted(p.terms - ATMOSPHERE_TERMS)} are not atmosphere")
        if p.name in names:
            raise ValueError(f"provider name {p.name!r} used twice")
        names.add(p.name)
    return out


_NOBODY = Absent(Missing.NO_ATMOSPHERE)


def walk_troposphere(providers: tuple[AtmosphereProvider, ...], where: Where) -> Zenith | Absent:
    """The zenith troposphere of the first provider that answers, in list order, or the
    ``Absent`` of the first provider asked (``no_atmosphere`` when no provider has the
    term). The providers passed over are not kept: the answer's ``origin`` says who
    answered, the reason why nobody did.
    [design: DESIGN 3.6 and rule A9: the caller ranks the providers; per
    term, the first answer wins; a provider whose ``terms`` lack the term is not
    asked. With the default list of one SPARTN provider the value and the reason
    are the store's. DESIGN 2.5: no list of the providers passed over]"""
    first = None
    for p in providers:
        if Term.TROPOSPHERE in p.terms:
            value = p.troposphere(where)
            if not isinstance(value, Absent):
                return value
            if first is None:
                first = value
    return _NOBODY if first is None else first


def walk_ionosphere(providers: tuple[AtmosphereProvider, ...], where: Where, sat: Sat,
                    look: Look) -> SlantTec | Absent:
    """The slant TEC of ``sat`` from the first provider that answers, in list order, or
    why not, as ``walk_troposphere``. [design: DESIGN 3.6 and rule A9]"""
    first = None
    for p in providers:
        if Term.IONOSPHERE in p.terms:
            value = p.ionosphere(where, sat, look)
            if not isinstance(value, Absent):
                return value
            if first is None:
                first = value
    return _NOBODY if first is None else first


@lru_cache(maxsize=32)
def _sigma(sf: int, index: int) -> float | None:
    """[spec: SPARTN Table 6.2 SF042, SF055, "<=" each step, 0 unknown, the last step a
    bound below] the 1-sigma a quality index states, None when it states none.
    [design: a function of eight indices per field, kept once evaluated]"""
    value = ssr.field(sf, index)
    return value.value if value.status == ssr.Status.OK else None


class SpartnAtmosphere:
    """The SPARTN atmosphere store presented as a provider, its numerics untouched.

    ``stream`` None (the default) takes the stream and SIOU issue of the
    query's satellite corrections (``Where.binding``): the atmosphere of the
    stream that gives orbit, clock and biases, under their SIOU (rule A1),
    which is the engine's behaviour before the seam. A named ``stream`` takes
    that stream's atmosphere: under the binding's SIOU issue when it is also
    the satellite stream (A1 again), otherwise under its own current SIOU
    issue (rule A10), since an SIOU compares only within one stream (S1).

    The area choice (A4) and the issue binding are the store's and stay here;
    the hydrostatic delay is brought to the receiver height by the UNB3m
    factor when ``height_reduction`` is set (DESIGN Q7, as
    ``Config.troposphere_height_reduction``).
    """

    kind = Kind.CORRECTION
    terms = ATMOSPHERE_TERMS

    def __init__(self, store: AtmosphereStore, stream: StreamId | None = None,
                 height_reduction: bool = True, name: str | None = None) -> None:
        """[design: DESIGN 3.6, the SPARTN provider; DESIGN Q7 for the reduction]"""
        self.store = store
        self.stream = stream
        self.height_reduction = bool(height_reduction)
        if name is None:
            name = "spartn" if stream is None else f"spartn:{stream.solution}/{stream.processor}"
        self.name = name
        self._origins: dict[StreamId, Origin] = {}

    def _issue(self, where: Where) -> tuple[StreamId, Issue] | Absent:
        """[spec: SPARTN 6.2, SF005 (rule A1); design: rule A10 for another stream]"""
        b = where.binding
        if self.stream is None:
            if b is None:
                return Absent(where.unbound or Missing.NO_CORRECTIONS)
            return b.stream, b.issue
        if b is not None and b.stream == self.stream:
            return b.stream, b.issue
        state = self.store.streams.state(self.stream)
        if state is None or state.siou.value is None:
            return Absent(Missing.NO_ATMOSPHERE)
        return self.stream, Issue(state.siou.value, state.siou.era)

    def area(self, where: Where) -> AreaChoice | Absent:
        """The store's area for the query (rules A1, A2, A4), once per epoch and issue.
        [spec: SPARTN 8.10 through ``AtmosphereStore.area``; design: DESIGN 2.8]"""
        issue = self._issue(where)
        if isinstance(issue, Absent):
            return issue
        key = ("spartn.area", id(self.store), issue)
        if key not in where.cache:
            where.cache[key] = self.store.area(issue[0], issue[1], where.latitude_deg, where.longitude_deg,
                                               where.t)
        return where.cache[key]

    def troposphere(self, where: Where) -> Zenith | Absent:
        """[spec: SPARTN 8.12 through ``AtmosphereStore.troposphere``; UNB3m eqs. (4), (6) for the
        height factor (DESIGN Q7)] the store's zenith troposphere, the hydrostatic part at
        the receiver height."""
        area = self.area(where)
        if isinstance(area, Absent):
            return area
        # [design: DESIGN 2.8, a function of the store and the epoch alone, so the
        # answer, height factor included, is kept once per epoch and area]
        key = ("spartn.zenith", id(self), area)
        answer = where.cache.get(key)
        if answer is None:
            answer = zenith = self.store.troposphere(area, where.t)
            if not isinstance(zenith, Absent):
                hydrostatic = zenith.hydrostatic
                if self.height_reduction:
                    hydrostatic = hydrostatic * troposphere.hydrostatic_height_factor(where.latitude, where.height,
                                                                                       where.doy)
                answer = Zenith(hydrostatic, zenith.wet, self._origin(area.stream), zenith.time,
                                zenith.age, _sigma(_SF_TROPOSPHERE_QUALITY, zenith.quality), zenith)
            where.cache[key] = answer
        return answer

    def ionosphere(self, where: Where, sat: Sat, look: Look) -> SlantTec | Absent:
        """[spec: SPARTN 8.11 through ``AtmosphereStore.ionosphere``; rules A1, A6, A7] the
        store's slant TEC of ``sat`` in the area of the query; ``look`` is not needed,
        the polynomial being per satellite at the receiver."""
        area = self.area(where)
        if isinstance(area, Absent):
            return area
        tec = self.store.ionosphere(area, sat, where.t)
        if isinstance(tec, Absent):
            return tec
        return SlantTec(sat, tec.slant, self._origin(area.stream), tec.time, tec.age,
                        _sigma(_SF_IONOSPHERE_QUALITY, tec.quality), tec)

    def _origin(self, stream: StreamId) -> Origin:
        """[design: DESIGN 3.6, the record's origin; one value per stream, kept]"""
        origin = self._origins.get(stream)
        if origin is None:
            origin = self._origins[stream] = Origin(self.name, self.kind, stream)
        return origin


class Unb3mTroposphere:
    """The UNB3m a-priori troposphere (``toolkit.unb3m``), for when no correction stream
    gives troposphere; its value is labelled a model.

    UNB3m's height is orthometric. ``undulation`` (latitude, longitude in
    radians to the geoid height N in metres) turns the receiver's ellipsoidal
    height into it; without one the ellipsoidal height is used. [design: no
    geoid model is part of ssrcore; |N| up to about 100 m moves the zenith
    delay by about 3 cm, inside the model's 4.9 cm, and a caller that has a
    geoid passes it]
    """

    kind = Kind.MODEL
    terms = frozenset({Term.TROPOSPHERE})

    def __init__(self, undulation: Callable[[float, float], float] | None = None, name: str = "unb3m") -> None:
        """[design: DESIGN 3.6, the a-priori troposphere]"""
        self.undulation = undulation
        self.name = name

    def troposphere(self, where: Where) -> Zenith | Absent:
        """[spec: UNB3m eqs. (1) to (11); Table 7 for sigma; design: a height beyond the
        lapse-rate profile is refused by the kernel, which DESIGN 3.4 names not_converged]"""
        height = where.height
        if self.undulation is not None:
            height = height - self.undulation(where.latitude, where.longitude)
        try:
            z = unb3m.zenith_delays(where.latitude, height, where.doy)
        except KernelError:
            return Absent(Missing.NOT_CONVERGED)
        return Zenith(z.hydrostatic, z.wet, Origin(self.name, self.kind), sigma=unb3m.TOTAL_SIGMA)

    def ionosphere(self, where: Where, sat: Sat, look: Look) -> SlantTec | Absent:
        """[design: UNB3m is a neutral-atmosphere model; ``terms`` says the provider never
        has ionosphere, so the walk never asks]"""
        return Absent(Missing.NO_ATMOSPHERE)


def build(config: Config, store: AtmosphereStore,
          callers: Mapping[str, AtmosphereProvider] | None = None) -> tuple[AtmosphereProvider, ...]:
    """The providers of ``config.atmosphere``, in its order: ``SpartnSource`` the SPARTN
    store (bound to the satellite stream when its stream is None), ``ModelSource`` the
    UNB3m troposphere, ``CallerSource`` the caller's object of that name, which must
    carry the same ``name``. Every ``CallerSource`` must have its object and every
    object its ``CallerSource``; then ``check``. Raises ValueError or TypeError.
    [design: DESIGN 2.1, 2.5, rule A9; the configuration names the sources and stays a
    frozen value, the engine builds the objects once]"""
    callers = dict(callers or {})
    named = [s.name for s in config.atmosphere if isinstance(s, CallerSource)]
    if set(named) != set(callers):
        raise ValueError(f"caller providers {sorted(callers)} do not match the CallerSources {sorted(named)}")
    out: list = []
    for source in config.atmosphere:
        if isinstance(source, SpartnSource):
            out.append(SpartnAtmosphere(store, source.stream, config.troposphere_height_reduction))
        elif isinstance(source, ModelSource):
            out.append(Unb3mTroposphere())
        else:
            provider = callers[source.name]
            if not isinstance(provider, AtmosphereProvider):
                raise TypeError(f"not an atmosphere provider: {type(provider).__name__}")
            if provider.name != source.name:
                raise ValueError(f"the provider passed as {source.name!r} is named {provider.name!r}")
            out.append(provider)
    return check(out)


def default_providers(store: AtmosphereStore, height_reduction: bool = True) -> tuple[AtmosphereProvider, ...]:
    """[design: DESIGN 3.6 and rule A9, the default: the SPARTN store bound to the satellite
    stream, alone, which is the engine's behaviour before the seam]"""
    return (SpartnAtmosphere(store, None, height_reduction),)
