"""The IGS SSR global ionosphere (IM201) as an atmosphere provider.

IGS SSR v1.00 4.5.1 and 8.4.1. The message gives the vertical TEC of up to
four thin layers as spherical harmonic expansions; for a query the provider
finds each layer's pierce point along the line of sight in the 6370 km
sphere of the document, evaluates the expansion there, maps it to slant TEC
with 1/sin(E + ψ) and sums the layers (``toolkit.ionosphere``). The engine
turns the slant TEC into metres per signal as for any provider (TM-23 eq.
(5.40), the 40.3/f² of IGS SSR 4.5.1).

It presents the interface of ``model.providers`` unchanged: ``Where`` gives
the receiver and the time, ``Look`` the ellipsoidal elevation and azimuth of
the line of sight in the frame of reception (DESIGN 2.4), from which the
spherical elevation and azimuth the document asks for are computed. The
caller makes one and passes it as a ``CallerSource`` (DESIGN 2.1); the
engine feeds it the messages (``feed``).

Steps the document leaves open, each decided where it is taken: the
normalization of the Legendre functions (kernel ``ionosphere``), the t of
the sun-fixed longitude (``ionosphere``), the quality as a 1-sigma
(``_sigma``), a coefficient the field marks as not available
(``VtecLayer``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ssrcore.model.atmosphere import Absent
from ssrcore.model.providers import Kind, Look, Origin, SlantTec, Where, Zenith
from ssrcore.model.types import FeedOutcome, FeedReason, FeedResult, GnssTime, MaxAge, Missing, Sat, StreamId, Term
from ssrcore.toolkit import KernelError, geodesy
from ssrcore.toolkit import ionosphere as _ionosphere

# [spec: IGS SSR v1.00 Table 7 IDF041, 0.05 TECU per step; 0 "unknown"; 511
# (25.55) "indicator exceeds 25.54 TECU", a bound and not a value]
_QUALITY_STEP = 0.05
_QUALITY_ABOVE = 511
# [spec: IGS SSR v1.00 4.5.1, t is the "SSR epoch time of computation epoch modulo 86400"]
_DAY = 86400


@dataclass(frozen=True, slots=True)
class VtecLayer:
    """One IM201 layer in physical units: ``height`` above the 6370 km sphere (m),
    degree N and order M, cosine C(n,m) for m = 0..M, n = m..N and sine S(n,m) for
    m = 1..M, n = m..N in TECU, in the order of the message (IGS SSR 8.4.1).

    A coefficient the field marks as not available is None. [spec: IGS SSR
    v1.00 Table 7 IDF039, IDF040, range ±163.835 TECU, "163.84 indicates data
    out of range or not available" (IDF039) and "-163.84" (IDF040); design: of
    the two only −163.84 TECU, the int16 −32768, can be sent, so that value
    marks a coefficient of either field; a layer with one is not used, since
    no part of an expansion stands for the whole]"""

    height: float
    degree: int
    order: int
    cosine: tuple[float | None, ...]
    sine: tuple[float | None, ...]

    @property
    def complete(self) -> bool:
        """[design: a layer is evaluated only with every coefficient, see the class]"""
        return None not in self.cosine and None not in self.sine


@dataclass(frozen=True, slots=True)
class VtecInput:
    """An IM201 message in physical units: its stream (provider and solution, the
    adapter's ``IgsStreamId``), epoch (GPS time), update interval (s), IOD SSR,
    quality (TECU, None for unknown or beyond the field) and layers.
    [spec: IGS SSR v1.00 Tables 24 to 27, IDF003 to IDF009, IDF035 to IDF041]"""

    stream: StreamId
    epoch: GnssTime
    interval: float
    iod_ssr: int
    quality: float | None
    layers: tuple[VtecLayer, ...]


def quality(index: int) -> float | None:
    """[spec: IGS SSR v1.00 Table 7 IDF041, 0.05 TECU per step, 0 "unknown", 511
    "exceeds 25.54 TECU"; design: a bound below is not a value, as for SPARTN's
    last quality step (``providers._sigma``)] the VTEC quality in TECU, or None."""
    if index <= 0 or index >= _QUALITY_ABOVE:
        return None
    return index * _QUALITY_STEP


class IgsVtecIonosphere:
    """The IGS SSR VTEC of one stream as an ionosphere provider (A9).

    ``stream`` None takes the VTEC of the query's satellite stream (the
    binding of ``Where``), under the same IOD SSR; a named stream (an
    ``IgsStreamId``) takes that stream's VTEC whatever gives the satellite
    corrections. ``max_age`` bounds the age of a message, and an age within
    its update interval is always allowed.
    [spec: IGS SSR v1.00 4, "The validity interval of SSR parameters is at
    least the SSR Update Interval"; IDF007, a change of IOD SSR is a change of
    the generating configuration; design: rule A1 read for this format, a
    satellite stream's own ionosphere under its own issue; rule A10 for a named
    stream; the age bound of DESIGN 2.5]
    """

    kind = Kind.CORRECTION
    terms = frozenset({Term.IONOSPHERE})

    def __init__(self, name: str = "igs_vtec", stream: StreamId | None = None,
                 max_age: float = MaxAge().atmosphere) -> None:
        """[design: DESIGN 3.6, a provider the caller makes and names]"""
        self.name = name
        self.stream = stream
        self.max_age = float(max_age)
        self._latest: dict[StreamId, list[VtecInput]] = {}
        self._origins: dict[StreamId, Origin] = {}

    # ---- feeding -----------------------------------------------------------------

    def feed(self, vtec: VtecInput) -> FeedResult:
        """Keep one IM201 message: the newest per stream, and the one before it, so a
        query between two epochs takes the one in force.
        [spec: IGS SSR v1.00 4, parameters apply from their epoch; design: DESIGN
        3.3 S8 for an equal message (duplicate) and an older one (superseded)]"""
        if self.stream is not None and vtec.stream != self.stream:
            return FeedResult(FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT,
                              detail=f"{self.name} takes the VTEC of {self.stream}, not of {vtec.stream}")
        kept = self._latest.setdefault(vtec.stream, [])
        if kept and vtec.epoch < kept[-1].epoch:
            return FeedResult(FeedOutcome.SUPERSEDED, FeedReason.OLDER, detail="IM201 older than the newest")
        if kept and vtec.epoch == kept[-1].epoch:
            if kept[-1] == vtec:
                return FeedResult(FeedOutcome.DUPLICATE, detail="IM201")
            kept[-1] = vtec
            return FeedResult(FeedOutcome.STORED, FeedReason.CONFLICT, detail="IM201 replaced at its epoch")
        kept.append(vtec)
        del kept[:-2]
        return FeedResult(FeedOutcome.STORED, detail=f"IM201, {len(vtec.layers)} layer(s)")

    # ---- the interface -----------------------------------------------------------

    def troposphere(self, where: Where) -> Zenith | Absent:
        """[design: IM201 carries no troposphere; ``terms`` says the walk never asks]"""
        return Absent(Missing.NO_ATMOSPHERE)

    def ionosphere(self, where: Where, sat: Sat, look: Look) -> SlantTec | Absent:
        """The slant TEC of ``sat`` along ``look`` at ``where``: the sum over the layers
        of VTEC at the pierce point (negative values replaced by 0) times 1/sin(E + ψ).
        [spec: IGS SSR v1.00 4.5.1, every step through ``toolkit.ionosphere``; design:
        the t of the sun-fixed longitude, "SSR epoch time of computation epoch modulo
        86400", is read as the query time, the computation epoch being the epoch the
        slant TEC is computed for, in GPS seconds of the day as the update intervals
        are (IDF004)]"""
        vtec = self._vtec(where)
        if isinstance(vtec, Absent):
            return vtec
        try:
            look_s = self._spherical(where, look)
            day = where.t.second % _DAY + where.t.fraction
            total, factor = 0.0, 0.0
            for layer in vtec.layers:
                p = _ionosphere.pierce(look_s, layer.height)
                v = _ionosphere.vtec(layer.degree, layer.order, layer.cosine, layer.sine, p.latitude, p.longitude, day)
                total += max(v, 0.0) * p.slant_factor
                factor = max(factor, p.slant_factor)
        except KernelError:
            return Absent(Missing.NOT_CONVERGED)
        sigma = None if vtec.quality is None else vtec.quality * factor
        return SlantTec(sat, total, self._origin(vtec.stream), vtec.epoch, where.t - vtec.epoch, sigma)

    # ---- the parts ---------------------------------------------------------------

    def _vtec(self, where: Where) -> VtecInput | Absent:
        """[spec: IGS SSR v1.00 4, IDF007; design: see the class] the message in force, or why none."""
        b = where.binding
        stream = self.stream
        if stream is None:
            if b is None:
                return Absent(where.unbound or Missing.NO_CORRECTIONS)
            stream = b.stream
        kept = self._latest.get(stream)
        if not kept:
            return Absent(Missing.NO_ATMOSPHERE)
        vtec = next((v for v in reversed(kept) if v.epoch <= where.t), None)
        if vtec is None or where.t - vtec.epoch > max(self.max_age, vtec.interval):
            return Absent(Missing.STALE)
        if b is not None and b.stream == stream and b.issue.value != vtec.iod_ssr:
            return Absent(Missing.SIOU_MISMATCH)
        if not vtec.layers or not all(layer.complete for layer in vtec.layers):
            return Absent(Missing.NO_ATMOSPHERE)
        return vtec

    @staticmethod
    def _spherical(where: Where, look: Look) -> _ionosphere.SphericalLook:
        """The receiver and line of sight in the sphere of IGS SSR 4.5.1, from the geodetic
        receiver and the ellipsoidal elevation and azimuth of the record.
        [spec: TM-23 B.1.2 (geodetic to ECEF), B.3 (the line of sight from its local
        east, north and up components); IGS SSR v1.00 4.5.1 through
        ``ionosphere.spherical``; design: DESIGN 2.8, the receiver once per epoch in
        ``Where.cache``]"""
        key = ("igs_vtec.receiver", where.latitude, where.longitude, where.height)
        receiver = where.cache.get(key)
        if receiver is None:
            receiver = where.cache[key] = geodesy.geodetic_to_ecef((where.latitude, where.longitude, where.height))
        ce = math.cos(look.elevation)
        enu = (ce * math.sin(look.azimuth), ce * math.cos(look.azimuth), math.sin(look.elevation))
        unit = geodesy.enu_to_ecef(enu, where.latitude, where.longitude)
        return _ionosphere.spherical(receiver, unit)

    def _origin(self, stream: StreamId) -> Origin:
        """[design: DESIGN 3.6, the record's origin names the stream, which
        rule A11 compares with the code biases' stream]"""
        origin = self._origins.get(stream)
        if origin is None:
            origin = self._origins[stream] = Origin(self.name, self.kind, stream)
        return origin
