"""Atmosphere providers and their walk (DESIGN 3.6, rules A9, A10).

The SPARTN provider is checked against the atmosphere store it wraps: the
same area, the same values and the same reasons, and the slant delays the
engine computes from them bit for bit the same as from the store directly
(the arithmetic of ``Engine._troposphere`` and ``Engine._ionosphere`` is
repeated here, so that going through a provider changes no number).
The walk is checked on lists in several orders, with the UNB3m model and an
outside provider written here. Inputs are typed and synthetic.
"""

from __future__ import annotations

import math

import pytest

from ssrcore.model.atmosphere import Absent, AtmosphereStore
from ssrcore.model.providers import (
    AtmosphereProvider,
    Binding,
    Kind,
    Look,
    Origin,
    SlantTec,
    SpartnAtmosphere,
    Unb3mTroposphere,
    Where,
    Zenith,
    build,
    check,
    default_providers,
    walk_ionosphere,
    walk_troposphere,
)
from ssrcore.model.stream import Issue, Streams
from ssrcore.model.types import (
    AreaDefinition,
    CallerSource,
    Config,
    ModelSource,
    SpartnSource,
    GadInput,
    GnssTime,
    HpacArea,
    HpacInput,
    IonosphereBlock,
    IonosphereSatellite,
    Missing,
    Sat,
    SpartnTag,
    StreamId,
    System,
    Term,
    TroposphereBlock,
)
from ssrcore.toolkit import ionosphere, troposphere, unb3m

A = StreamId(1, 1)      # the satellite stream (orbit, clock, biases) in these tests
B = StreamId(2, 7)      # another provider's stream, atmosphere only
T = GnssTime(2300, 100000.0)
G5, G9 = Sat.gps(5), Sat.gps(9)
AREA = AreaDefinition(1, 50.0, 10.0, 3, 4, 1.0, 1.5)       # 48 to 50 N, 10 to 14.5 E
INSIDE = (49.3, 11.7, 350.0)                                # degrees, degrees, metres
OUTSIDE = (45.0, 11.7, 350.0)
LOOK = Look(math.radians(35.0), math.radians(120.0))
DOY = 150.25
F_L1, F_L2 = 1575.42e6, 1227.60e6


def gad(stream, siou=7, aiou=1, area=AREA):
    return GadInput(stream, SpartnTag(0, True), siou, aiou, (area,))


def hpac(stream, siou=7, aiou=1, hydrostatic=2.25, t00=0.1, c00=10.0, grid=True):
    residuals = tuple(0.001 * k for k in range(12)) if grid else None
    tropo = TroposphereBlock(2, 2, hydrostatic, (t00, 0.01, -0.02, 0.004), residuals)
    iono = IonosphereBlock(2, (IonosphereSatellite(G5, 5, (c00, 0.4, -0.8, 0.02)),
                               IonosphereSatellite(G9, 4, (c00 + 3.0, -0.1, 0.2, 0.0))))
    return HpacInput(System.GPS, stream, SpartnTag(0, True), siou, aiou,
                     (HpacArea(1, 12 if grid else 0, tropo, iono),))


def store_with(*streams_siou):
    streams = Streams()
    store = AtmosphereStore(streams, 120.0)
    for stream, siou in streams_siou:
        store.feed_gad(gad(stream, siou), T)
        store.feed_hpac(hpac(stream, siou), T)
    return store


def binding(store, stream=A, siou=7):
    return Binding(stream, Issue(siou, store.streams.state(stream).siou.era))


def where(position=INSIDE, bound=None, unbound=None, dt=10.0, cache=None):
    lat, lon, h = position
    return Where(T + dt, math.radians(lat), math.radians(lon), lat, lon, h, DOY, bound, unbound,
                 {} if cache is None else cache)


# ---- the SPARTN provider returns what the store returns -------------------------------------


def test_spartn_provider_returns_exactly_what_the_store_returns():
    store = store_with((A, 7))
    w = where(bound=binding(store))
    p = SpartnAtmosphere(store)
    area = store.area(A, w.binding.issue, w.latitude_deg, w.longitude_deg, w.t)
    assert p.area(w) == area
    z_store = store.troposphere(area, w.t)
    z = p.troposphere(w)
    assert z.spartn == z_store
    assert z.hydrostatic == z_store.hydrostatic * troposphere.hydrostatic_height_factor(w.latitude, w.height, DOY)
    assert z.wet == z_store.wet and (z.time, z.age) == (z_store.time, z_store.age)
    assert z.origin == Origin("spartn", Kind.CORRECTION, A)
    assert z.sigma == 0.020                                      # SF042 index 2: <= 0.020 m
    raw = SpartnAtmosphere(store, height_reduction=False).troposphere(where(bound=binding(store)))
    assert raw.hydrostatic == z_store.hydrostatic
    for sat in (G5, G9):
        tec_store = store.ionosphere(area, sat, w.t)
        tec = p.ionosphere(w, sat, LOOK)
        assert tec.spartn == tec_store and tec.tec == tec_store.slant and tec.sat == sat
        assert tec.origin == Origin("spartn", Kind.CORRECTION, A)
    assert p.ionosphere(w, G5, LOOK).sigma == 0.28               # SF055 index 5: <= 0.28 TECU


def test_slant_delays_from_the_provider_are_bit_identical_to_the_engine_today():
    # The arithmetic of Engine._troposphere and Engine._ionosphere, from the store
    # directly and from the provider.
    store = store_with((A, 7))
    lat, lon, h = INSIDE
    elevation = math.radians(23.0)
    w = where(bound=binding(store))
    area = store.area(A, w.binding.issue, lat, lon, w.t)
    z0 = store.troposphere(area, w.t)
    m = troposphere.niell_mapping(w.latitude, h, DOY, elevation)
    today = troposphere.slant(z0.hydrostatic * troposphere.hydrostatic_height_factor(w.latitude, h, DOY),
                              z0.polynomial, z0.grid or 0.0, m.hydrostatic, m.wet)
    z = SpartnAtmosphere(store).troposphere(w)
    assert troposphere.slant(z.hydrostatic, z.spartn.polynomial, z.spartn.grid or 0.0, m.hydrostatic, m.wet) == today
    assert troposphere.slant(z.hydrostatic, z.wet, 0.0, m.hydrostatic, m.wet) == today
    tec0 = store.ionosphere(area, G5, w.t)
    tec = SpartnAtmosphere(store).ionosphere(w, G5, LOOK)
    for f in (F_L1, F_L2):
        assert ionosphere.delay(tec.spartn.polynomial, f, tec.spartn.grid or 0.0) == \
            ionosphere.delay(tec0.polynomial, f, tec0.grid or 0.0)


def test_spartn_provider_gives_the_store_reasons():
    store = store_with((A, 7))
    p = SpartnAtmosphere(store)
    outside = where(OUTSIDE, bound=binding(store))
    assert p.troposphere(outside) == Absent(Missing.OUTSIDE_AREA)
    assert p.ionosphere(outside, G5, LOOK) == Absent(Missing.OUTSIDE_AREA)
    wrong_issue = Binding(A, Issue(8, 0))
    assert p.troposphere(where(bound=wrong_issue)) == Absent(Missing.SIOU_MISMATCH)
    old = where(bound=binding(store), dt=121.0)                  # max_age.atmosphere 120 s
    assert store.area(A, old.binding.issue, *INSIDE[:2], old.t) == Absent(Missing.STALE)
    assert p.troposphere(old) == p.ionosphere(old, G5, LOOK) == Absent(Missing.STALE)
    assert p.ionosphere(where(bound=binding(store)), Sat.gps(17), LOOK) == Absent(Missing.NO_ATMOSPHERE)
    # Without satellite corrections the selection's reason, as the engine gives it today.
    assert p.troposphere(where(unbound=Missing.STALE)) == Absent(Missing.STALE)
    assert p.ionosphere(where(unbound=Missing.DNU), G5, LOOK) == Absent(Missing.DNU)


def test_troposphere_and_ionosphere_share_one_area_per_epoch():
    # Rule A4: one area serves both terms of a query; DESIGN 2.8: chosen once per epoch.
    class Counting(AtmosphereStore):
        calls = 0

        def area(self, *args):
            Counting.calls += 1
            return super().area(*args)

    store = Counting(Streams(), 120.0)
    store.feed_gad(gad(A), T)
    store.feed_hpac(hpac(A), T)
    cache: dict = {}
    p = SpartnAtmosphere(store)
    for sat in (G5, G9):
        w = where(bound=binding(store), cache=cache)
        p.troposphere(w)
        p.ionosphere(w, sat, LOOK)
    assert Counting.calls == 1


# ---- another stream's atmosphere (rule A10) -------------------------------------------------


def test_atmosphere_of_another_stream_under_its_own_siou():
    store = store_with((B, 3))                                   # B sends GAD and HPAC only
    w = where(bound=Binding(A, Issue(7, 0)))                     # A's corrections, A's SIOU
    assert SpartnAtmosphere(store).troposphere(w) == Absent(Missing.NO_ATMOSPHERE)   # A has none
    p = SpartnAtmosphere(store, B)
    assert p.name == "spartn:2/7"
    z = p.troposphere(w)
    area = store.area(B, Issue(3, 0), *INSIDE[:2], w.t)
    assert z.spartn == store.troposphere(area, w.t) and z.origin == Origin("spartn:2/7", Kind.CORRECTION, B)
    assert p.ionosphere(w, G5, LOOK).spartn == store.ionosphere(area, G5, w.t)
    # B's SIOU changes: the provider follows the stream's current issue.
    store.feed_gad(gad(B, 4), T + 5.0)
    store.feed_hpac(hpac(B, 4, t00=0.2), T + 5.0)
    z4 = p.troposphere(where(bound=Binding(A, Issue(7, 0))))
    assert z4.spartn.polynomial != z.spartn.polynomial
    assert z4.spartn == store.troposphere(store.area(B, Issue(4, 0), *INSIDE[:2], T + 10.0), T + 10.0)
    # A named stream that is also the satellite stream keeps the binding's issue (A1).
    both = store_with((A, 7))
    named, bound_p = SpartnAtmosphere(both, A, name="a"), SpartnAtmosphere(both)
    w = where(bound=binding(both))
    assert named.troposphere(w).spartn == bound_p.troposphere(w).spartn
    assert named.troposphere(where(bound=Binding(A, Issue(8, 0)))) == Absent(Missing.SIOU_MISMATCH)
    assert SpartnAtmosphere(both, StreamId(9, 9)).troposphere(w) == Absent(Missing.NO_ATMOSPHERE)


# ---- the UNB3m provider ---------------------------------------------------------------------


def test_unb3m_provider_is_the_toolkit_model_labelled_a_model():
    p = Unb3mTroposphere()
    w = where(OUTSIDE)
    z = p.troposphere(w)
    expected = unb3m.zenith_delays(w.latitude, w.height, DOY)
    assert (z.hydrostatic, z.wet) == expected
    assert z.origin == Origin("unb3m", Kind.MODEL) and z.time is None and z.age is None
    assert z.sigma == unb3m.TOTAL_SIGMA and z.spartn is None
    assert p.terms == frozenset({Term.TROPOSPHERE})
    geoid = Unb3mTroposphere(lambda lat, lon: 47.0)
    assert geoid.troposphere(w).hydrostatic == unb3m.zenith_delays(w.latitude, w.height - 47.0, DOY).hydrostatic
    assert p.troposphere(where((45.0, 11.7, 60_000.0))) == Absent(Missing.NOT_CONVERGED)


# ---- an outside provider, written here ------------------------------------------------------


class ThinShellMap:
    """A caller's ionosphere: a vertical TEC on a single layer at 450 km, mapped by
    1/cos of the zenith angle at the pierce point. Ready only from ``ready``."""

    kind = Kind.CORRECTION
    terms = frozenset({Term.IONOSPHERE})

    def __init__(self, vtec: float, ready: GnssTime, name: str = "shell") -> None:
        self.vtec, self.ready, self.name = vtec, ready, name

    def troposphere(self, where):
        return Absent(Missing.NO_ATMOSPHERE)

    def ionosphere(self, where, sat, look):
        if where.t < self.ready:
            return Absent(Missing.STALE)                 # not yet
        sin_z = 6371.0 / (6371.0 + 450.0) * math.cos(look.elevation)
        return SlantTec(sat, self.vtec / math.sqrt(1.0 - sin_z * sin_z), Origin(self.name, self.kind),
                        self.ready, where.t - self.ready, 3.0)


class FixedZenith:
    """A caller's troposphere: zenith delays measured at the site (a meteorological sensor)."""

    kind = Kind.CORRECTION
    terms = frozenset({Term.TROPOSPHERE})
    name = "sensor"

    def troposphere(self, where):
        return Zenith(2.21, 0.13, Origin(self.name, self.kind), where.t, 0.0, 0.01)

    def ionosphere(self, where, sat, look):
        return Absent(Missing.NO_ATMOSPHERE)


def test_outside_providers_have_the_interface():
    assert isinstance(ThinShellMap(10.0, T), AtmosphereProvider)
    assert isinstance(FixedZenith(), AtmosphereProvider)
    assert isinstance(Unb3mTroposphere(), AtmosphereProvider)
    assert isinstance(SpartnAtmosphere(store_with()), AtmosphereProvider)
    assert not isinstance(object(), AtmosphereProvider)


# ---- the walk (rule A9) ---------------------------------------------------------------------


def test_default_list_is_the_bound_spartn_store_alone():
    store = store_with((A, 7))
    (p,) = check(default_providers(store, True))
    assert isinstance(p, SpartnAtmosphere) and p.stream is None and p.height_reduction
    w = where(bound=binding(store))
    assert walk_troposphere((p,), w) == p.troposphere(w)
    outside = walk_troposphere((p,), where(OUTSIDE, bound=binding(store)))
    assert outside == Absent(Missing.OUTSIDE_AREA)


def test_the_first_answer_wins_and_a_term_never_had_is_not_asked():
    store = store_with((A, 7))
    shell = ThinShellMap(12.0, T)
    providers = check((SpartnAtmosphere(store), shell, Unb3mTroposphere()))
    # Inside the area SPARTN answers both terms.
    w = where(bound=binding(store))
    tw, iw = walk_troposphere(providers, w), walk_ionosphere(providers, w, G5, LOOK)
    assert tw.origin.provider == "spartn" and iw.origin.provider == "spartn"
    # Outside it the troposphere falls to the model; the shell map answers the
    # ionosphere; UNB3m, which never has ionosphere, is not asked for it.
    w = where(OUTSIDE, bound=binding(store))
    tw, iw = walk_troposphere(providers, w), walk_ionosphere(providers, w, G5, LOOK)
    assert tw.origin == Origin("unb3m", Kind.MODEL)
    assert iw.origin == Origin("shell", Kind.CORRECTION)
    assert iw.tec == pytest.approx(12.0 / math.sqrt(1 - (6371 / 6821 * math.cos(LOOK.elevation)) ** 2))
    # Before the shell map is ready ("not yet") nothing answers the ionosphere: the
    # reason is the first provider's, not the shell map's (stale).
    early = where(OUTSIDE, bound=binding(store), dt=-5.0)
    assert shell.ionosphere(early, G5, LOOK) == Absent(Missing.STALE)
    assert walk_ionosphere(providers, early, G5, LOOK) == Absent(Missing.OUTSIDE_AREA)


def test_order_is_the_callers():
    store = store_with((A, 7))
    w = where(bound=binding(store))
    model_first = check((Unb3mTroposphere(), FixedZenith(), SpartnAtmosphere(store)))
    assert walk_troposphere(model_first, w).origin.provider == "unb3m"
    assert walk_ionosphere(model_first, w, G5, LOOK).origin.provider == "spartn"
    sensor_first = check((FixedZenith(), SpartnAtmosphere(store), Unb3mTroposphere()))
    assert walk_troposphere(sensor_first, w).hydrostatic == 2.21
    # Two SPARTN streams: one provider's satellite corrections (A, no atmosphere here)
    # and another's regional atmosphere (B).
    regional = store_with((B, 3))
    two = check((SpartnAtmosphere(regional), SpartnAtmosphere(regional, B), Unb3mTroposphere()))
    tw = walk_troposphere(two, where(bound=Binding(A, Issue(7, 0))))
    assert tw.origin == Origin("spartn:2/7", Kind.CORRECTION, B)


def test_an_empty_list_and_an_unbound_query():
    assert walk_troposphere((), where()) == Absent(Missing.NO_ATMOSPHERE)
    assert walk_ionosphere((), where(), G5, LOOK) == Absent(Missing.NO_ATMOSPHERE)
    only_model = check((Unb3mTroposphere(),))
    assert walk_ionosphere(only_model, where(), G5, LOOK) == walk_ionosphere((), where(), G5, LOOK)
    # Without satellite corrections the model still gives a troposphere.
    store = store_with((A, 7))
    providers = check((SpartnAtmosphere(store), Unb3mTroposphere()))
    assert walk_troposphere(providers, where(unbound=Missing.NO_CORRECTIONS)).origin.kind == Kind.MODEL


def test_check_refuses_what_cannot_be_walked():
    store = store_with()
    with pytest.raises(ValueError, match="used twice"):
        check((SpartnAtmosphere(store), SpartnAtmosphere(store)))
    with pytest.raises(TypeError):
        check((object(),))

    class Bad(FixedZenith):
        kind = "model"
    with pytest.raises(TypeError):
        check((Bad(),))

    class Wide(FixedZenith):
        terms = frozenset({Term.TROPOSPHERE, Term.CODE_BIAS})
    with pytest.raises(ValueError, match="not atmosphere"):
        check((Wide(),))


# ---- the list built from the configuration (P4) ---------------------------------------------


def test_build_gives_one_provider_per_source_in_order():
    store = store_with((A, 7), (B, 3))
    (only,) = build(Config(), store)
    (default,) = default_providers(store, True)
    assert type(only) is type(default) and (only.stream, only.height_reduction, only.name) == \
        (default.stream, default.height_reduction, default.name)
    assert not build(Config(troposphere_height_reduction=False), store)[0].height_reduction
    sensor = FixedZenith()
    config = Config(satellite_stream=A, atmosphere=(SpartnSource(), SpartnSource(B), ModelSource(),
                                                    CallerSource("sensor")))
    providers = build(config, store, {"sensor": sensor})
    assert [p.name for p in providers] == ["spartn", "spartn:2/7", "unb3m", "sensor"]
    assert providers[1].stream == B and providers[1].store is store
    assert isinstance(providers[2], Unb3mTroposphere) and providers[3] is sensor
    assert build(Config(atmosphere=()), store) == ()


def test_build_refuses_what_does_not_match():
    store = store_with()
    config = Config(atmosphere=(SpartnSource(), CallerSource("sensor")))
    with pytest.raises(ValueError, match="do not match"):
        build(config, store)                                    # a CallerSource without its object
    with pytest.raises(ValueError, match="do not match"):
        build(Config(), store, dict(sensor=FixedZenith()))      # an object without its CallerSource
    with pytest.raises(ValueError, match="named"):
        build(Config(atmosphere=(CallerSource("other"),)), store, dict(other=FixedZenith()))
    with pytest.raises(ValueError, match="used twice"):
        build(Config(atmosphere=(SpartnSource(), SpartnSource())), store)
    with pytest.raises(TypeError):
        build(Config(atmosphere=(CallerSource("x"),)), store, dict(x=object()))
