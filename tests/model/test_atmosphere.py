"""Atmosphere store and SPARTN GAD/HPAC adapter: rules A1 to A8 (DESIGN 3.3).

Every message is synthetic: built as a ``gnsswire.spartn`` table from chosen
physical values, encoded to a frame, parsed and decoded back, and converted
by ``adapters.spartn_atmosphere``. The issues of update (rules S3, S4, A2)
are those of the real stream store, ``model.stream.Streams``.
"""

from __future__ import annotations

import math

import pytest
from gnsswire import spartn

from ssrcore.adapters import spartn_atmosphere as adapter
from ssrcore.model.atmosphere import (
    Absent,
    AreaChoice,
    AtmosphereStore,
    IonosphereTec,
    TroposphereZenith,
)
from ssrcore.model.stream import AIOU_LIFETIME as AIOU_OUTAGE
from ssrcore.model.stream import Issue, Streams
from ssrcore.model.types import FeedOutcome, FeedReason, GnssTime, Missing, Sat, StreamId, System

# ---- the stream store's SIOU issue, as the selection would name it ------------------------


class SiouOf:
    """The SIOU ``Issue`` a selection would use: the value with the stream's current era."""

    def __init__(self, streams: Streams) -> None:
        self.streams = streams

    def token(self, stream, siou):
        s = self.streams.state(stream)
        return Issue(siou, 0 if s is None else s.siou.era)


# ---- encoding physical values (SPARTN 6.3, Table 6.2) ------------------------------------

# SF: (resolution, range minimum, offset of the notes)
FIELDS = {
    32: (0.1, -90.0, 0.0), 33: (0.1, -180.0, 0.0), 36: (0.1, 0.1, 0.0), 37: (0.1, 0.1, 0.0),
    43: (0.004, -0.508, 2.3),
    45: (0.004, -0.252, 0.252), 46: (0.001, -0.063, 0.0), 47: (0.0002, -0.051, 0.0),
    48: (0.004, -1.020, 0.252), 49: (0.001, -0.255, 0.0), 50: (0.0002, -0.2046, 0.0),
    52: (0.004, -0.124, 0.0), 53: (0.004, -0.508, 0.0),
    57: (0.04, -81.88, 0.0), 58: (0.008, -16.376, 0.0), 59: (0.002, -8.190, 0.0),
    60: (0.04, -327.64, 0.0), 61: (0.008, -65.528, 0.0), 62: (0.002, -32.766, 0.0),
    64: (0.04, -0.28, 0.0), 65: (0.04, -2.52, 0.0), 66: (0.04, -20.44, 0.0), 67: (0.04, -327.64, 0.0),
}
INVALID = {52: 0x3F, 53: 0xFF, 64: 0xF, 65: 0x7F, 66: 0x3FF, 67: 0x3FFF}


def enc(sf: int, value: float | None) -> int:
    """[spec: SPARTN 6.3] Vencoded = (Voriginal − Rngmin) / Res, after the note's offset."""
    if value is None:
        return INVALID[sf]
    res, lo, offset = FIELDS[sf]
    return round((value - offset - lo) / res)


STREAM = StreamId(5, 3)
T0 = GnssTime(2300, 100000.0)
SECONDS_2010 = 1_000_000_000


def _wire(message_type: int, subtype: int, body, stream: StreamId, tag: int) -> spartn.Message:
    """Encode to a frame and decode it back, as a receiver would."""
    d = spartn.Description(subtype=subtype, time_tag_type=1, time_tag=tag,
                           solution_id=stream.solution, processor_id=stream.processor,
                           encryption_id=None, sequence=None, auth_indicator=None, auth_length=None)
    raw = spartn.encode(spartn.Message(message_type, 0, 2, d, body))
    frame = spartn.parse(raw)
    assert isinstance(frame, spartn.Frame)
    return spartn.decode(frame)


def gad_area(area_id=1, lat=50.0, lon=10.0, nlat=3, nlon=4, dlat=1.0, dlon=1.5):
    return spartn.GadArea(area_id, enc(32, lat), enc(33, lon), nlat - 1, nlon - 1, enc(36, dlat), enc(37, dlon))


def gad_message(*areas, siou=7, aiou=1, stream=STREAM, tag=SECONDS_2010):
    areas = areas or (gad_area(),)
    return _wire(2, 0, spartn.Gad(siou, aiou, 0, len(areas) - 1, tuple(areas)), stream, tag)


def tropo(equation_type=2, coefficients=(0.1, 0.01, -0.02, 0.004), hydrostatic=2.3, quality=2,
          large=False, residuals=None, large_residuals=False):
    n = {0: 1, 1: 3, 2: 4}.get(equation_type, 4)
    fields = (48, 49, 49, 50) if large else (45, 46, 46, 47)
    raw = [enc(f, c) for f, c in zip(fields, coefficients)][:n] + [None] * (4 - n)
    block = (spartn.TropoLarge if large else spartn.TropoSmall)(*raw)
    sf = 53 if large_residuals else 52
    return spartn.Troposphere(
        equation_type, quality, enc(43, hydrostatic), int(large),
        None if large else block, block if large else None,
        None if residuals is None else int(large_residuals),
        None if residuals is None else tuple(enc(sf, r) for r in residuals))


def iono_sat(coefficients=(10.0, 0.4, -0.8, 0.02), equation_type=2, quality=5, large=False,
             residuals=None, residual_size=1):
    n = {0: 1, 1: 3, 2: 4}[equation_type]
    fields = (60, 61, 61, 62) if large else (57, 58, 58, 59)
    raw = [enc(f, c) for f, c in zip(fields, coefficients)][:n] + [None] * (4 - n)
    block = (spartn.IonoLarge if large else spartn.IonoSmall)(*raw)
    sf = (64, 65, 66, 67)[residual_size]
    return spartn.IonoSatellite(
        quality, int(large), None if large else block, block if large else None,
        None if residuals is None else residual_size,
        None if residuals is None else tuple(enc(sf, r) for r in residuals))


def mask_of(system: System, numbers, mask_size=0):
    width = adapter.SATELLITE_MASK_BITS[system][mask_size]
    return sum(1 << (width - n) for n in numbers)


def hpac_area(area_id=1, troposphere=None, sats=None, system=System.GPS, equation_type=2,
              grid_points=0, grid=False):
    """An HPAC area: ``sats`` maps satellite number to an ``iono_sat``."""
    iono = None
    if sats:
        numbers = sorted(sats)
        iono = spartn.Ionosphere(equation_type, 0, mask_of(system, numbers),
                                 tuple(sats[n] for n in numbers))
    blocks = 2 if grid else 1
    return spartn.HpacArea(area_id, grid_points, blocks if troposphere else 0, blocks if iono else 0,
                           troposphere, iono)


def hpac_message(*areas, system=System.GPS, siou=7, aiou=1, stream=STREAM, tag=SECONDS_2010):
    areas = areas or (hpac_area(troposphere=tropo(), sats={5: iono_sat()}),)
    return _wire(1, int(system), spartn.Hpac(siou, aiou, 0, len(areas) - 1, tuple(areas)), stream, tag)


def make_store(max_age=120.0):
    streams = Streams()
    return AtmosphereStore(streams, max_age), SiouOf(streams)


def feed_gad(store, t, **kw):
    return store.feed_gad(adapter.from_gad(gad_message(*kw.pop("areas", ()), **kw)), t)


def feed_hpac(store, t, *areas, **kw):
    return store.feed_hpac(adapter.from_hpac(hpac_message(*areas, **kw)), t)


def centroid(lat, lon, nlat, nlon, dlat, dlon):
    """[spec: SPARTN 8.10 figure, DESIGN Q8] centre of the node grid."""
    return lat - (nlat - 1) * dlat / 2, lon + (nlon - 1) * dlon / 2


# ---- adapter ---------------------------------------------------------------------------


def test_gad_adapter_scales_table_6_24():
    m = gad_message(gad_area(9, lat=50.0, lon=10.0, nlat=3, nlon=4, dlat=1.0, dlon=1.5),
                    gad_area(10, lat=-30.5, lon=-179.9, nlat=1, nlon=8, dlat=3.2, dlon=0.1),
                    siou=300, aiou=11, tag=1234)
    g = adapter.from_gad(m)
    assert (g.siou, g.aiou, g.stream, g.tag.value, g.tag.full) == (300, 11, STREAM, 1234, True)
    a, b = g.areas
    assert a.area_id == 9 and a.latitude_nodes == 3 and a.longitude_nodes == 4
    assert a.latitude == pytest.approx(50.0, abs=1e-12) and a.longitude == pytest.approx(10.0, abs=1e-12)
    assert a.latitude_spacing == pytest.approx(1.0) and a.longitude_spacing == pytest.approx(1.5)
    assert b.latitude == pytest.approx(-30.5) and b.longitude == pytest.approx(-179.9)
    assert (b.latitude_nodes, b.longitude_nodes) == (1, 8)
    assert b.latitude_spacing == pytest.approx(3.2) and b.longitude_spacing == pytest.approx(0.1)


def test_gad_adapter_range_ends_and_undefined_values():
    # SF032 raw 0 is −90°, SF033 raw 0 is −180°: an area at the south-west end is
    # still refused because its grid reaches below −90° only if it has rows.
    g = adapter.from_gad(gad_message(gad_area(1, lat=-89.0, lon=-180.0, nlat=1, nlon=1)))
    assert g.areas[0].longitude == -180.0
    with pytest.raises(ValueError):          # SF032 = 2047: 114.7°, outside ±90°
        adapter.from_gad(gad_message(spartn.GadArea(1, 2047, 0, 0, 0, 0, 0)))
    with pytest.raises(ValueError):          # three rows of 3.2° south of −85°
        adapter.from_gad(gad_message(gad_area(1, lat=-85.0, nlat=3, dlat=3.2)))
    with pytest.raises(TypeError):
        adapter.from_gad(hpac_message())
    with pytest.raises(TypeError):
        adapter.from_hpac(gad_message())


def test_hpac_adapter_troposphere_offsets_and_residuals():
    t = tropo(2, (0.1, 0.012, -0.02, 0.004), hydrostatic=2.1, quality=3,
              residuals=(0.012, None, -0.124, 0.0, 0.1, 0.12))
    h = adapter.from_hpac(hpac_message(hpac_area(4, troposphere=t, grid_points=6, grid=True),
                                       system=System.GALILEO, siou=12, aiou=3))
    assert (h.system, h.siou, h.aiou) == (System.GALILEO, 12, 3)
    (area,) = h.areas
    assert area.area_id == 4 and area.grid_points == 6 and area.ionosphere is None
    tb = area.troposphere
    assert tb.equation_type == 2 and tb.quality == 3
    assert tb.hydrostatic == pytest.approx(2.1, abs=1e-12)        # SF043 offset from 2.3 m
    assert tb.coefficients == pytest.approx((0.1, 0.012, -0.02, 0.004), abs=1e-12)  # T00 + 0.252 m
    assert tb.residuals[1] is None                                 # 0x3F invalid
    assert tb.residuals[0] == pytest.approx(0.012) and tb.residuals[2] == pytest.approx(-0.124)


def test_hpac_adapter_large_blocks_and_equation_types():
    t = tropo(1, (-0.7, 0.2, -0.1), large=True, residuals=(0.5, None), large_residuals=True)
    sats = {1: iono_sat((100.0,), equation_type=0, large=True, residuals=(300.0, None), residual_size=3),
            32: iono_sat((-5.0,), equation_type=0, residuals=(0.2, -0.28), residual_size=0)}
    h = adapter.from_hpac(hpac_message(hpac_area(2, troposphere=t, sats=sats, equation_type=0,
                                                 grid_points=2, grid=True)))
    tb = h.areas[0].troposphere
    assert tb.coefficients == pytest.approx((-0.7, 0.2, -0.1), abs=1e-12)
    assert tb.residuals[0] == pytest.approx(0.5) and tb.residuals[1] is None   # 0xFF invalid
    io = h.areas[0].ionosphere
    assert io.equation_type == 0
    assert [s.sat for s in io.satellites] == [Sat.gps(1), Sat.gps(32)]
    assert io.satellites[0].coefficients == pytest.approx((100.0,))
    assert io.satellites[0].residuals[0] == pytest.approx(300.0) and io.satellites[0].residuals[1] is None
    assert io.satellites[1].residuals == pytest.approx((0.2, -0.28))


def test_satellite_masks_of_each_system():
    # [spec: SF011, SF012, SF093 to SF095] position 0 is the leftmost bit, number = position + 1
    assert adapter.satellites(System.GPS, 0, 1 << 31) == (Sat.gps(1),)
    assert adapter.satellites(System.GLONASS, 3, 1) == (Sat.glonass(63),)
    assert adapter.satellites(System.QZSS, 0, 0b1000000001) == (Sat.qzss(1), Sat.qzss(10))
    assert adapter.satellites(System.BEIDOU, 1, mask_of(System.BEIDOU, (2, 46), 1)) == \
        (Sat.beidou(2), Sat.beidou(46))
    with pytest.raises(ValueError):
        adapter.satellites(System.GPS, 0, 1 << 32)


def test_rule_A8_undefined_equation_types_never_reach_the_store():
    # gnsswire refuses SF041 and SF054 of 3 to 7 and SF040 = 3, on encoding and decoding.
    with pytest.raises(ValueError):
        hpac_message(hpac_area(troposphere=tropo(equation_type=3)))
    bad_iono = spartn.Ionosphere(5, 0, mask_of(System.GPS, (1,)), (iono_sat(),))
    with pytest.raises(ValueError):
        hpac_message(spartn.HpacArea(1, 0, 0, 1, None, bad_iono))
    with pytest.raises(ValueError):
        hpac_message(spartn.HpacArea(1, 0, 3, 0, tropo(), None))
    # A body that bypassed the codec is refused by the adapter as well.
    good = hpac_message()
    forged = spartn.Message(good.message_type, good.eaf, good.crc_type, good.description,
                            spartn.Hpac(7, 1, 0, 0, (hpac_area(troposphere=spartn.Troposphere(
                                4, 0, 0, 0, spartn.TropoSmall(0, 0, 0, 0), None, None, None)),)))
    with pytest.raises(ValueError, match="A8"):
        adapter.from_hpac(forged)


# ---- the store: rules A1 to A7 ------------------------------------------------------------

LAT, LON = 49.0, 12.0     # inside the default area (50 N to 48 N, 10 E to 14.5 E)


def ready(store, siou, t=T0, **kw):
    feed_gad(store, t, **kw)
    kw.pop("areas", None)
    feed_hpac(store, t + 1.0, **kw)
    return siou.token(STREAM, kw.get("siou", 7))


def test_rule_A1_hpac_combines_with_gad_of_same_area_aiou_and_siou():
    store, siou = make_store()
    token = ready(store, siou)
    t = T0 + 5.0
    choice = store.area(STREAM, token, LAT, LON, t)
    assert isinstance(choice, AreaChoice) and choice.area_id == 1
    assert isinstance(store.troposphere(choice, t), TroposphereZenith)
    assert isinstance(store.ionosphere(choice, Sat.gps(5), t), IonosphereTec)
    # Another SIOU: nothing combines (S3); the reason names it.
    assert store.area(STREAM, siou.token(STREAM, 8), LAT, LON, t) == Absent(Missing.SIOU_MISMATCH)


def test_rule_A1_aiou_mismatch():
    store, siou = make_store()
    feed_gad(store, T0, aiou=1)
    feed_hpac(store, T0 + 1.0, aiou=2)
    assert store.area(STREAM, siou.token(STREAM, 7), LAT, LON, T0 + 2.0) == \
        Absent(Missing.AIOU_MISMATCH)


def test_rule_A1_siou_mismatch_between_gad_and_hpac():
    store, siou = make_store()
    feed_gad(store, T0, siou=7)
    feed_hpac(store, T0 + 1.0, siou=8)
    assert store.area(STREAM, siou.token(STREAM, 7), LAT, LON, T0 + 2.0) == \
        Absent(Missing.SIOU_MISMATCH)
    assert store.area(STREAM, siou.token(STREAM, 8), LAT, LON, T0 + 2.0) == \
        Absent(Missing.SIOU_MISMATCH)


def test_rule_A1_hpac_of_another_area_id():
    store, siou = make_store()
    feed_gad(store, T0)
    feed_hpac(store, T0 + 1.0, hpac_area(2, troposphere=tropo()))
    assert store.area(STREAM, siou.token(STREAM, 7), LAT, LON, T0 + 2.0) == \
        Absent(Missing.NO_ATMOSPHERE)
    assert store.area(StreamId(1, 1), siou.token(STREAM, 7), LAT, LON, T0) == \
        Absent(Missing.NO_ATMOSPHERE)


def test_rule_A1_aiou_change_switches_to_the_new_definition():
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0, aiou=1, areas=(gad_area(1, lat=50.0, lon=10.0),))
    feed_hpac(store, T0 + 1.0, aiou=1)
    # The provider moves area 1 one degree east under AIOU 2; the new GAD arrives first.
    feed_gad(store, T0 + 30.0, aiou=2, areas=(gad_area(1, lat=50.0, lon=11.0),))
    choice = store.area(STREAM, token, LAT, LON, T0 + 31.0)
    assert choice.definition.longitude == pytest.approx(10.0)       # old pair still combines
    feed_hpac(store, T0 + 32.0, aiou=2)
    choice = store.area(STREAM, token, LAT, LON, T0 + 33.0)
    assert choice.definition.longitude == pytest.approx(11.0)       # the newest data wins
    assert store.troposphere(choice, T0 + 33.0).time == T0 + 32.0
    # Outside the new definition but inside the old one, the old pair still serves.
    choice = store.area(STREAM, token, LAT, 10.5, T0 + 33.0)
    assert choice.definition.longitude == pytest.approx(10.0)


def test_rule_A2_outage_of_450_s_makes_every_aiou_new():
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0)
    feed_hpac(store, T0 + AIOU_OUTAGE - 1.0)         # 449 s: same AIOU, combines
    assert isinstance(store.area(STREAM, token, LAT, LON, T0 + AIOU_OUTAGE), AreaChoice)

    store, siou = make_store(max_age=400.0)
    feed_gad(store, T0)
    feed_hpac(store, T0 + AIOU_OUTAGE)               # 450 s: equal AIOU, but new
    # The definition of the old generation can no longer combine with anything and
    # is forgotten, so there is no area at all.
    assert store.area(STREAM, token, LAT, LON, T0 + AIOU_OUTAGE + 1.0) == \
        Absent(Missing.NO_ATMOSPHERE)
    feed_gad(store, T0 + AIOU_OUTAGE + 2.0)          # the definition repeated after the outage
    assert isinstance(store.area(STREAM, token, LAT, LON, T0 + AIOU_OUTAGE + 3.0), AreaChoice)


def test_rule_A2_no_use_after_an_outage_up_to_the_query():
    store, siou = make_store(max_age=440.0)
    token = ready(store, siou)
    last = T0 + 1.0
    assert isinstance(store.area(STREAM, token, LAT, LON, last + 440.0), AreaChoice)
    assert store.area(STREAM, token, LAT, LON, last + AIOU_OUTAGE) == Absent(Missing.STALE)


def test_rule_A2_message_older_than_the_outage_is_rejected():
    store, siou = make_store()
    ready(store, siou, t=T0 + 1000.0)
    r = feed_hpac(store, T0)
    assert (r.outcome, r.reason) == (FeedOutcome.REJECTED, FeedReason.STALE)


def test_rule_S4_through_the_stream_state_siou_gap():
    # The SIOU token is the stream's (S4): after 510 s without a message the same SIOU
    # is new, so data before and after do not combine.
    store, siou = make_store(max_age=400.0)
    feed_gad(store, T0)
    old = siou.token(STREAM, 7)
    feed_hpac(store, T0 + 520.0)
    new = siou.token(STREAM, 7)
    assert old != new
    t = T0 + 521.0
    assert store.area(STREAM, new, LAT, LON, t) == Absent(Missing.NO_ATMOSPHERE)


def test_rule_A3_areas_keyed_by_stream_and_id():
    store, siou = make_store()
    other = StreamId(5, 4)                       # same solution, another processor
    feed_gad(store, T0, stream=other, areas=(gad_area(1, lat=50.0, lon=10.0),))
    feed_hpac(store, T0 + 1.0, stream=other)
    feed_gad(store, T0, areas=(gad_area(1, lat=20.0, lon=100.0),))
    feed_hpac(store, T0 + 1.0)
    token = siou.token(STREAM, 7)
    assert store.area(STREAM, token, LAT, LON, T0 + 2.0) == Absent(Missing.OUTSIDE_AREA)
    choice = store.area(other, siou.token(other, 7), LAT, LON, T0 + 2.0)
    assert choice.stream == other and choice.definition.latitude == pytest.approx(50.0)


def test_rule_A4_inside_the_area_edges_included():
    store, siou = make_store()
    token = ready(store, siou)
    t = T0 + 2.0
    for lat, lon in ((50.0, 10.0), (48.0, 14.5), (49.5, 12.25)):
        assert isinstance(store.area(STREAM, token, lat, lon, t), AreaChoice)
    for lat, lon in ((50.05, 12.0), (47.9, 12.0), (49.0, 9.9), (49.0, 14.6)):
        assert store.area(STREAM, token, lat, lon, t) == Absent(Missing.OUTSIDE_AREA)


def test_rule_A4_not_beyond_85_degrees():
    store, siou = make_store()
    token = ready(store, siou, areas=(gad_area(1, lat=87.0, lon=0.0, nlat=2, nlon=2, dlat=3.0, dlon=3.0),))
    t = T0 + 2.0
    assert isinstance(store.area(STREAM, token, 85.0, 1.0, t), AreaChoice)
    assert store.area(STREAM, token, 86.0, 1.0, t) == Absent(Missing.LATITUDE_LIMIT)


def test_rule_A4_overlapping_areas_nearest_centroid():
    store, siou = make_store()
    a1 = gad_area(1, lat=50.0, lon=10.0, nlat=3, nlon=3, dlat=1.0, dlon=1.0)   # centroid 49 N, 11 E
    a2 = gad_area(2, lat=51.0, lon=11.0, nlat=3, nlon=3, dlat=1.0, dlon=1.0)   # centroid 50 N, 12 E
    feed_gad(store, T0, areas=(a1, a2))
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo(0, (0.10,))),
              hpac_area(2, troposphere=tropo(0, (0.20,)), sats={5: iono_sat((20.0,), 0)},
                        equation_type=0))
    token, t = siou.token(STREAM, 7), T0 + 2.0
    # Both areas contain 49.4 N 11.4 E; area 1's centroid is nearer.
    choice = store.area(STREAM, token, 49.4, 11.4, t)
    assert choice.area_id == 1
    assert store.troposphere(choice, t).polynomial == pytest.approx(0.10)
    # The same area serves the ionosphere, even though only area 2 has it.
    assert store.ionosphere(choice, Sat.gps(5), t) == Absent(Missing.NO_ATMOSPHERE)
    # Nearer area 2's centroid.
    choice = store.area(STREAM, token, 49.8, 11.9, t)
    assert choice.area_id == 2
    assert store.troposphere(choice, t).polynomial == pytest.approx(0.20)
    assert store.ionosphere(choice, Sat.gps(5), t).polynomial == pytest.approx(20.0)
    # Equidistant (49.5 N, 11.5 E on the local plane up to cos φ): the distance is
    # δφ² + (δλ cos φ)², equal here, so the lower area ID is chosen.
    assert store.area(STREAM, token, 49.5, 11.5, t).area_id == 1


def test_rule_A4_distance_uses_longitude_scaled_by_latitude():
    store, siou = make_store()
    a1 = gad_area(1, lat=62.0, lon=0.0, nlat=5, nlon=5, dlat=1.0, dlon=1.0)    # centroid 60 N, 2 E
    a2 = gad_area(2, lat=63.0, lon=1.0, nlat=3, nlon=3, dlat=1.5, dlon=1.0)    # centroid 61.5 N, 2 E
    feed_gad(store, T0, areas=(a1, a2))
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo()), hpac_area(2, troposphere=tropo()))
    # At 60.8 N, 3.0 E: area 1 is 0.8° north and 1° of longitude (0.49° at 60.8°) away,
    # area 2 0.7° south and 0.49°; area 2 is nearer, though not in raw degrees of latitude alone.
    d1 = 0.8 ** 2 + (1.0 * math.cos(math.radians(60.8))) ** 2
    d2 = 0.7 ** 2 + (1.0 * math.cos(math.radians(60.8))) ** 2
    assert d2 < d1
    assert store.area(STREAM, siou.token(STREAM, 7), 60.8, 3.0, T0 + 2.0).area_id == 2


def test_rule_A4_area_without_usable_data_is_passed_over():
    store, siou = make_store()
    a1 = gad_area(1, lat=50.0, lon=10.0, nlat=3, nlon=3, dlat=1.0, dlon=1.0)
    a2 = gad_area(2, lat=51.0, lon=11.0, nlat=3, nlon=3, dlat=1.0, dlon=1.0)
    feed_gad(store, T0, areas=(a1, a2))
    feed_hpac(store, T0 + 1.0, hpac_area(2, troposphere=tropo()))
    assert store.area(STREAM, siou.token(STREAM, 7), 49.4, 11.4, T0 + 2.0).area_id == 2


def test_rule_A5_latest_received_troposphere_of_an_epoch():
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0)
    t1 = T0 + 10.0
    feed_hpac(store, t1, hpac_area(1, troposphere=tropo(0, (0.10,))), system=System.GPS)
    feed_hpac(store, t1, hpac_area(1, troposphere=tropo(0, (0.12,))), system=System.GALILEO)
    choice = store.area(STREAM, token, LAT, LON, t1 + 1.0)
    assert store.troposphere(choice, t1 + 1.0).polynomial == pytest.approx(0.12)
    # An older epoch arriving later does not replace it.
    r = feed_hpac(store, t1 - 5.0, hpac_area(1, troposphere=tropo(0, (0.30,))), system=System.BEIDOU)
    assert (r.outcome, r.reason) == (FeedOutcome.SUPERSEDED, FeedReason.OLDER)
    assert store.troposphere(choice, t1 + 1.0).polynomial == pytest.approx(0.12)
    # An HPAC without a troposphere block leaves it as it is (SPARTN 8.4 item 2).
    feed_hpac(store, t1 + 2.0, hpac_area(1, sats={3: iono_sat()}), system=System.GPS)
    z = store.troposphere(choice, t1 + 3.0)
    assert z.polynomial == pytest.approx(0.12) and z.time == t1 and z.age == pytest.approx(3.0)


def test_a_query_before_the_newest_item_uses_the_item_of_its_time():
    # DESIGN 2.5: an item applies from its time tag until the next, as for OCB blocks;
    # a query between two HPAC epochs gets the earlier, not "stale" for the later.
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0)
    t1 = T0 + 10.0
    feed_hpac(store, t1, hpac_area(1, troposphere=tropo(0, (0.10,)), sats={5: iono_sat((10.0,), 0)},
                                   equation_type=0))
    feed_hpac(store, t1 + 5.0, hpac_area(1, troposphere=tropo(0, (0.14,)), sats={5: iono_sat((12.0,), 0)},
                                         equation_type=0))
    for t, value, tec in ((t1 + 2.0, 0.10, 10.0), (t1 + 5.0, 0.14, 12.0), (t1 + 7.0, 0.14, 12.0)):
        choice = store.area(STREAM, token, LAT, LON, t)
        z, i = store.troposphere(choice, t), store.ionosphere(choice, Sat.gps(5), t)
        assert z.polynomial == pytest.approx(value) and i.polynomial == pytest.approx(tec)
        assert z.age == i.age == pytest.approx(t - (t1 if t < t1 + 5.0 else t1 + 5.0))
    # Before the first item there is none.
    assert store.area(STREAM, token, LAT, LON, t1 - 1.0) == Absent(Missing.NO_ATMOSPHERE)


def test_rule_A6_ionosphere_per_satellite_of_its_constellation():
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0)
    r = feed_hpac(store, T0 + 1.0, hpac_area(1, sats={5: iono_sat((10.0,), 0), 7: iono_sat((12.0,), 0)},
                                             equation_type=0), system=System.GPS)
    assert r.outcome == FeedOutcome.STORED and r.sats == (Sat.gps(5), Sat.gps(7)) and r.areas == (1,)
    feed_hpac(store, T0 + 2.0, hpac_area(1, sats={5: iono_sat((30.0,), 0)}, system=System.GALILEO,
                                         equation_type=0), system=System.GALILEO)
    feed_hpac(store, T0 + 3.0, hpac_area(1, sats={5: iono_sat((11.0,), 0)}, equation_type=0),
              system=System.GPS)
    t = T0 + 4.0
    choice = store.area(STREAM, token, LAT, LON, t)
    assert store.ionosphere(choice, Sat.gps(5), t).polynomial == pytest.approx(11.0)
    assert store.ionosphere(choice, Sat.galileo(5), t).polynomial == pytest.approx(30.0)
    # GPS 7 was not in the latest GPS message; its own block is still the latest for it.
    i7 = store.ionosphere(choice, Sat.gps(7), t)
    assert i7.polynomial == pytest.approx(12.0) and i7.age == pytest.approx(3.0)
    assert store.ionosphere(choice, Sat.gps(9), t) == Absent(Missing.NO_ATMOSPHERE)
    assert store.ionosphere(choice, Sat.glonass(5), t) == Absent(Missing.NO_ATMOSPHERE)


def test_polynomials_about_the_area_centroid():
    store, siou = make_store()
    feed_gad(store, T0)
    tc = (0.1, 0.012, -0.02, 0.004)
    ic = (10.0, 0.4, -0.8, 0.02)
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo(2, tc, hydrostatic=2.2, quality=4),
                                         sats={5: iono_sat(ic, quality=6)}))
    t = T0 + 3.0
    lat, lon = 48.3, 13.9
    choice = store.area(STREAM, siou.token(STREAM, 7), lat, lon, t)
    phi0, lam0 = centroid(50.0, 10.0, 3, 4, 1.0, 1.5)
    dphi, dlam = lat - phi0, lon - lam0
    z = store.troposphere(choice, t)
    # [spec: SPARTN 8.12] T_P = T00 + T01 δφ + T10 δλ + T11 δφ δλ
    assert z.polynomial == pytest.approx(tc[0] + tc[1] * dphi + tc[2] * dlam + tc[3] * dphi * dlam, abs=1e-12)
    assert z.hydrostatic == pytest.approx(2.2, abs=1e-12) and z.quality == 4
    assert z.grid == 0.0 and z.grid_missing is None and z.wet == z.polynomial   # no grid block
    assert z.age == pytest.approx(2.0)
    i = store.ionosphere(choice, Sat.gps(5), t)
    # [spec: SPARTN 8.11] I_P = C00 + C01 δφ + C10 δλ + C11 δφ δλ
    assert i.polynomial == pytest.approx(ic[0] + ic[1] * dphi + ic[2] * dlam + ic[3] * dphi * dlam, abs=1e-9)
    assert i.grid == 0.0 and i.slant == i.polynomial and i.quality == 6


def bilinear(values, lat, lon, top=50.0, west=10.0, nlon=4, dlat=1.0, dlon=1.5):
    """Independent bilinear interpolation on the node grid (rows north to south, west to east)."""
    y = (top - lat) / dlat
    x = (lon - west) / dlon
    r, c = min(int(y), 1), min(int(x), nlon - 2)
    fy, fx = y - r, x - c
    v = lambda rr, cc: values[rr * nlon + cc]  # noqa: E731
    return ((1 - fy) * ((1 - fx) * v(r, c) + fx * v(r, c + 1))
            + fy * ((1 - fx) * v(r + 1, c) + fx * v(r + 1, c + 1)))


def test_rule_A7_grid_residuals_bilinear():
    store, siou = make_store()
    feed_gad(store, T0)
    tres = tuple(0.004 * k - 0.02 for k in range(12))
    ires = tuple(0.04 * k - 0.2 for k in range(12))
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo(residuals=tres),
                                         sats={5: iono_sat(residuals=ires)}, grid_points=12, grid=True))
    t = T0 + 2.0
    for lat, lon in ((49.3, 11.2), (48.0, 14.5), (50.0, 10.0), (48.7, 13.1)):
        choice = store.area(STREAM, siou.token(STREAM, 7), lat, lon, t)
        z = store.troposphere(choice, t)
        assert z.grid == pytest.approx(bilinear(tres, lat, lon), abs=1e-12)
        assert z.wet == pytest.approx(z.polynomial + z.grid)
        i = store.ionosphere(choice, Sat.gps(5), t)
        assert i.grid == pytest.approx(bilinear(ires, lat, lon), abs=1e-9)
        assert i.slant == pytest.approx(i.polynomial + i.grid)


def test_rule_A7_invalid_node_leaves_the_polynomial_alone():
    store, siou = make_store()
    feed_gad(store, T0)
    tres = [0.012] * 12
    tres[5] = None                        # node (1, 1)
    ires = [0.4] * 12
    ires[11] = None                       # node (2, 3), the south-east corner
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo(residuals=tres),
                                         sats={5: iono_sat(residuals=ires)}, grid_points=12, grid=True))
    t = T0 + 2.0
    token = siou.token(STREAM, 7)
    near = store.area(STREAM, token, 49.5, 11.0, t)           # cell with node 5
    z = store.troposphere(near, t)
    assert z.grid is None and z.grid_missing == Missing.INVALID_GRID_NODE
    assert z.wet == z.polynomial
    assert store.ionosphere(near, Sat.gps(5), t).grid == pytest.approx(0.4)
    corner = store.area(STREAM, token, 48.0, 14.5, t)         # node 11 with weight 1
    i = store.ionosphere(corner, Sat.gps(5), t)
    assert i.grid is None and i.grid_missing == Missing.INVALID_GRID_NODE
    assert i.slant == i.polynomial
    far = store.area(STREAM, token, 49.9, 14.0, t)            # cell of nodes 2, 3, 6, 7
    assert store.troposphere(far, t).grid == pytest.approx(0.012)


def test_rule_A7_residual_count_not_the_area_node_count():
    store, siou = make_store()
    feed_gad(store, T0)
    feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo(residuals=(0.0,) * 6),
                                         grid_points=6, grid=True))
    choice = store.area(STREAM, siou.token(STREAM, 7), LAT, LON, T0 + 2.0)
    z = store.troposphere(choice, T0 + 2.0)
    assert z.grid is None and z.grid_missing == Missing.GRID_MISMATCH


def test_ages_limit_hpac_data():
    store, siou = make_store(max_age=60.0)
    token = ready(store, siou)
    hpac_time = T0 + 1.0
    choice = store.area(STREAM, token, LAT, LON, hpac_time + 60.0)
    assert store.troposphere(choice, hpac_time + 60.0).age == pytest.approx(60.0)
    assert store.troposphere(choice, hpac_time + 60.5) == Absent(Missing.STALE)
    assert store.area(STREAM, token, LAT, LON, hpac_time + 61.0) == Absent(Missing.STALE)
    # Data apply from their time tag, not before it: before the first there is none.
    assert store.area(STREAM, token, LAT, LON, hpac_time - 0.5) == Absent(Missing.NO_ATMOSPHERE)
    with pytest.raises(ValueError):
        AtmosphereStore(siou.streams, 450.0)


def test_gad_without_age_limit_while_its_aiou_is_carried():
    # The GAD is sent once; HPAC with the same AIOU keeps it in use beyond 450 s.
    store, siou = make_store()
    token = siou.token(STREAM, 7)
    feed_gad(store, T0)
    for k in range(1, 12):
        feed_hpac(store, T0 + 60.0 * k)
    t = T0 + 661.0
    choice = store.area(STREAM, token, LAT, LON, t)
    assert isinstance(choice, AreaChoice)
    assert store.troposphere(choice, t).age == pytest.approx(1.0)


def test_feed_results():
    store, siou = make_store()
    r = feed_gad(store, T0, areas=(gad_area(1), gad_area(2, lat=30.0)))
    assert (r.outcome, r.areas) == (FeedOutcome.STORED, (1, 2))
    r = feed_gad(store, T0 - 1.0, areas=(gad_area(1),))
    assert (r.outcome, r.reason) == (FeedOutcome.SUPERSEDED, FeedReason.OLDER)
    r = feed_hpac(store, T0)
    assert (r.outcome, r.areas, r.sats) == (FeedOutcome.STORED, (1,), (Sat.gps(5),))


def test_replay_gives_the_same_answers():
    def run():
        store, siou = make_store()
        feed_gad(store, T0, areas=(gad_area(1), gad_area(2, lat=51.0, lon=11.0)))
        feed_hpac(store, T0 + 1.0, hpac_area(1, troposphere=tropo()), hpac_area(2, troposphere=tropo(),
                                                                                 sats={5: iono_sat()}))
        choice = store.area(STREAM, siou.token(STREAM, 7), 49.5, 12.0, T0 + 2.0)
        return choice, store.troposphere(choice, T0 + 2.0), store.ionosphere(choice, Sat.gps(5), T0 + 2.0)

    assert run() == run()


# ---- live-service facts (DESIGN 5.1) ------------------------------------------------------


def test_live_area_with_constant_terms_only():
    """A service may send areas whose polynomials have the constant term only
    (SF041 and SF054 equation type 0) and no grid: troposphere and ionosphere are
    then the same everywhere in the area (SPARTN 8.11, 8.12)."""
    store, siou = make_store()
    feed_gad(store, T0)
    area = hpac_area(troposphere=tropo(equation_type=0, coefficients=(0.148,)),
                     sats={5: iono_sat((12.4,), equation_type=0), 9: iono_sat((-3.2,), equation_type=0)},
                     equation_type=0)
    assert feed_hpac(store, T0 + 1.0, area).outcome == FeedOutcome.STORED
    token = siou.token(STREAM, 7)
    t = T0 + 2.0
    for lat, lon in ((50.0, 10.0), (49.0, 12.0), (48.0, 14.5), (48.3, 10.2)):   # corners and inside
        choice = store.area(STREAM, token, lat, lon, t)
        assert isinstance(choice, AreaChoice)
        tz = store.troposphere(choice, t)
        assert tz.polynomial == pytest.approx(0.148, abs=1e-12) and tz.grid == 0.0 and tz.wet == tz.polynomial
        assert tz.hydrostatic == pytest.approx(2.3, abs=1e-12)
        assert store.ionosphere(choice, Sat.gps(5), t).slant == pytest.approx(12.4, abs=1e-12)
        assert store.ionosphere(choice, Sat.gps(9), t).slant == pytest.approx(-3.2, abs=1e-12)
