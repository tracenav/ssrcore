"""The engine end to end (DESIGN 2, 3.4).

Every input is synthetic and goes through the wire: RTCM 1019 and 1020 and
SPARTN GAD, HPAC and OCB are encoded with gnsswire from chosen values, framed
by ``adapters.framing.Framer`` and fed with receipt times. The broadcast
elements are the public 2024-01-01 records of tests/test_joins.py (G05 and
R01), received near GPS week 2295, 87 400 s, by a receiver that sees both at
about 70° elevation.

Each record term is compared with the toolkit's value for the same inputs,
computed here from the decoded messages; every ``Missing`` reason is produced
by a scenario below (``test_every_missing_reason_is_produced``); and a replay
of the same bytes gives the same records.
"""

from __future__ import annotations

import dataclasses
import math
import time
from types import MappingProxyType

import pytest
from gnsswire import rtcm3, spartn

from ssrcore.adapters import spartn_atmosphere, spartn_ocb
from ssrcore.adapters.framing import Framer
from ssrcore.engine import Engine
from ssrcore.model.atmosphere import Absent
from ssrcore.model.ephemeris import broadcast_state, glonass_orbit
from ssrcore.model.providers import Kind, Origin, Zenith
from ssrcore.model.selection import choose, transmit
from ssrcore.model.types import (
    SPEED_OF_LIGHT,
    CallerSource,
    Config,
    ModelSource,
    SpartnSource,
    FeedOutcome,
    FeedReason,
    GnssTime,
    KeplerEphemeris,
    LeapSeconds,
    MaxAge,
    Missing,
    OcbHeaderInput,
    Sat,
    Signal,
    StreamId,
    System,
    Term,
)
from ssrcore.record import IonosphereTerm, MissingTermError, TroposphereTerm
from ssrcore.toolkit import area, astro, geodesy, ionosphere, orbits, relativity, ssr, tides, timescale, troposphere
from ssrcore.toolkit import unb3m
from ssrcore.toolkit import windup as wu

PI_ICD = 3.1415926535898     # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2
LEAP = 18                    # GPS − UTC on 2024-01-01
T = GnssTime(2295, 87400.0)  # the epoch of the queries
TS = T - 2.0                 # the SPARTN messages' time
LAT, LON, HEIGHT = -51.4, -162.7, 100.0
RECEIVER = geodesy.geodetic_to_ecef((math.radians(LAT), math.radians(LON), HEIGHT))
STREAM = StreamId(1, 0)
SIOU, AIOU = 3, 1
G05, R01 = Sat.gps(5), Sat.glonass(1)
L1C, L2W, L5Q = Signal.GPS_L1C, Signal.GPS_L2W, Signal.GPS_L5Q
R1C = Signal.GLONASS_L1C

# ---- RTCM ephemerides (resolutions of RTCM Tables 3.5-21, 3.5-22) ------------------------

G05_ELEMENTS = dict(
    toc=93600.0, af0=-1.562642864883e-04, af1=-1.477928890381e-12, af2=0.0,
    toe=93600.0, sqrt_a=5153.783403397, e=0.005790898343548, m0=-2.494229256397,
    dn=4.453756945425e-09, omega=1.235180510432, omega0=-0.6299845609285,
    omega_dot=-8.155696860368e-09, i0=0.9693533956084, idot=-3.928735076157e-11,
    cuc=-1.747161149979e-06, cus=3.31737101078e-06, crc=319.34375, crs=-34.125,
    cic=1.098960638046e-07, cis=4.097819328308e-08)
R01_STATE = dict(  # tb index 13, GPS 2295 87 318 s
    position=(-11175834.47266, -9508662.597656, -20866544.92188),
    velocity=(2857.489585876, -39.36195373535, -1511.75403595),
    lunisolar=(-2.793967723846e-06, -1.862645149231e-06, 0.0),
    tau_n=-7.346831262112e-05, gamma_n=9.094947017729e-13)


def q(value, exponent, semicircle=False):
    """The broadcast integer of a value at resolution 2^exponent (angles from rad)."""
    return round((value / PI_ICD if semicircle else value) / 2.0**exponent)


def rtcm(cls, **fields) -> bytes:
    """[spec: RTCM 4] a frame of ``cls`` with every field 0 but those given."""
    values = {f.name: 0 for f in dataclasses.fields(cls) if f.name != "extension"}
    values.update(number=int(cls.__name__[3:]), **fields)
    return rtcm3.encode(cls(**values))


def frame_1019(iode=33) -> bytes:
    r = G05_ELEMENTS
    return rtcm(rtcm3.Msg1019, satellite=5, week=2295 % 1024, iode=iode, iodc=iode,
                toe=int(r["toe"] // 16), toc=int(r["toc"] // 16), af0=q(r["af0"], -31), af1=q(r["af1"], -43),
                sqrt_a=q(r["sqrt_a"], -19), e=q(r["e"], -33), m0=q(r["m0"], -31, True),
                delta_n=q(r["dn"], -43, True), omega=q(r["omega"], -31, True),
                omega0=q(r["omega0"], -31, True), omega_dot=q(r["omega_dot"], -43, True),
                i0=q(r["i0"], -31, True), idot=q(r["idot"], -43, True),
                cuc=q(r["cuc"], -29), cus=q(r["cus"], -29), cic=q(r["cic"], -29), cis=q(r["cis"], -29),
                crc=q(r["crc"], -5), crs=q(r["crs"], -5))


def frame_1020() -> bytes:
    r, km = R01_STATE, 1000.0
    return rtcm(rtcm3.Msg1020, satellite=1, channel=1 + 7, tb=13,
                x=q(r["position"][0] / km, -11), y=q(r["position"][1] / km, -11),
                z=q(r["position"][2] / km, -11),
                x_dot=q(r["velocity"][0] / km, -20), y_dot=q(r["velocity"][1] / km, -20),
                z_dot=q(r["velocity"][2] / km, -20),
                x_ddot=q(r["lunisolar"][0] / km, -30), y_ddot=q(r["lunisolar"][1] / km, -30),
                z_ddot=q(r["lunisolar"][2] / km, -30),
                tau=q(r["tau_n"], -30), gamma=q(r["gamma_n"], -40))


# ---- SPARTN (SPARTN v2.0.3 6.3, Table 6.2) -----------------------------------------------

SAT_BITS = {System.GPS: 32, System.GLONASS: 24}      # SF011, SF012 size 0
BIAS_BITS = {System.GPS: 6, System.GLONASS: 5}       # SF025, SF026 size 0
# SF: (resolution, range minimum, offset of the notes)
FIELDS = {
    20: (0.002, -16.382, 0.0), 29: (0.02, -20.46, 0.0),
    32: (0.1, -90.0, 0.0), 33: (0.1, -180.0, 0.0), 36: (0.1, 0.1, 0.0), 37: (0.1, 0.1, 0.0),
    43: (0.004, -0.508, 2.3), 45: (0.004, -0.252, 0.252), 46: (0.001, -0.063, 0.0), 47: (0.0002, -0.051, 0.0),
    52: (0.004, -0.124, 0.0), 57: (0.04, -81.88, 0.0), 58: (0.008, -16.376, 0.0), 59: (0.002, -8.190, 0.0),
    65: (0.04, -2.52, 0.0),
}
INVALID = {52: 0x3F, 65: 0x7F}


def enc(sf: int, value: float | None) -> int:
    """[spec: SPARTN 6.3] Vencoded = (Voriginal − Rngmin) / Res, after the note's offset."""
    if value is None:
        return INVALID[sf]
    res, lo, offset = FIELDS[sf]
    return round((value - offset - lo) / res)


def description(system: System, t: GnssTime, subtype: int, stream: StreamId = STREAM) -> spartn.Description:
    tag = timescale.spartn_tag_from_gps(t.week, t.seconds, system, LEAP if system == System.GLONASS else None)
    return spartn.Description(subtype, 1, tag, stream.solution, stream.processor, None, None, None, None)


def sat(number, *, dnu=False, continuity=7, orbit=None, clock=None, phase=None, code=None):
    """A satellite block: orbit (iode, r, a, c, SF021 code); clock (SF022 code, m, SF024);
    phase [(position, fix, SF015 code, m)]; code [(position, m)]."""
    return dict(number=number, dnu=dnu, continuity=continuity, orbit=orbit, clock=clock, phase=phase, code=code)


def ocb(system: System, t: GnssTime, sats, *, siou=SIOU, eos=False, ephemeris_type=0,
        stream=STREAM) -> bytes:
    """[spec: SPARTN 6.5.1, Tables 6.3 to 6.13] an OCB frame."""
    mask, blocks = 0, []
    for s in sorted(sats, key=lambda s: s["number"]):
        mask |= 1 << (SAT_BITS[system] - s["number"])
        if s["dnu"]:
            blocks.append(spartn.OcbSatellite(1, None, None, None, None, None))
            continue
        orbit = clock = bias = None
        if s["orbit"] is not None:
            iode, r, a, c, y = s["orbit"]
            orbit = spartn.Orbit(iode, enc(20, r), enc(20, a), enc(20, c), y)
        if s["clock"] is not None:
            k, m, ure = s["clock"]
            clock = spartn.Clock(k, enc(20, m), ure)
        if s["phase"] is not None or s["code"] is not None:
            phase, code, width = sorted(s["phase"] or ()), sorted(s["code"] or ()), BIAS_BITS[system]
            bias = spartn.Bias(0, sum(1 << (width - 1 - p) for p, *_ in phase),
                               tuple(spartn.PhaseBias(f, k, enc(20, v)) for _, f, k, v in phase),
                               0, sum(1 << (width - 1 - p) for p, _ in code), tuple(enc(29, v) for _, v in code))
        present = (orbit is not None) << 2 | (clock is not None) << 1 | (bias is not None)
        blocks.append(spartn.OcbSatellite(0, present, s["continuity"], orbit, clock, bias))
    body = spartn.Ocb(siou, int(eos), 0, 1, 0, ephemeris_type, 0, mask, tuple(blocks))
    return spartn.encode(spartn.Message(0, 0, 0, description(system, t, int(system), stream), body))


def gad(t: GnssTime, *, siou=SIOU, aiou=AIOU, lat=-50.0, lon=-164.0, stream=STREAM) -> bytes:
    """[spec: SPARTN 6.5.3, Table 6.24] one 3 × 4 area of 1° spacing, north-west corner (lat, lon)."""
    a = spartn.GadArea(1, enc(32, lat), enc(33, lon), 2, 3, enc(36, 1.0), enc(37, 1.0))
    return spartn.encode(spartn.Message(2, 0, 2, description(System.GPS, t, 0, stream),
                                        spartn.Gad(siou, aiou, 0, 0, (a,))))


TROPO = (2.31, (0.1, 0.01, -0.02, 0.004))                     # SF043, SF045 to SF047
TRES = tuple(0.004 * k - 0.02 for k in range(12))             # SF052, 12 nodes
IONO = (10.0, 0.4, -0.8, 0.02)                                # SF057 to SF059
IRES = tuple(0.04 * k - 0.2 for k in range(12))               # SF065


def hpac(system: System, t: GnssTime, numbers=(), *, siou=SIOU, aiou=AIOU, troposphere_block=True,
         tres=TRES, ires=IRES, nodes=12, stream=STREAM) -> bytes:
    """[spec: SPARTN 6.5.2, Tables 6.14 to 6.23] one HPAC area, ID 1, with the grid of ``nodes``."""
    tropo = None
    if troposphere_block:
        tropo = spartn.Troposphere(2, 2, enc(43, TROPO[0]), 0,
                                   spartn.TropoSmall(*(enc(f, c) for f, c in zip((45, 46, 46, 47), TROPO[1]))),
                                   None, 0, tuple(enc(52, r) for r in tres))
    iono = None
    if numbers:
        blocks = tuple(spartn.IonoSatellite(5, 0, spartn.IonoSmall(*(enc(f, c) for f, c in
                                                                      zip((57, 58, 58, 59), IONO))),
                                            None, 1, tuple(enc(65, r) for r in ires)) for _ in numbers)
        mask = sum(1 << (SAT_BITS[system] - n) for n in numbers)
        iono = spartn.Ionosphere(2, 0, mask, blocks)
    a = spartn.HpacArea(1, nodes, 2 if tropo else 0, 2 if iono else 0, tropo, iono)
    body = spartn.Hpac(siou, aiou, 0, 0, (a,))
    return spartn.encode(spartn.Message(1, 0, 2, description(system, t, int(system), stream), body))


GPS_SAT = sat(5, orbit=(33, 0.5, -0.25, 0.124, 10), clock=(7, 1.234, 2),
              phase=[(0, 1, 7, 0.01), (1, 0, 7, -0.03)], code=[(0, 1.2), (1, -0.4)])
GLO_SAT = sat(1, orbit=(13, -0.3, 0.2, 0.1, 0), clock=(7, -0.5, 3), phase=[(0, 1, 7, 0.02)], code=[(0, 0.8)])


def stream(t: GnssTime = TS, gps_sat=GPS_SAT, **grid) -> list[tuple[bytes, GnssTime]]:
    """The synthetic stream: ephemerides, then one SPARTN set at ``t``, each with its
    receipt time; ``grid`` changes the GPS HPAC's residuals."""
    return [(frame_1019(), T - 30.0), (frame_1020(), T - 29.0),
            (gad(t), t + 0.1),
            (hpac(System.GPS, t, (5,), **grid), t + 0.2),
            (hpac(System.GLONASS, t, (1,), troposphere_block=False), t + 0.3),
            (ocb(System.GLONASS, t, [GLO_SAT]), t + 0.4),
            (ocb(System.GPS, t, [gps_sat], eos=True), t + 0.5)]


def feed(engine: Engine, pieces, framer: Framer | None = None, timed: bool = True) -> list:
    """Every piece through one framer and into the engine; the feed results."""
    framer = framer or Framer()
    results = []
    for piece, received in pieces:
        r = received if timed else None
        for message in framer.push("ip", piece, r):
            results.append(engine.feed(message, r))
    return results


def engine_with(pieces=None, config: Config | None = None, leap=True) -> Engine:
    e = Engine(config)
    if leap:
        assert e.feed(LeapSeconds(LEAP)).outcome == FeedOutcome.STORED
    feed(e, stream() if pieces is None else pieces)
    return e


def decoded(frame: bytes) -> spartn.Message:
    f = spartn.parse(frame)
    assert isinstance(f, spartn.Frame)
    return spartn.decode(f)


# ---- the record, term by term ---------------------------------------------------------


def test_every_term_of_a_gps_record_is_the_toolkits_value():
    e = engine_with()
    assert e.status().active is not None
    rec = e.query(G05, T, RECEIVER)
    assert rec.missing == {}, rec.missing
    assert rec.signals == (L1C, L2W)            # every signal with a bias
    assert rec.corrected

    # the inputs, decoded independently of the engine
    eph = e.ephemerides.issues(G05)[0].ephemeris
    block = spartn_ocb.ocb_input(decoded(ocb(System.GPS, TS, [GPS_SAT], eos=True))).satellites[0]
    o, c = block.orbit, block.clock
    code = {x.signal: x.value for x in block.code_biases}
    phase = {x.signal: x.value for x in block.phase_biases}
    hp = spartn_atmosphere.from_hpac(decoded(hpac(System.GPS, TS, (5,)))).areas[0]

    # transmit time: SPARTN 8.8 position, the light time (Ashby 2003 section 2)
    tau = T - rec.t_tx
    state, clock = broadcast_state(eph, rec.t_tx)
    d_o = (o.radial, o.along, o.cross)
    x = ssr.corrected_position(state.position, state.velocity, d_o)
    sagnac = relativity.sagnac(x, RECEIVER)
    rho = math.dist(x, RECEIVER)
    assert abs((rho + sagnac) / SPEED_OF_LIGHT - tau) < 2e-12
    assert 0.06 < tau < 0.09

    s = rec.satellite
    assert s.position == x
    assert s.velocity == state.velocity
    assert s.clock == ssr.corrected_clock(clock, c.correction)              # SPARTN 8.7
    assert s.clock_drift == pytest.approx(eph.af1, abs=1e-20)
    assert s.broadcast.position == state.position and s.broadcast.clock == clock
    assert (s.broadcast.iod, s.broadcast.message, s.broadcast.toe) == (33, 1019, eph.toe)
    assert s.broadcast.age == T - eph.toe

    assert rec.orbit.ecef == ssr.orbit_correction(state.position, state.velocity, d_o)
    assert (rec.orbit.radial, rec.orbit.along, rec.orbit.cross, rec.orbit.iode) == (o.radial, o.along, o.cross, 33)
    assert rec.orbit.yaw == ssr.field(21, 10).value and rec.orbit.age == T - TS
    assert rec.clock_correction.correction == c.correction
    assert rec.clock_correction.ure == ssr.field(24, 2).value and rec.clock_correction.ure_index == 2

    # geometry in the frame of reception (TM-23 B.3)
    geodetic = geodesy.ecef_to_geodetic(RECEIVER)
    arrival = geodesy.rotate_z(x, geodesy.WGS84_OMEGA_E * tau)
    unit, elevation, azimuth = geodesy.line_of_sight(arrival, RECEIVER, geodetic)
    g = rec.geometry
    assert (g.unit, g.elevation, g.azimuth) == (unit, elevation, azimuth)
    assert g.range == pytest.approx(rho, abs=1e-8)
    assert g.range_rate == geodesy.range_rate(x, state.velocity, RECEIVER)[1]
    assert math.degrees(elevation) == pytest.approx(70.7, abs=0.2)

    # modelled terms
    m = rec.modelled
    assert m.earth_rotation == sagnac
    assert m.shapiro == relativity.shapiro(x, RECEIVER)
    assert m.relativity_clock == relativity.clock(state.position, state.velocity)
    sky = astro.sky(T.week, T.seconds, -float(LEAP))
    tide = tides.solid_earth_tide(RECEIVER, sky)
    assert m.tide_displacement == tide
    assert m.solid_tide == -(unit[0] * tide[0] + unit[1] * tide[1] + unit[2] * tide[2])
    east, north, _ = geodesy.enu_basis(geodetic[0], geodetic[1])
    fx, fy, _ = wu.satellite_frame(x, state.velocity, math.radians(rec.orbit.yaw))
    assert m.windup == wu.windup(x, RECEIVER, east, north, fx, fy)

    # biases (SPARTN 8.9)
    assert {k: v.value for k, v in rec.code_bias.items()} == code
    assert {k: v.value for k, v in rec.phase_bias.items()} == phase
    assert rec.phase_bias[L1C].fixed and not rec.phase_bias[L2W].fixed
    assert rec.phase_bias[L1C].discontinuities == 0
    assert rec.frequency == {L1C: L1C.frequency(), L2W: L2W.frequency()}

    # troposphere (SPARTN 8.12, Niell, DESIGN Q7)
    gad_area = area.Area(-50.0, -164.0, 3, 4, 1.0, 1.0)
    lat, lon = math.degrees(geodetic[0]), math.degrees(geodetic[1])
    dlat, dlon = area.offsets(gad_area, lat, lon)
    nodes, weights = area.grid(gad_area, lat, lon)
    tb = hp.troposphere
    poly = area.polynomial(tb.equation_type, tb.coefficients, dlat, dlon)
    grid = area.interpolate(weights, [tb.residuals[n] for n in nodes])
    day = T.second // 86400
    doy = timescale.day_of_year(T.week * 7 + day) + (T.seconds - day * 86400.0) / 86400.0
    hydro = tb.hydrostatic * troposphere.hydrostatic_height_factor(geodetic[0], geodetic[2], doy)
    mapping = troposphere.niell_mapping(geodetic[0], geodetic[2], doy, elevation)
    tr = rec.troposphere
    assert (tr.polynomial, tr.grid, tr.hydrostatic_zenith) == (poly, grid, hydro)
    assert tr.slant == troposphere.slant(hydro, poly, grid, mapping.hydrostatic, mapping.wet)
    assert tr.age == T - TS and tr.area_id == 1

    # ionosphere (SPARTN 8.11)
    ib = hp.ionosphere.satellites[0]
    tec = area.polynomial(hp.ionosphere.equation_type, ib.coefficients, dlat, dlon)
    tec_grid = area.interpolate(weights, [ib.residuals[n] for n in nodes])
    io = rec.ionosphere
    assert (io.polynomial, io.grid) == (tec, tec_grid)
    assert io.code_delay == {sg: ionosphere.delay(tec, sg.frequency(), tec_grid)[0] for sg in (L1C, L2W)}

    # quality (DESIGN 2.5)
    sq = rec.quality.satellite
    assert sq.ure_index == 2 and not sq.fixed and sq.yaw_present and sq.discontinuities == 0
    st = rec.quality.stream
    assert (st.solution_id, st.processor_id, st.siou, st.aiou, st.end_of_set) == (1, 0, SIOU, AIOU, True)


def test_glonass_record_uses_the_channel_and_no_relativistic_clock_term():
    e = engine_with()
    rec = e.query(R01, T, RECEIVER, (R1C,))
    assert rec.missing == {}, rec.missing
    eph = e.ephemerides.issues(R01)[0].ephemeris
    block = spartn_ocb.ocb_input(decoded(ocb(System.GLONASS, TS, [GLO_SAT]))).satellites[0]
    state, clock = broadcast_state(eph, rec.t_tx)
    o = block.orbit
    x = ssr.corrected_position(state.position, state.velocity, (o.radial, o.along, o.cross))
    tau = T - rec.t_tx
    assert abs((math.dist(x, RECEIVER) + relativity.sagnac(x, RECEIVER)) / SPEED_OF_LIGHT - tau) < 2e-12
    assert rec.satellite.position == x
    assert rec.satellite.clock == ssr.corrected_clock(clock, block.clock.correction)
    assert rec.satellite.clock_drift == orbits.glonass_clock_drift(glonass_orbit(eph))
    assert rec.satellite.frequency_channel == 1 and rec.frequency == {R1C: R1C.frequency(1)}
    assert rec.modelled.relativity_clock == 0.0
    assert rec.troposphere is not None          # the area's troposphere, from the GPS HPAC
    assert rec.ionosphere.code_delay[R1C] == ionosphere.delay(rec.ionosphere.polynomial, R1C.frequency(1),
                                                               rec.ionosphere.grid)[0]
    assert rec.code_bias[R1C].value == block.code_biases[0].value


def test_range_rate_is_the_derivative_of_the_range():
    e = engine_with()
    h = 0.5
    before, after = e.query(G05, T - h, RECEIVER), e.query(G05, T + h, RECEIVER)
    rate = e.query(G05, T, RECEIVER).geometry.range_rate
    # d/dt_rx of ρ(t_tx) is ρ̇ (1 − ρ̇/c): 2 mm/s at 800 m/s
    assert (after.geometry.range - before.geometry.range) / (2 * h) == pytest.approx(rate, abs=0.01)


# ---- replay and caching (DESIGN 2.1, 2.8) -----------------------------------------------


def test_a_replay_of_the_same_bytes_gives_identical_records():
    a, b = engine_with(), engine_with()
    for sat_, signals in ((G05, ()), (R01, (R1C,)), (Sat.gps(9), (L1C,))):
        for dt in (0.0, 0.5, 7.0):
            assert a.query(sat_, T + dt, RECEIVER, signals) == b.query(sat_, T + dt, RECEIVER, signals)


def test_a_replay_without_receipt_times_places_everything_by_the_stream():
    """The ephemerides arrive before any time is known: held (T5), released by the
    first SPARTN tag; the record is that of the timed feed, receipt times aside."""
    pieces = stream()
    untimed = Engine()
    untimed.feed(LeapSeconds(LEAP))
    results = feed(untimed, pieces, timed=False)
    assert [r.outcome for r in results[:2]] == [FeedOutcome.HELD, FeedOutcome.HELD]
    assert results[0].reason == FeedReason.NO_TIME
    assert untimed.status().held == 0
    timed = engine_with(pieces)
    for sat_ in (G05, R01):
        x, y = untimed.query(sat_, T, RECEIVER), timed.query(sat_, T, RECEIVER)
        assert x == y


def test_queries_without_a_feed_reuse_the_selection():
    e = engine_with()
    e.query(G05, T, RECEIVER)
    n = e.selections
    for dt in (0.0, 1.0, 10.0):
        e.query(G05, T + dt, RECEIVER)
    assert e.selections == n                    # the choice holds until a block ages (DESIGN 2.8)
    e.query(G05, T + 40.0, RECEIVER)            # past the clock's 30 s: chosen again
    assert e.selections == n + 1
    e.feed(LeapSeconds(LEAP))
    e.query(G05, T + 40.0, RECEIVER)
    assert e.selections == n + 2                # every feed invalidates


def test_query_cost(record_property):
    """DESIGN 2.8: the cost of a query once the epoch's shared work is done."""
    e = engine_with()
    sats = (G05, R01)
    e.query(G05, T, RECEIVER)
    n = 2000
    start = time.perf_counter()
    for k in range(n):
        e.query(sats[k % 2], T, RECEIVER)
    per_query = (time.perf_counter() - start) / n
    record_property("query_seconds", per_query)
    assert per_query < 5e-3


# ---- phase-bias discontinuity (rule C7) -------------------------------------------------


def test_a_phase_bias_gap_longer_than_its_continuity_counts_a_discontinuity():
    short = sat(5, orbit=(33, 0.5, -0.25, 0.124, 10), clock=(7, 1.234, 2), phase=[(0, 1, 2, 0.01)])  # 5 s
    e = engine_with(stream(TS - 8.0, short))
    assert e.query(G05, T, RECEIVER, (L1C,)).phase_bias[L1C].discontinuities == 0
    feed(e, [(ocb(System.GPS, TS, [short], eos=True), TS + 0.5)])
    rec = e.query(G05, T, RECEIVER, (L1C,))
    assert rec.phase_bias[L1C].discontinuities == 1
    assert rec.quality.satellite.discontinuities == 0      # the satellite's own arc is unbroken


# ---- feed results (DESIGN 2.1, 2.2) -----------------------------------------------------


def test_feed_results_say_what_each_input_did():
    e = Engine()
    r = e.feed(rtcm3.decode(rtcm3.parse(frame_1020())), T)
    assert (r.outcome, r.reason) == (FeedOutcome.HELD, FeedReason.NO_LEAP_SECONDS)
    r = e.feed(LeapSeconds(LEAP))
    assert r.outcome == FeedOutcome.STORED and "released" in r.detail
    assert e.ephemerides.glonass_channel(R01) == 1
    r = e.feed(rtcm3.decode(rtcm3.parse(rtcm(rtcm3.Msg1005))), T)
    assert (r.outcome, r.reason) == (FeedOutcome.REJECTED, FeedReason.NOT_AN_INPUT)
    bad = rtcm(rtcm3.Msg1019, satellite=5, week=2295 % 1024, iode=33, iodc=34, sqrt_a=q(5153.7, -19))
    r = e.feed(rtcm3.decode(rtcm3.parse(bad)), T)
    assert (r.outcome, r.reason) == (FeedOutcome.REJECTED, FeedReason.INCONSISTENT)
    r = e.feed(OcbHeaderInput(0, STREAM, None, SIOU, True, 0), T)
    assert r.outcome in (FeedOutcome.COUNTED, FeedOutcome.REJECTED)
    with pytest.raises(TypeError):
        e.feed(object())


# ---- every missing reason -----------------------------------------------------------------


def _engine(pieces, config=None, leap=True):
    return engine_with(pieces, config, leap)


def _unhealthy():
    e = Engine()
    e.feed(LeapSeconds(LEAP))
    feed(e, stream()[:1])
    eph = e.ephemerides.issues(G05)[0].ephemeris
    other = Engine()
    other.feed(dataclasses.replace(eph, health=1), T)
    return other.query(G05, T, RECEIVER)


def _not_converged():
    e = Engine()
    feed(e, stream()[:1])
    eph: KeplerEphemeris = e.ephemerides.issues(G05)[0].ephemeris
    other = Engine()
    other.feed(dataclasses.replace(eph, e=1.5), T)
    return other.query(G05, T, RECEIVER)


def _iode_change():
    first = sat(5, orbit=(33, 0.5, -0.25, 0.124, 10), clock=(7, 1.234, 2))
    later = sat(5, clock=(0, 1.3, 2))                # SF022 code 0: the IODE held for 0 s
    pieces = stream(TS - 20.0, first) + [(ocb(System.GPS, TS, [later], eos=True), TS + 0.5)]
    return _engine(pieces).query(G05, T, RECEIVER)


def _engine_query(pieces=None, sat_=G05, t=T, receiver=RECEIVER, signals=(), config=None, leap=True):
    return _engine(stream() if pieces is None else pieces, config, leap).query(sat_, t, receiver, signals)


def _receiver(lat, lon):
    return geodesy.geodetic_to_ecef((math.radians(lat), math.radians(lon), HEIGHT))


SCENARIOS = {
    "no_ephemeris": (lambda: _engine_query(stream()[2:]), (Term.SATELLITE, None), Missing.NO_EPHEMERIS),
    "outside_validity": (lambda: _engine_query(t=T + 4 * 3600.0), (Term.SATELLITE, None),
                         Missing.OUTSIDE_VALIDITY),
    "iode_mismatch": (lambda: _engine_query(stream(gps_sat=sat(5, orbit=(34, 0.5, 0, 0, 0), clock=(7, 1.0, 2)))),
                      (Term.ORBIT, None), Missing.IODE_MISMATCH),
    "ephemeris_type": (lambda: _engine_query(stream()[:4] + [(ocb(System.GPS, TS, [GPS_SAT], eos=True,
                                                                  ephemeris_type=1), TS + 0.5)]),
                       (Term.ORBIT, None), Missing.EPHEMERIS_TYPE),
    "unhealthy": (_unhealthy, (Term.SATELLITE, None), Missing.UNHEALTHY),
    "stale": (lambda: _engine_query(t=T + 40.0), (Term.CLOCK_CORRECTION, None), Missing.STALE),
    "outside_area": (lambda: _engine_query(receiver=_receiver(-46.4, -157.7)), (Term.TROPOSPHERE, None),
                     Missing.OUTSIDE_AREA),
    "latitude_limit": (lambda: _engine_query(receiver=_receiver(-86.0, 180.0)), (Term.TROPOSPHERE, None),
                       Missing.LATITUDE_LIMIT),
    "dnu": (lambda: _engine_query(stream(gps_sat=sat(5, dnu=True))), (Term.ORBIT, None), Missing.DNU),
    "siou_mismatch": (lambda: _engine_query(stream()[:2] + [(gad(TS, siou=SIOU + 1), TS + 0.1),
                                                            (hpac(System.GPS, TS, (5,), siou=SIOU + 1), TS + 0.2),
                                                            (ocb(System.GPS, TS, [GPS_SAT], eos=True), TS + 0.5)]),
                      (Term.TROPOSPHERE, None), Missing.SIOU_MISMATCH),
    "aiou_mismatch": (lambda: _engine_query(stream()[:3] + [(hpac(System.GPS, TS, (5,), aiou=AIOU + 1), TS + 0.2),
                                                            (ocb(System.GPS, TS, [GPS_SAT], eos=True), TS + 0.5)]),
                      (Term.TROPOSPHERE, None), Missing.AIOU_MISMATCH),
    "no_corrections": (lambda: _engine_query(stream()[:2]), (Term.ORBIT, None), Missing.NO_CORRECTIONS),
    "iode_change": (_iode_change, (Term.ORBIT, None), Missing.IODE_CHANGE),
    "no_bias": (lambda: _engine_query(signals=(L1C, L5Q)), (Term.CODE_BIAS, L5Q), Missing.NO_BIAS),
    "no_atmosphere": (lambda: _engine_query(stream()[:3] + stream()[5:]), (Term.TROPOSPHERE, None),
                      Missing.NO_ATMOSPHERE),
    "invalid_grid_node": (lambda: _engine_query(stream(tres=(None,) * 12)), (Term.TROPOSPHERE_GRID, None),
                          Missing.INVALID_GRID_NODE),
    "grid_mismatch": (lambda: _engine_query(stream(ires=IRES[:6], tres=TRES[:6], nodes=6)),
                      (Term.IONOSPHERE_GRID, None), Missing.GRID_MISMATCH),
    "no_leap_seconds": (lambda: _engine_query(sat_=R01, leap=False), (Term.SATELLITE, None),
                        Missing.NO_LEAP_SECONDS),
    "below_elevation": (lambda: _engine_query(config=Config(elevation_cutoff=math.radians(80.0))),
                        (Term.TROPOSPHERE, None), Missing.BELOW_ELEVATION),
    "not_converged": (_not_converged, (Term.GEOMETRY, None), Missing.NOT_CONVERGED),
    "no_frequency": (lambda: _engine_query(sat_=Sat.glonass(2), signals=(R1C,)), (Term.FREQUENCY, R1C),
                     Missing.NO_FREQUENCY),
}


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_missing_reason(name):
    make, key, reason = SCENARIOS[name]
    rec = make()
    assert rec.missing.get(key) == reason, rec.missing
    term, signal = key
    field_name = term.value.split(".")[0]
    if field_name in ("code_bias", "phase_bias", "frequency"):
        assert signal not in getattr(rec, field_name)
    elif "." not in term.value:
        assert getattr(rec, field_name) is None


def test_every_missing_reason_is_produced():
    assert {reason for _, _, reason in SCENARIOS.values()} == set(Missing)
    assert set(SCENARIOS) == {m.value for m in Missing}


def test_tide_without_leap_seconds_is_missing_but_the_rest_is_there():
    rec = _engine_query(leap=False)
    assert rec.missing == {(Term.SOLID_TIDE, None): Missing.NO_LEAP_SECONDS}
    assert rec.modelled.solid_tide is None and rec.modelled.windup is not None


def test_broadcast_state_without_corrections_is_marked_uncorrected():
    rec = _engine_query(stream()[:2])
    assert not rec.corrected and rec.satellite is not None
    eph = None
    state, clock = None, None
    tau = T - rec.t_tx
    e = Engine()
    feed(e, stream()[:1])
    eph = e.ephemerides.issues(G05)[0].ephemeris
    state, clock = broadcast_state(eph, rec.t_tx)
    assert rec.satellite.position == state.position and rec.satellite.clock == clock
    assert abs((math.dist(state.position, RECEIVER) + relativity.sagnac(state.position, RECEIVER))
               / SPEED_OF_LIGHT - tau) < 2e-12


# ---- streams by class and corrections with rates (DESIGN 3.3, S10 and C11) ------------------

B_STREAM = StreamId(2, 7)
B_GPS_SAT = sat(5, orbit=(33, 0.6, -0.2, 0.1, 10), clock=(7, 2.0, 2), code=[(0, 1.0), (1, -0.3)])


def test_choose_takes_the_named_stream_and_not_the_active_one():
    """Rule S10: a named satellite stream is used even when S2 elects another."""
    pieces = stream() + [(ocb(System.GPS, TS + 10.0, [B_GPS_SAT], eos=True, stream=B_STREAM), TS + 10.5)]
    e = engine_with(pieces)
    assert e.streams.active == B_STREAM
    t = T + 10.0
    named = choose(e.streams, e.corrections, e.ephemerides, G05, t, MaxAge(), stream=STREAM)
    active = choose(e.streams, e.corrections, e.ephemerides, G05, t, MaxAge())
    assert named.stream == STREAM and named.selection.stream == STREAM
    assert named.selection.clock.value.correction == pytest.approx(1.234, abs=1e-9)
    assert active.stream == B_STREAM and active.selection.clock.value.correction == pytest.approx(2.0, abs=1e-9)
    quiet = choose(e.streams, e.corrections, e.ephemerides, G05, t, MaxAge(), stream=StreamId(9, 9))
    assert quiet.selection is None and quiet.corrections_missing == Missing.NO_CORRECTIONS
    assert quiet.broadcast is not None                         # the broadcast state, no other stream


def test_transmit_evaluates_an_orbit_correction_with_rates_at_t_tx():
    """Rule C11: X = X_b − δX(δO + δȮ·(t_tx − t0)); without rates, δO as sent."""
    e = engine_with()
    choice = e._choice(G05, T)
    o = choice.selection.orbit.value
    t0 = TS + 2.5
    rated = dataclasses.replace(o, rate=(0.01, -0.02, 0.03), reference=t0)
    tx = transmit(choice.broadcast, rated, RECEIVER, T)
    d_o = ssr.orbit_at((o.radial, o.along, o.cross), rated.rate, tx.t_tx - t0)
    state, _ = choice.broadcast.state(tx.t_tx)
    assert tx.orbit == d_o
    assert tx.position == ssr.corrected_position(state.position, state.velocity, d_o)
    assert tx.correction == ssr.orbit_correction(state.position, state.velocity, d_o)
    assert abs((tx.range + tx.sagnac) / SPEED_OF_LIGHT - (T - tx.t_tx)) < 2e-12
    plain = transmit(choice.broadcast, o, RECEIVER, T)
    assert plain.orbit == (o.radial, o.along, o.cross)
    state, _ = choice.broadcast.state(plain.t_tx)
    assert plain.position == ssr.corrected_position(state.position, state.velocity, plain.orbit)
    assert plain.position != tx.position


# ---- the atmosphere providers in the engine (DESIGN 3.6, rule A9) -----------------------------


def _b_atmosphere(t: GnssTime) -> list[tuple[bytes, GnssTime]]:
    """Stream B's GAD and HPAC only (B sends no OCB here)."""
    return [(gad(t, stream=B_STREAM), t + 0.1),
            (hpac(System.GPS, t, (5,), stream=B_STREAM), t + 0.2)]


def _a_satellites(t: GnssTime) -> list[tuple[bytes, GnssTime]]:
    """Stream A's (STREAM's) OCB only, for G05."""
    return [(ocb(System.GPS, t, [GPS_SAT], eos=True), t + 0.5)]


TWO_STREAMS = Config(satellite_stream=STREAM, atmosphere=(SpartnSource(), SpartnSource(B_STREAM), ModelSource()))


def test_satellite_corrections_from_one_stream_and_atmosphere_from_another():
    """Rules S10, A9, A10, A11: A gives orbit, clock and biases, B the atmosphere."""
    pieces = stream()[:2] + _b_atmosphere(TS) + _a_satellites(TS)
    e = engine_with(pieces, TWO_STREAMS)
    rec = e.query(G05, T, RECEIVER)
    assert rec.corrected and rec.quality.stream.solution_id == STREAM.solution
    assert rec.troposphere.origin == Origin("spartn:2/7", Kind.CORRECTION, B_STREAM)
    assert rec.ionosphere.origin == Origin("spartn:2/7", Kind.CORRECTION, B_STREAM)
    assert rec.ionosphere.with_biases is False                   # A11: not the biases' stream
    assert rec.troposphere.sigma == ssr.field(42, 2).value and rec.ionosphere.sigma == ssr.field(55, 5).value
    # B's values are those its own engine gives, A's corrections those of A alone.
    alone = engine_with(stream()[:2] + _b_atmosphere(TS) + [(ocb(System.GPS, TS, [GPS_SAT], eos=True,
                                                                   stream=B_STREAM), TS + 0.5)])
    b = alone.query(G05, T, RECEIVER)
    assert b.ionosphere.with_biases is True
    for name in ("slant", "hydrostatic_zenith", "polynomial", "grid", "quality", "area_id", "age", "wet_zenith"):
        assert getattr(rec.troposphere, name) == getattr(b.troposphere, name), name
    assert rec.ionosphere.code_delay == b.ionosphere.code_delay
    assert e.status().satellite_stream == STREAM
    assert e.status().atmosphere == ("spartn", "spartn:2/7", "unb3m")


def test_when_the_atmosphere_stream_goes_quiet_the_walk_falls_to_the_model():
    pieces = stream()[:2] + _b_atmosphere(TS) + _a_satellites(TS) + _a_satellites(TS + 125.0)
    e = engine_with(pieces, TWO_STREAMS)
    t = T + 125.0                                              # B's items are 127 s old (max 120 s)
    rec = e.query(G05, t, RECEIVER)
    assert rec.corrected
    tr = rec.troposphere
    assert tr.origin == Origin("unb3m", Kind.MODEL) and tr.sigma == unb3m.TOTAL_SIGMA
    assert (tr.polynomial, tr.grid, tr.quality, tr.area_id, tr.age) == (None,) * 5
    lat, _, h = geodesy.ecef_to_geodetic(RECEIVER)
    z = unb3m.zenith_delays(lat, h, e._epoch.doy)
    m = troposphere.niell_mapping(lat, h, e._epoch.doy, rec.geometry.elevation)
    assert (tr.hydrostatic_zenith, tr.wet_zenith) == (z.hydrostatic, z.wet)
    assert tr.slant == troposphere.slant(z.hydrostatic, z.wet, 0.0, m.hydrostatic, m.wet)
    assert (Term.TROPOSPHERE_GRID, None) not in rec.missing
    # The ionosphere has no model: the reason is that of the first provider with the term,
    # the satellite stream's own atmosphere, which A never sent.
    assert rec.ionosphere is None and rec.why(Term.IONOSPHERE) == Missing.NO_ATMOSPHERE
    with pytest.raises(MissingTermError):
        rec.pseudorange_model(L1C)
    # Before B went quiet, B answered.
    assert e.query(G05, T, RECEIVER).troposphere.origin.provider == "spartn:2/7"


class _Sensor:
    """A caller's troposphere: zenith delays measured at the site."""

    kind = Kind.CORRECTION
    terms = frozenset({Term.TROPOSPHERE})
    name = "sensor"

    def troposphere(self, where):
        return Zenith(2.21, 0.13, Origin(self.name, self.kind), where.t, 0.0, 0.01)

    def ionosphere(self, where, sat, look):
        return Absent(Missing.NO_ATMOSPHERE)


def test_a_callers_provider_answers_first_when_ranked_first():
    config = Config(atmosphere=(CallerSource("sensor"), SpartnSource()))
    e = Engine(config, providers=dict(sensor=_Sensor()))
    e.feed(LeapSeconds(LEAP))
    feed(e, stream())
    rec = e.query(G05, T, RECEIVER)
    tr = rec.troposphere
    lat, _, h = geodesy.ecef_to_geodetic(RECEIVER)
    m = troposphere.niell_mapping(lat, h, e._epoch.doy, rec.geometry.elevation)
    assert tr.origin == Origin("sensor", Kind.CORRECTION) and tr.sigma == 0.01
    assert tr.slant == troposphere.slant(2.21, 0.13, 0.0, m.hydrostatic, m.wet)
    assert (tr.hydrostatic_zenith, tr.wet_zenith, tr.age, tr.polynomial) == (2.21, 0.13, 0.0, None)
    assert rec.ionosphere.origin == Origin("spartn", Kind.CORRECTION, STREAM) and rec.ionosphere.with_biases
    assert rec.why(Term.TROPOSPHERE_GRID) is None
    with pytest.raises(ValueError):
        Engine(config)                                         # the CallerSource has no object
    with pytest.raises(ValueError):
        Engine(Config(), providers=dict(sensor=_Sensor()))     # the object has no CallerSource


def test_a_clock_polynomial_and_orbit_rates_are_evaluated_at_t_tx():
    """Rule C11 end to end, from a typed OCB input that carries rates."""
    e = engine_with(stream()[:-1])
    typed = spartn_ocb.ocb_input(decoded(ocb(System.GPS, TS, [GPS_SAT], eos=True)))
    s = typed.satellites[0]
    t0 = TS + 2.5
    rated = dataclasses.replace(typed, satellites=(dataclasses.replace(
        s, orbit=dataclasses.replace(s.orbit, rate=(0.001, -0.004, 0.002), reference=t0),
        clock=dataclasses.replace(s.clock, rate=0.01, acceleration=-2e-4, reference=t0)),))
    assert e.feed(rated, TS + 0.5).outcome == FeedOutcome.STORED
    rec = e.query(G05, T, RECEIVER)
    plain = engine_with().query(G05, T, RECEIVER)
    dt = rec.t_tx - t0
    c = rated.satellites[0].clock
    d_c = ssr.clock_at(c.correction, 0.01, -2e-4, dt)
    assert rec.clock_correction.correction == d_c and d_c != c.correction
    assert (rec.clock_correction.rate, rec.clock_correction.acceleration, rec.clock_correction.reference) == \
        (0.01, -2e-4, t0)
    eph = e.ephemerides.issues(G05)[0].ephemeris
    state, clock = broadcast_state(eph, rec.t_tx)
    assert rec.satellite.clock == ssr.corrected_clock(clock, d_c)
    drift = plain.satellite.clock_drift
    assert rec.satellite.clock_drift == ssr.corrected_clock(
        e._choice(G05, T).broadcast.drift(rec.t_tx), ssr.clock_rate_at(0.01, -2e-4, dt))
    assert rec.satellite.clock_drift == pytest.approx(drift - (0.01 - 4e-4 * dt) / SPEED_OF_LIGHT, abs=1e-18)
    o = rated.satellites[0].orbit
    d_o = ssr.orbit_at((o.radial, o.along, o.cross), o.rate, dt)
    assert (rec.orbit.radial, rec.orbit.along, rec.orbit.cross) == d_o
    assert rec.orbit.ecef == ssr.orbit_correction(state.position, state.velocity, d_o)
    assert rec.satellite.position == ssr.corrected_position(state.position, state.velocity, d_o)
    assert (rec.orbit.rate, rec.orbit.reference) == (o.rate, t0) and rec.orbit.age == T - TS
    assert (plain.orbit.rate, plain.clock_correction.reference) == (None, None)


# ---- the default configuration changes nothing (DESIGN 3.6, the default list) ---------------
#
# ``_Before`` is the engine with the atmosphere code of commit 1e9da6e, the last
# before the seam, written out here so that the comparison needs no file of
# that time: the area once per epoch and issue, the store's zenith and
# ionosphere, the height factor at the receiver, the slant from the polynomial
# and grid. Everything else is the engine's own code, which for a SPARTN
# stream (no rates) passes the values as sent, as the term-by-term tests above
# check against the toolkit. Every record field that existed before the seam
# must be equal; the fields the seam added are checked for what they say.

_ADDED = {"wet_zenith", "origin", "sigma", "with_biases", "rate", "acceleration", "reference"}


class _Before(Engine):
    def _atmosphere(self, choice, sat, t_rx, rcv, elevation, azimuth, frequency, epoch, missing):
        selection = choice.selection
        if selection is None:
            missing[(Term.TROPOSPHERE, None)] = missing[(Term.IONOSPHERE, None)] = choice.corrections_missing
            return None, None
        area = self.atmosphere.area(selection.stream, selection.issue, rcv.latitude_deg, rcv.longitude_deg, t_rx)
        if isinstance(area, Absent):
            missing[(Term.TROPOSPHERE, None)] = missing[(Term.IONOSPHERE, None)] = area.reason
            return None, None
        return (self._troposphere_before(area, t_rx, rcv, elevation, epoch, missing),
                self._ionosphere_before(area, sat, t_rx, frequency, missing))

    def _troposphere_before(self, area, t_rx, rcv, elevation, epoch, missing):
        zenith = self.atmosphere.troposphere(area, t_rx)
        if isinstance(zenith, Absent):
            missing[(Term.TROPOSPHERE, None)] = zenith.reason
            return None
        if not elevation > 0.0:
            missing[(Term.TROPOSPHERE, None)] = Missing.BELOW_ELEVATION
            return None
        lat, _, height = rcv.geodetic
        mapping = troposphere.niell_mapping(lat, height, epoch.doy, elevation)
        hydrostatic = zenith.hydrostatic
        if self.config.troposphere_height_reduction:
            hydrostatic = hydrostatic * troposphere.hydrostatic_height_factor(lat, height, epoch.doy)
        if zenith.grid is None:
            missing[(Term.TROPOSPHERE_GRID, None)] = zenith.grid_missing
        slant = troposphere.slant(hydrostatic, zenith.polynomial, zenith.grid or 0.0,
                                  mapping.hydrostatic, mapping.wet)
        return TroposphereTerm(slant, hydrostatic, mapping.hydrostatic, mapping.wet, zenith.polynomial,
                               zenith.grid, zenith.quality, zenith.area_id, zenith.age)

    def _ionosphere_before(self, area, sat, t_rx, frequency, missing):
        tec = self.atmosphere.ionosphere(area, sat, t_rx)
        if isinstance(tec, Absent):
            missing[(Term.IONOSPHERE, None)] = tec.reason
            return None
        if tec.grid is None:
            missing[(Term.IONOSPHERE_GRID, None)] = tec.grid_missing
        grid = tec.grid or 0.0
        delays = {s: ionosphere.delay(tec.polynomial, f, grid)[0] for s, f in frequency.items()}
        for (term, s) in list(missing):
            if term == Term.FREQUENCY:
                missing[(Term.IONOSPHERE, s)] = Missing.NO_FREQUENCY
        return IonosphereTerm(tec.slant, tec.polynomial, tec.grid, delays, tec.quality, tec.area_id, tec.age)


def _old_fields(value):
    """The value with every field the seam added left out, recursively."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {f.name: _old_fields(getattr(value, f.name)) for f in dataclasses.fields(value)
                if f.name not in _ADDED}
    if isinstance(value, (dict, MappingProxyType)):
        return {k: _old_fields(v) for k, v in value.items()}
    if isinstance(value, tuple):
        return tuple(_old_fields(v) for v in value)
    return value


_SCENES = {
    "full": (None, None, RECEIVER),
    "no corrections": (stream()[:2], None, RECEIVER),
    "siou mismatch": (stream()[:2] + [(gad(TS, siou=SIOU + 1), TS + 0.1),
                                      (hpac(System.GPS, TS, (5,), siou=SIOU + 1), TS + 0.2),
                                      (ocb(System.GPS, TS, [GPS_SAT], eos=True), TS + 0.5)], None, RECEIVER),
    "aiou mismatch": (stream()[:3] + [(hpac(System.GPS, TS, (5,), aiou=AIOU + 1), TS + 0.2),
                                      (ocb(System.GPS, TS, [GPS_SAT], eos=True), TS + 0.5)], None, RECEIVER),
    "no leap seconds": ("no leap", None, RECEIVER),
    "no atmosphere": (stream()[:3] + stream()[5:], None, RECEIVER),
    "invalid grid": (stream(tres=(None,) * 12), None, RECEIVER),
    "grid mismatch": (stream(ires=IRES[:6], tres=TRES[:6], nodes=6), None, RECEIVER),
    "outside area": (None, None, _receiver(-46.4, -157.7)),
    "latitude limit": (None, None, _receiver(-86.0, 180.0)),
    "cutoff": (None, Config(elevation_cutoff=math.radians(80.0)), RECEIVER),
    "no reduction": (None, Config(troposphere_height_reduction=False), RECEIVER),
    "below horizon": (None, Config(elevation_cutoff=-math.pi / 2 + 1e-9), _receiver(51.4, 17.3)),
}


@pytest.mark.parametrize("scene", sorted(_SCENES))
def test_the_default_configuration_gives_the_records_of_the_engine_before_the_seam(scene):
    pieces, config, receiver = _SCENES[scene]
    leap = pieces != "no leap"
    pieces = stream() if pieces in (None, "no leap") else pieces
    now, before = Engine(config), _Before(config)
    for e in (now, before):
        if leap:
            e.feed(LeapSeconds(LEAP))
        feed(e, pieces)
    for sat_, signals in ((G05, ()), (R01, (R1C,)), (G05, (L1C, L2W, L5Q)), (Sat.gps(9), (L1C,)),
                          (Sat.glonass(2), (R1C,))):
        for dt in (0.0, 0.5, 7.0, 40.0):
            a, b = now.query(sat_, T + dt, receiver, signals), before.query(sat_, T + dt, receiver, signals)
            assert _old_fields(a) == _old_fields(b), (scene, sat_, dt)
            if a.troposphere is not None:
                assert a.troposphere.origin == Origin("spartn", Kind.CORRECTION, STREAM)
                assert a.troposphere.wet_zenith == a.troposphere.polynomial + (a.troposphere.grid or 0.0)
            if a.ionosphere is not None:
                assert a.ionosphere.origin == Origin("spartn", Kind.CORRECTION, STREAM)
                assert a.ionosphere.with_biases is True        # every scene's satellites have biases
            if a.clock_correction is not None:
                # No rates: orbit, clock and drift are the values as sent and the
                # broadcast drift, independently of the engine's code for C11.
                choice = now._choice(sat_, T + dt)
                o, c = choice.selection.orbit.value, choice.selection.clock.value
                assert (a.clock_correction.rate, a.orbit.rate) == (None, None)
                assert (a.orbit.radial, a.orbit.along, a.orbit.cross) == (o.radial, o.along, o.cross)
                assert a.clock_correction.correction == c.correction
                assert a.satellite.clock_drift == choice.broadcast.drift(a.t_tx)
