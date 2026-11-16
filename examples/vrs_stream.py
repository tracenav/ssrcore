"""Turn a synthetic input into the RTCM stream of a virtual base station.

The input is built here with gnsswire from chosen values and goes through the
wire as a receiver and a correction service would deliver it:

- RTCM 1019 ephemerides of seven GPS satellites (the public broadcast record
  of G05 on 2024-01-01, and six orbits made from it by turning its mean
  anomaly and node, so that all seven are above the horizon of the site) and
  the 1020 of GLONASS R01;
- SPARTN GAD and HPAC every 30 s and OCB every 5 s, whose clock corrections
  change a little at every update;
- events: G05 changes its ephemeris issue (a new 1019 with IODE 34, then
  corrections naming it), the L2W phase bias of G09 starts a new arc
  (SF015 = 0, its value moved), and the rover drives 6 km east, so that the
  station moves once.

Every piece goes through ``Framer`` into an ``Engine`` with its receipt time;
a ``VirtualBase`` on the engine gives, at every epoch, the messages of the
station, encoded with gnsswire into one RTCM stream.

Run it with ``uv run python examples/vrs_stream.py [seconds] [--write PATH]``:
it prints what the stream holds and, with ``--write``, writes the stream to
PATH.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from dataclasses import dataclass, field

from gnsswire import rtcm3, spartn

from ssrcore import Engine, Framer, GnssTime, LeapSeconds, System
from ssrcore.model.types import VrsConfig
from ssrcore.toolkit import geodesy, timescale
from ssrcore.vrs.station import VirtualBase

START = GnssTime(2295, 87400.0)    # 2024-01-01 00:16:22 GPS time, the first epoch
LEAP = 18                           # GPS − UTC on that day (IERS Bulletin C)
SITE = (math.radians(-51.4), math.radians(-162.7), 100.0)
PI_ICD = 3.1415926535898            # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2

# [spec: IGS BRDC 2024-01-01, G05 at toe 93 600 s] the public broadcast elements.
G05 = dict(toc=93600.0, af0=-1.562642864883e-04, af1=-1.477928890381e-12, toe=93600.0, sqrt_a=5153.783403397,
           e=0.005790898343548, m0=-2.494229256397, dn=4.453756945425e-09, omega=1.235180510432,
           omega0=-0.6299845609285, omega_dot=-8.155696860368e-09, i0=0.9693533956084,
           idot=-3.928735076157e-11, cuc=-1.747161149979e-06, cus=3.31737101078e-06, crc=319.34375,
           crs=-34.125, cic=1.098960638046e-07, cis=4.097819328308e-08)
# [design: synthetic satellites] PRN: (turn of the mean anomaly, turn of the node), rad;
# each is above 15° at the site for ten minutes from START.
GPS_ORBITS = {5: (0.0, 0.0), 9: (-0.35, -0.5), 12: (-0.7, 1.0), 17: (0.0, -1.0), 20: (0.0, 1.5),
              24: (-0.35, 1.0), 28: (-1.05, 1.0)}
# [spec: IGS BRDC 2024-01-01, R01 at tb index 13] the public broadcast state.
R01 = dict(position=(-11175834.47266, -9508662.597656, -20866544.92188),
           velocity=(2857.489585876, -39.36195373535, -1511.75403595),
           lunisolar=(-2.793967723846e-06, -1.862645149231e-06, 0.0), tau_n=-7.346831262112e-05,
           gamma_n=9.094947017729e-13)


@dataclass(frozen=True)
class Scenario:
    """The events of the synthetic input, in seconds from ``START``.

    [design: synthetic input for DESIGN 5.1; each field one event the
    virtual base station's tests need]
    """

    seconds: int = 300
    issue_change: float | None = 100.0        # G05: 1019 IODE 34 received; corrections name it 20 s later
    bias_arc: float | None = 200.0            # G09 L2W: SF015 = 0 and the value moved by 0.1 m
    move: float | None = 150.0                # the rover is 6 km east from here on
    siou_change: float | None = None          # every message from here on carries SIOU 4
    float_from: float | None = None           # G12's phase biases flagged float from here on
    dnu: tuple[float, float] | None = None    # G17 flagged DNU at the first time, back at the second
    rate: float = 1.0
    msm: int = 7
    float_phase: bool = True
    extra: dict = field(default_factory=dict)


def _q(value: float, exponent: int, semicircle: bool = False) -> int:
    # [spec: RTCM Table 3.5-21] the integer of a value at resolution 2^exponent
    return round((value / PI_ICD if semicircle else value) / 2.0 ** exponent)


def _angle(value: float) -> float:
    # [design: a turned angle back into [−π, π), the range of the semicircle fields]
    return (value + math.pi) % (2.0 * math.pi) - math.pi


def frame_1019(prn: int, iode: int = 33) -> bytes:
    """[spec: RTCM Table 3.5-21] the 1019 of a GPS satellite of ``GPS_ORBITS``."""
    dm0, dnode = GPS_ORBITS[prn]
    r = G05
    fields = {f.name: 0 for f in dataclasses.fields(rtcm3.Msg1019) if f.name != "extension"}
    fields.update(
        number=1019, satellite=prn, week=2295 % 1024, iode=iode, iodc=iode, toe=int(r["toe"] // 16),
        toc=int(r["toc"] // 16), af0=_q(r["af0"], -31), af1=_q(r["af1"], -43), sqrt_a=_q(r["sqrt_a"], -19),
        e=_q(r["e"], -33), m0=_q(_angle(r["m0"] + dm0), -31, True), delta_n=_q(r["dn"], -43, True),
        omega=_q(r["omega"], -31, True), omega0=_q(_angle(r["omega0"] + dnode), -31, True),
        omega_dot=_q(r["omega_dot"], -43, True), i0=_q(r["i0"], -31, True), idot=_q(r["idot"], -43, True),
        cuc=_q(r["cuc"], -29), cus=_q(r["cus"], -29), cic=_q(r["cic"], -29), cis=_q(r["cis"], -29),
        crc=_q(r["crc"], -5), crs=_q(r["crs"], -5))
    return rtcm3.encode(rtcm3.Msg1019(**fields))


def frame_1020() -> bytes:
    """[spec: RTCM Table 3.5-22] the 1020 of R01, frequency channel 1."""
    r, km = R01, 1000.0
    fields = {f.name: 0 for f in dataclasses.fields(rtcm3.Msg1020) if f.name != "extension"}
    fields.update(number=1020, satellite=1, channel=1 + 7, tb=13,
                  x=_q(r["position"][0] / km, -11), y=_q(r["position"][1] / km, -11),
                  z=_q(r["position"][2] / km, -11), x_dot=_q(r["velocity"][0] / km, -20),
                  y_dot=_q(r["velocity"][1] / km, -20), z_dot=_q(r["velocity"][2] / km, -20),
                  x_ddot=_q(r["lunisolar"][0] / km, -30), y_ddot=_q(r["lunisolar"][1] / km, -30),
                  z_ddot=_q(r["lunisolar"][2] / km, -30), tau=_q(r["tau_n"], -30), gamma=_q(r["gamma_n"], -40))
    return rtcm3.encode(rtcm3.Msg1020(**fields))


# ---- SPARTN (v2.0.3 6.3, Table 6.2) ---------------------------------------------------

# [spec: SPARTN Table 6.2] SF: (resolution, range minimum, offset of the notes)
_SF = {20: (0.002, -16.382, 0.0), 29: (0.02, -20.46, 0.0), 32: (0.1, -90.0, 0.0), 33: (0.1, -180.0, 0.0),
       36: (0.1, 0.1, 0.0), 37: (0.1, 0.1, 0.0), 43: (0.004, -0.508, 2.3), 45: (0.004, -0.252, 0.252),
       46: (0.001, -0.063, 0.0), 47: (0.0002, -0.051, 0.0), 57: (0.04, -81.88, 0.0), 58: (0.008, -16.376, 0.0),
       59: (0.002, -8.190, 0.0)}
# [spec: SPARTN Table 6.2 SF011, SF012 size 0; SF025, SF026 size 0] mask widths.
_SAT_BITS = {System.GPS: 32, System.GLONASS: 24}
_BIAS_BITS = {System.GPS: 6, System.GLONASS: 5}


def _sf(sf: int, value: float) -> int:
    # [spec: SPARTN 6.3] Vencoded = (Voriginal − Rngmin) / Res, after the note's offset
    res, lo, offset = _SF[sf]
    return round((value - offset - lo) / res)


def _description(system: System, t: GnssTime, subtype: int) -> spartn.Description:
    # [spec: SPARTN Table 7.1 TF009] a 32-bit time tag in the system's time, solution 1, processor 0
    tag = timescale.spartn_tag_from_gps(t.week, t.seconds, system, LEAP if system == System.GLONASS else None)
    return spartn.Description(subtype, 1, tag, 1, 0, None, None, None, None)


def _gad(t: GnssTime, siou: int) -> bytes:
    # [spec: SPARTN Table 6.24] one area of 3 × 4 nodes, 1° apart, north-west corner 50° S 164° W
    area = spartn.GadArea(1, _sf(32, -50.0), _sf(33, -164.0), 2, 3, _sf(36, 1.0), _sf(37, 1.0))
    return spartn.encode(spartn.Message(2, 0, 2, _description(System.GPS, t, 0), spartn.Gad(siou, 1, 0, 0, (area,))))


def _hpac(system: System, t: GnssTime, numbers, siou: int, troposphere: bool) -> bytes:
    # [spec: SPARTN Tables 6.14 to 6.23] troposphere (SF043 2.31 m) and per satellite a
    # small ionosphere polynomial (about 10 TECU), no grid
    tropo = None
    if troposphere:
        tropo = spartn.Troposphere(2, 2, _sf(43, 2.31), 0,
                                   spartn.TropoSmall(_sf(45, 0.1), _sf(46, 0.01), _sf(46, -0.02), _sf(47, 0.004)),
                                   None, None, None)
    blocks = tuple(spartn.IonoSatellite(5, 0, spartn.IonoSmall(_sf(57, 10.0 + 0.5 * k), _sf(58, 0.4),
                                                               _sf(58, -0.8), _sf(59, 0.02)), None, None, None)
                   for k, _ in enumerate(numbers))
    mask = sum(1 << (_SAT_BITS[system] - n) for n in numbers)
    area = spartn.HpacArea(1, 0, 1 if tropo else 0, 1, tropo, spartn.Ionosphere(2, 0, mask, blocks))
    return spartn.encode(spartn.Message(1, 0, 2, _description(system, t, int(system)),
                                        spartn.Hpac(siou, 1, 0, 0, (area,))))


def _ocb(system: System, t: GnssTime, blocks: dict, siou: int, end_of_set: bool) -> bytes:
    """[spec: SPARTN 6.5.1, Tables 6.3 to 6.13] an OCB frame; ``blocks`` maps a satellite
    number to None (DNU) or (iode, clock m, [(position, fix, SF015, m)], [(position, m)])."""
    mask, out = 0, []
    for number in sorted(blocks):
        mask |= 1 << (_SAT_BITS[system] - number)
        block = blocks[number]
        if block is None:
            out.append(spartn.OcbSatellite(1, None, None, None, None, None))
            continue
        iode, clock, phase, code = block
        width = _BIAS_BITS[system]
        bias = spartn.Bias(0, sum(1 << (width - 1 - p) for p, *_ in phase),
                           tuple(spartn.PhaseBias(f, k, _sf(20, v)) for _, f, k, v in phase),
                           0, sum(1 << (width - 1 - p) for p, _ in code), tuple(_sf(29, v) for _, v in code))
        out.append(spartn.OcbSatellite(0, 0b111, 7, spartn.Orbit(iode, _sf(20, 0.3), _sf(20, -0.2), _sf(20, 0.1), None),
                                       spartn.Clock(7, _sf(20, clock), 2), bias))
    body = spartn.Ocb(siou, int(end_of_set), 0, 0, 0, 0, 0, mask, tuple(out))
    return spartn.encode(spartn.Message(0, 0, 0, _description(system, t, int(system)), body))


# [design: synthetic biases, metres] GPS phase and code biases at mask positions 0 (1C),
# 1 (2W), 2 (2L), and 3 (5Q) for the satellites that have it; GLONASS 1C and 2C.
_GPS_PHASE = {0: 0.012, 1: -0.034, 2: 0.051, 3: -0.021}
_GPS_CODE = {0: 1.2, 1: -0.4, 2: 0.6, 3: 0.9}
_L5_SATELLITES = (5, 12, 24)


def _gps_block(prn: int, k: int, at: float, s: Scenario) -> tuple | None:
    """[design: the synthetic GPS block of update ``k`` at ``at`` seconds]"""
    if s.dnu is not None and s.dnu[0] <= at < s.dnu[1] and prn == 17:
        return None
    iode = 34 if prn == 5 and s.issue_change is not None and at >= s.issue_change + 20.0 else 33
    clock = 1.0 + 0.1 * prn / 32 + 0.0004 * k            # a clock that changes by 0.4 mm per update
    fixed = 0 if s.float_from is not None and at >= s.float_from and prn == 12 else 1
    positions = (0, 1, 2, 3) if prn in _L5_SATELLITES else (0, 1, 2)
    phase = []
    for p in positions:
        value, continuity = _GPS_PHASE[p], 7
        if prn == 9 and p == 1 and s.bias_arc is not None and at >= s.bias_arc:
            value += 0.1
            continuity = 0 if at < s.bias_arc + 5.0 else 7     # SF015 = 0 in the first block of the new arc
        phase.append((p, fixed, continuity, value))
    return iode, clock, phase, [(p, _GPS_CODE[p]) for p in positions]


def stream(s: Scenario) -> list[tuple[bytes, GnssTime]]:
    """The synthetic input as (bytes, receipt time) pieces, in receipt order.

    [design: DESIGN 5.1, synthetic messages built with gnsswire's encode]
    """
    pieces = [(frame_1019(prn), START - 30.0) for prn in GPS_ORBITS] + [(frame_1020(), START - 29.0)]
    if s.issue_change is not None:
        pieces.append((frame_1019(5, iode=34), START + s.issue_change + 0.05))
    k = 0
    for at in range(-10, s.seconds + 1, 5):
        tag = START + float(at)
        siou = 4 if s.siou_change is not None and at >= s.siou_change else 3
        if at % 30 == 20:                       # GAD and HPAC every 30 s from −10 s
            pieces += [(_gad(tag, siou), tag + 0.1),
                       (_hpac(System.GPS, tag, sorted(GPS_ORBITS), siou, True), tag + 0.2),
                       (_hpac(System.GLONASS, tag, (1,), siou, False), tag + 0.3)]
        glonass = {1: (13, -0.5 + 0.0003 * k, [(0, 1, 7, 0.02), (1, 1, 7, -0.015)], [(0, 0.8), (1, 0.3)])}
        gps = {prn: _gps_block(prn, k, at, s) for prn in GPS_ORBITS}
        pieces += [(_ocb(System.GLONASS, tag, glonass, siou, False), tag + 0.4),
                   (_ocb(System.GPS, tag, gps, siou, True), tag + 0.5)]
        k += 1
    return sorted(pieces, key=lambda p: (p[1].week, p[1].second, p[1].fraction))


def rover(s: Scenario, at: float) -> tuple[float, float, float]:
    """[design: the rover at the site, then 6 km east from ``s.move``] ECEF metres."""
    site = geodesy.geodetic_to_ecef(SITE)
    if s.move is None or at < s.move:
        return site
    d = geodesy.enu_to_ecef((6000.0, 0.0, 0.0), SITE[0], SITE[1])
    return (site[0] + d[0], site[1] + d[1], site[2] + d[2])


@dataclass
class Run:
    """What a run gives: the engine and station after it, the messages of each epoch
    and the whole stream. [design: for the tests and the reader]"""

    engine: Engine
    station: VirtualBase
    epochs: list = field(default_factory=list)      # (t, rover, station, messages)
    wire: bytes = b""


def run(s: Scenario | None = None, on_epoch=None) -> Run:
    """Feed the synthetic input and produce the station's epochs from ``START`` for
    ``s.seconds``; every input is fed before the first epoch at or after its receipt.
    ``on_epoch(t, engine, station, messages)``, when given, is called after each
    epoch, before the next input is fed.

    [design: DESIGN 2.6, a station on the engine, one call per epoch; DESIGN
    2.1, I/O at the edge: bytes and times in, messages out]
    """
    s = s or Scenario()
    engine, framer = Engine(), Framer()
    engine.feed(LeapSeconds(LEAP))
    config = VrsConfig(rate=s.rate, msm=s.msm, float_phase=s.float_phase, **s.extra)
    station = VirtualBase(engine, config)
    pieces = stream(s)
    out = Run(engine, station)
    wire = []
    i = 0
    period = 1.0 / s.rate
    for n in range(int(round(s.seconds * s.rate)) + 1):
        t = START + n * period
        while i < len(pieces) and pieces[i][1] - t <= 0.0:
            piece, received = pieces[i]
            for message in framer.push("ip", piece, received):
                engine.feed(message, received)
            i += 1
        here = rover(s, t - START)
        messages = station.epoch(t, here)
        out.epochs.append((t, here, station.station, messages))
        if on_epoch is not None:
            on_epoch(t, engine, station, messages)
        wire.extend(rtcm3.encode(m) for m in messages)
    out.wire = b"".join(wire)
    return out


if __name__ == "__main__":
    seconds = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 300
    result = run(Scenario(seconds=seconds))
    counts: dict[int, int] = {}
    for _, _, _, messages in result.epochs:
        for m in messages:
            counts[m.number] = counts.get(m.number, 0) + 1
    print(f"{len(result.epochs)} epochs, {len(result.wire)} bytes of RTCM 3")
    for number in sorted(counts):
        print(f"  {number}: {counts[number]} messages")
    if "--write" in sys.argv[:-1]:
        with open(sys.argv[sys.argv.index("--write") + 1], "wb") as f:
            f.write(result.wire)
