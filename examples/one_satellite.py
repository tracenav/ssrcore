"""Feed a short synthetic stream to the engine and print one satellite's record.

The stream is built here with gnsswire from chosen values: an RTCM 1019
ephemeris of GPS 5 (the public broadcast record of 2024-01-01 18:00 GPS
time), then one SPARTN set: a GAD area, an HPAC message and an OCB message
with orbit, clock and biases. Its bytes go through ``Framer`` as a receiver
would deliver them, each piece with its receipt time, and the record of G05
for a receiver in that area is printed term by term.

Run it with ``uv run python examples/one_satellite.py``.
"""

from __future__ import annotations

import dataclasses
import math

from gnsswire import rtcm3, spartn

from ssrcore import Engine, Framer, GnssTime, LeapSeconds, LineOfSight, Sat, System
from ssrcore.toolkit import geodesy, timescale

T = GnssTime(2295, 87400.0)        # 2024-01-01 00:16:22 GPS time, the epoch of the query
SPARTN_TIME = T - 2.0              # the SPARTN set's time tag
RECEIVER = geodesy.geodetic_to_ecef((math.radians(-51.4), math.radians(-162.7), 100.0))
PI_ICD = 3.1415926535898           # IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2


def _q(value: float, exponent: int, semicircle: bool = False) -> int:
    # [spec: RTCM Table 3.5-21] the integer of a value at resolution 2^exponent
    return round((value / PI_ICD if semicircle else value) / 2.0 ** exponent)


def _sf(value: float, resolution: float, range_min: float) -> int:
    # [spec: SPARTN 6.3] Vencoded = (Voriginal − Rngmin) / Res
    return round((value - range_min) / resolution)


def _rtcm_1019() -> bytes:
    # [spec: RTCM Table 3.5-21] the broadcast elements of G05, IGS BRDC 2024-01-01
    fields = {f.name: 0 for f in dataclasses.fields(rtcm3.Msg1019) if f.name != "extension"}
    fields.update(
        number=1019, satellite=5, week=2295 % 1024, iode=33, iodc=33, toe=93600 // 16, toc=93600 // 16,
        af0=_q(-1.562642864883e-04, -31), af1=_q(-1.477928890381e-12, -43),
        sqrt_a=_q(5153.783403397, -19), e=_q(0.005790898343548, -33), m0=_q(-2.494229256397, -31, True),
        delta_n=_q(4.453756945425e-09, -43, True), omega=_q(1.235180510432, -31, True),
        omega0=_q(-0.6299845609285, -31, True), omega_dot=_q(-8.155696860368e-09, -43, True),
        i0=_q(0.9693533956084, -31, True), idot=_q(-3.928735076157e-11, -43, True),
        cuc=_q(-1.747161149979e-06, -29), cus=_q(3.31737101078e-06, -29), cic=_q(1.098960638046e-07, -29),
        cis=_q(4.097819328308e-08, -29), crc=_q(319.34375, -5), crs=_q(-34.125, -5))
    return rtcm3.encode(rtcm3.Msg1019(**fields))


def _description(subtype: int) -> spartn.Description:
    # [spec: SPARTN Table 7.1 TF009] a 32-bit time tag, solution 1, processor 0
    tag = timescale.spartn_tag_from_gps(SPARTN_TIME.week, SPARTN_TIME.seconds, System.GPS)
    return spartn.Description(subtype, 1, tag, 1, 0, None, None, None, None)


def _gad() -> bytes:
    # [spec: SPARTN Table 6.24] one area of 3 × 4 nodes, 1° apart, north-west corner 50° S 164° W
    area = spartn.GadArea(1, _sf(-50.0, 0.1, -90.0), _sf(-164.0, 0.1, -180.0), 2, 3,
                          _sf(1.0, 0.1, 0.1), _sf(1.0, 0.1, 0.1))
    return spartn.encode(spartn.Message(2, 0, 2, _description(0), spartn.Gad(3, 1, 0, 0, (area,))))


def _hpac() -> bytes:
    # [spec: SPARTN Tables 6.14 to 6.23] troposphere (SF043 2.31 m, small coefficients) and
    # the ionosphere of G05 (SF057 10 TECU, small coefficients), no grid
    tropo = spartn.Troposphere(2, 2, _sf(2.31 - 2.3, 0.004, -0.508), 0,
                               spartn.TropoSmall(_sf(0.1 - 0.252, 0.004, -0.252), _sf(0.01, 0.001, -0.063),
                                                 _sf(-0.02, 0.001, -0.063), _sf(0.004, 0.0002, -0.051)),
                               None, None, None)
    iono = spartn.Ionosphere(2, 0, 1 << (32 - 5), (spartn.IonoSatellite(
        5, 0, spartn.IonoSmall(_sf(10.0, 0.04, -81.88), _sf(0.4, 0.008, -16.376), _sf(-0.8, 0.008, -16.376),
                               _sf(0.02, 0.002, -8.190)), None, None, None),))
    area = spartn.HpacArea(1, 0, 1, 1, tropo, iono)
    return spartn.encode(spartn.Message(1, 0, 2, _description(0), spartn.Hpac(3, 1, 0, 0, (area,))))


def _ocb() -> bytes:
    # [spec: SPARTN Tables 6.3 to 6.13] G05: orbit (IODE 33; 0.5, −0.25, 0.124 m), clock 1.234 m,
    # L1C and L2W phase biases (L1C fixed), L1C and L2W code biases
    m = lambda v: _sf(v, 0.002, -16.382)   # noqa: E731, SF020
    bias = spartn.Bias(0, 0b110000, (spartn.PhaseBias(1, 7, m(0.01)), spartn.PhaseBias(0, 7, m(-0.03))),
                       0, 0b110000, (_sf(1.2, 0.02, -20.46), _sf(-0.4, 0.02, -20.46)))
    block = spartn.OcbSatellite(0, 0b111, 7, spartn.Orbit(33, m(0.5), m(-0.25), m(0.124), None),
                                spartn.Clock(7, m(1.234), 2), bias)
    body = spartn.Ocb(3, 1, 0, 0, 0, 0, 0, 1 << (32 - 5), (block,))
    return spartn.encode(spartn.Message(0, 0, 0, _description(0), body))


def stream() -> list[tuple[bytes, GnssTime]]:
    """The synthetic stream as (bytes, receipt time) pieces.

    [design: the smallest stream that gives a complete record: an ephemeris,
    then a SPARTN set whose OCB message ends it (SF010)]
    """
    return [(_rtcm_1019(), T - 30.0), (_gad(), SPARTN_TIME + 0.1), (_hpac(), SPARTN_TIME + 0.2),
            (_ocb(), SPARTN_TIME + 0.5)]


def run() -> tuple[Engine, LineOfSight]:
    """Feed ``stream()`` through a framer and query G05 at ``T``.

    [design: DESIGN 2.1, 2.2: framed bytes and receipt times in, one record out]
    """
    engine, framer = Engine(), Framer()
    # GPS − UTC as a receiver's NAV-TIMELS gives it; the tide needs it (UT1, DESIGN 2.4)
    engine.feed(LeapSeconds(18))
    for piece, received in stream():
        for message in framer.push("ip", piece, received):
            engine.feed(message, received)
    return engine, engine.query(Sat.gps(5), T, RECEIVER)


def _m(value: float | None, unit: str = "m") -> str:
    # [design: an absent term prints as absent; its reason is listed under "missing"]
    return "absent" if value is None else f"{value:.4f} {unit}"


def describe(los: LineOfSight) -> list[str]:
    """The record as lines of text, one term per line.

    [design: DESIGN 2.4, each term with its unit; absent terms with their reason]
    """
    lines = [f"{los.sat} received {los.t_rx}, corrected={los.corrected}"]
    if not los.available:
        return lines + [f"  missing: {term.value} ({reason.value})" for (term, _), reason in los.missing.items()]
    g, s, m = los.geometry, los.satellite, los.modelled
    lines += [
        f"  transmit time      {los.t_tx}",
        f"  elevation/azimuth  {math.degrees(g.elevation):.3f} / {math.degrees(g.azimuth):.3f} deg",
        f"  range, rate        {g.range:.4f} m, {g.range_rate:.4f} m/s",
        f"  position           {s.position[0]:.4f} {s.position[1]:.4f} {s.position[2]:.4f} m",
        f"  clock              {s.clock * 1e6:.6f} us (drift {s.clock_drift:.3e})",
        f"  orbit correction   {los.orbit.radial:.3f} {los.orbit.along:.3f} {los.orbit.cross:.3f} m RAC,"
        f" age {los.orbit.age:.1f} s",
        f"  clock correction   {los.clock_correction.correction:.3f} m, age {los.clock_correction.age:.1f} s",
        f"  Sagnac, Shapiro    {m.earth_rotation:.4f} m, {m.shapiro:.4f} m",
        f"  dt_rel             {m.relativity_clock * 1e9:.4f} ns",
        f"  solid tide         {_m(m.solid_tide)}",
        f"  wind-up            {_m(m.windup, 'cycle')}",
        f"  troposphere        {los.troposphere.slant:.4f} m (hydrostatic zenith "
        f"{los.troposphere.hydrostatic_zenith:.4f} m)",
    ]
    for signal in los.signals:
        lines.append(f"  {signal.code:4s}  code bias {los.code_bias[signal].value:+.3f} m, "
                     f"phase bias {los.phase_bias[signal].value:+.3f} m"
                     f" ({'fixed' if los.phase_bias[signal].fixed else 'float'}), "
                     f"ionosphere {los.ionosphere.code_delay[signal]:.4f} m")
    lines += [f"  missing: {term.value} ({reason.value})" for (term, _), reason in los.missing.items()]
    return lines


if __name__ == "__main__":
    _, record = run()
    print("\n".join(describe(record)))
