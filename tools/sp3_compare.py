"""Corrected and broadcast satellite states of the engine against an SP3 product.

    uv run python tools/sp3_compare.py --eph EPH.bin --ssr SSR.bin --sp3 PRODUCT.sp3[.gz]
        [--sp3 ANOTHER.sp3.gz] [--atx igs20.atx] [--stream PROVIDER:SOLUTION]
        [--leap 18] [--satellites G]

Replays recorded IGS SSR corrections and broadcast ephemerides through the
engine, each recording a gnsswire capture (``name.bin`` and ``name.bin.idx``,
a line "unix_time byte_count" per chunk as it arrived; recordings and
products are never committed), and at each epoch of each SP3 product within
the recording queries every satellite the product gives, with the
corrections received by then. The state at the time of transmission is
moved to the product's epoch by the broadcast velocity and clock drift.

It prints, per product and system, over the satellites the stream corrects:
the orbit differences in radial, along and cross track, corrected and
broadcast; the clock differences after one offset per epoch and system
(the mean over those satellites) is removed, corrected, broadcast, and
with the stream's clock correction negated (the reading rule IG5 rejects);
and radial minus clock, the part of the orbit and clock error a user on
the ground sees. With ``--satellites`` it prints per satellite of those
systems the mean, spread and drift of each over the epochs, the stream's
clock correction and the broadcast issue used. Given two products, it
first prints how well they agree with each other at their common epochs.

SP3 positions are of the centre of mass. The stream's orbits are of an
antenna phase centre (the stream is declared APC, rule IG11), so given an
ANTEX file the product position is moved to the phase centre of the
reference frequency of IGS SSR Table 2 with the nominal attitude; without
one the radial difference holds each satellite's offset.

[design: a tool, not a test, because its inputs are recordings and products
that are not distributed; the statistics are those an orbit and clock
comparison needs, nothing fitted]
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import math
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ssrcore.adapters.framing import Framer
from ssrcore.adapters.igs_ssr import IgsReference, ReferencePoint
from ssrcore.engine import Engine
from ssrcore.model.types import Config, GnssTime, LeapSeconds, MaxAge, Sat, System
from ssrcore.toolkit import astro, geodesy

C = 299792458.0
_GPS_UNIX = 315964800                        # 1980-01-06 00:00:00 UTC in Unix seconds
_GPS0 = datetime.datetime(1980, 1, 6)
# [spec: SP3-d (IGS, Hilla 2016), the system letter of a satellite identifier]
_SYSTEMS = {"G": System.GPS, "R": System.GLONASS, "E": System.GALILEO, "C": System.BEIDOU, "J": System.QZSS}
_LETTERS = {v: k for k, v in _SYSTEMS.items()}
# [spec: IGS SSR v1.00 4.1 Table 2, the reference frequency of an APC stream; ANTEX 1.4
# frequency codes]
_REFERENCE = {System.GPS: "G01", System.GLONASS: "R01", System.GALILEO: "E01", System.BEIDOU: "C02",
              System.QZSS: "J01"}
# [design: a receiver anywhere; the satellite state depends on it only through t_tx]
_RECEIVER = geodesy.geodetic_to_ecef((math.radians(-35.3), math.radians(149.1), 600.0))
# [design: corrections received up to this long after an epoch are those of the epoch
# (the stream's latency); the query is then as a real-time user's]
_LATENCY = 5.0


# ---- inputs --------------------------------------------------------------------------------


def gps_time(unix: float, leap: int) -> GnssTime:
    """[design: as tests/test_igs_ssr_recording.py, GPS time counted from 1980-01-06
    00:00:00 UTC with GPS − UTC = ``leap`` s over the recording] the GPS time of a
    receipt stamp (Unix seconds)."""
    s = unix - _GPS_UNIX + leap
    week = int(s // 604800)
    return GnssTime(week, s - week * 604800)


def chunks(path: Path, transport: str) -> list[tuple[float, str, bytes]]:
    """[design: the gnsswire capture format of the module docstring] (unix time,
    transport, bytes) per arrival; bytes after the last stamp go at its time."""
    data = path.read_bytes()
    out, pos, last = [], 0, 0.0
    for line in Path(str(path) + ".idx").read_text().splitlines():
        if not line.strip():
            continue
        t, n = line.split()
        last = float(t)
        out.append((last, transport, data[pos:pos + int(n)]))
        pos += int(n)
    if pos < len(data):
        out.append((last, transport, data[pos:]))
    return out


def _calendar(fields: list[str]) -> GnssTime:
    """[spec: SP3-d, the epoch header record "*  yyyy mm dd hh mm ss.ssssssss" in the
    file's time system] the epoch as a GPS week and seconds."""
    y, mo, d, h, mi = (int(x) for x in fields[:5])
    s = (datetime.datetime(y, mo, d, h, mi) - _GPS0).total_seconds() + float(fields[5])
    week = int(s // 604800)
    return GnssTime(week, s - week * 604800)


def read_sp3(path: Path, start: GnssTime, end: GnssTime) -> dict[GnssTime, dict[Sat, tuple]]:
    """[spec: SP3-d (IGS, Hilla 2016), the position and clock record: positions in km,
    clocks in µs, 999999.999999 a bad clock, 0.000000 a bad position; the time system
    in columns 10-12 of the first "%c" line] the
    product's (position m, clock s or None) per epoch in [start, end] and satellite.
    Only GPS time is read, so that no leap second enters."""
    opener = gzip.open if path.suffix == ".gz" else open
    out: dict[GnssTime, dict[Sat, tuple]] = {}
    current, scale = None, None
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("%c") and scale is None:
                scale = line[9:12]
                if scale != "GPS":
                    raise ValueError(f"{path}: time system {scale!r}; this tool reads GPS time only")
            elif line.startswith("*"):
                t = _calendar(line[1:].split())
                current = t if start <= t <= end else None
                if current is not None:
                    out[current] = {}
            elif line.startswith("P") and current is not None:
                system = _SYSTEMS.get(line[1])
                if system is None:
                    continue
                xyz = tuple(float(line[4 + 14 * k:18 + 14 * k]) * 1000.0 for k in range(3))
                if xyz == (0.0, 0.0, 0.0):
                    continue
                clk = float(line[46:60])
                out[current][Sat(system, int(line[2:4]))] = (xyz, None if clk > 999999.0 else clk * 1e-6)
    return out


@dataclass(frozen=True)
class Offset:
    """[spec: ANTEX 1.4] A satellite antenna's phase centre offset of one frequency (m,
    satellite body axes x, y, z) and its validity."""

    start: datetime.datetime
    end: datetime.datetime
    xyz: tuple[float, float, float]


def read_atx(path: Path) -> dict[tuple[str, str], list[Offset]]:
    """[spec: ANTEX 1.4 (IGS, Rothacher and Schmid): a satellite antenna block "TYPE /
    SERIAL NO" with the PRN in columns 21-23, "VALID FROM", "VALID UNTIL", and per
    frequency "NORTH / EAST / UP", for a satellite the x, y, z offsets in mm] the
    offsets per (PRN, frequency)."""
    out: dict[tuple[str, str], list[Offset]] = defaultdict(list)
    prn = None
    start = end = None
    freq = None
    for line in path.read_text(errors="replace").splitlines():
        label = line[60:].strip()
        if label == "TYPE / SERIAL NO":
            prn = line[20:23].strip() if line[20:23].strip()[:1] in _SYSTEMS else None
            start, end = datetime.datetime.min, datetime.datetime.max
        elif prn is None:
            continue
        elif label in ("VALID FROM", "VALID UNTIL"):
            p = line[:43].split()
            t = datetime.datetime(*(int(x) for x in p[:5])) + datetime.timedelta(seconds=float(p[5]))
            start, end = (t, end) if label == "VALID FROM" else (start, t)
        elif label == "START OF FREQUENCY":
            freq = line[3:6]
        elif label == "NORTH / EAST / UP" and freq is not None:
            x, y, z = (float(v) * 1e-3 for v in line[:30].split())
            out[(prn, freq)].append(Offset(start, end, (x, y, z)))
            freq = None
        elif label == "END OF ANTENNA":
            prn = None
    return out


def phase_centre(sat: Sat, t: GnssTime, xyz, offsets, leap: int):
    """[spec: the IGS axes the ANTEX file's header names (Montenbruck et al. 2015, GNSS
    satellite geometry and attitude models): z towards the Earth's centre, y along
    z × (Sun − satellite), x completing the right-handed set, under nominal yaw
    steering; IGS SSR 4.1, the offset of the reference frequency of Table 2. design:
    the nominal attitude, as no attitude is an input; x and y offsets lie across the
    radial, so the radial and the clock comparisons do not depend on it] the antenna
    phase centre of a centre-of-mass position, or None when the file has no offset
    for the satellite at ``t``."""
    when = _GPS0 + datetime.timedelta(weeks=t.week, seconds=t.seconds)
    found = [o for o in offsets.get((f"{_LETTERS[sat.system]}{sat.number:02d}", _REFERENCE[sat.system]), ())
             if o.start <= when < o.end]
    if not found:
        return None
    sun = astro.sky(t.week, t.seconds, -float(leap)).sun
    r = math.sqrt(sum(v * v for v in xyz))
    ez = tuple(-v / r for v in xyz)
    es = tuple(s - v for s, v in zip(sun, xyz))
    ey = _unit(_cross(ez, es))
    ex = _cross(ey, ez)
    dx, dy, dz = found[0].xyz
    return tuple(p + dx * a + dy * b + dz * c for p, a, b, c in zip(xyz, ex, ey, ez))


def _cross(a, b):
    """[design: vector helper]"""
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a):
    """[design: vector helper]"""
    n = math.sqrt(sum(v * v for v in a))
    return tuple(v / n for v in a)


def rac(dx, position, velocity) -> tuple[float, float, float]:
    """[spec: IGS SSR v1.00 4.1, e_along = ṙ/|ṙ|, e_cross = r × ṙ/|r × ṙ|, e_radial =
    e_along × e_cross] the radial, along and cross components of ``dx``."""
    along = _unit(velocity)
    cross = _unit(_cross(position, velocity))
    radial = _cross(along, cross)
    return tuple(sum(p * q for p, q in zip(dx, e)) for e in (radial, along, cross))


# ---- the replay ----------------------------------------------------------------------------


@dataclass(frozen=True)
class State:
    """[design: what the comparison needs of a record] The engine's satellite at a
    product epoch: corrected and broadcast position (m)
    and clock (s), velocity, the stream's clock correction δC (m, the engine's sign,
    t_s = t_b − δC/c) and orbit correction (radial, along, cross m), and the broadcast
    issue used."""

    position: tuple
    clock: float
    broadcast_position: tuple
    broadcast_clock: float
    velocity: tuple
    clock_correction: float
    orbit_correction: tuple
    iod: int


def replay(eph: list[Path], ssr: list[Path], epochs: list[GnssTime], sats: dict[GnssTime, set],
           stream: tuple[int, int], leap: int) -> dict[tuple[GnssTime, Sat], State]:
    """[spec: DESIGN 2.2, 2.3; design: inputs in arrival order through the framer, as
    tests/test_igs_ssr_recording.py feeds them; an epoch is queried once every input
    received within ``_LATENCY`` of it has been fed; ages as the stream's update
    intervals need (IG2)] the corrected states at each epoch."""
    pieces = sorted([c for p in eph for c in chunks(p, "eph")] + [c for p in ssr for c in chunks(p, "ssr")],
                    key=lambda c: c[0])
    engine = Engine(Config(max_age=MaxAge(orbit=120.0, clock=30.0, bias=300.0)))
    engine.feed(LeapSeconds(leap))
    engine.feed(IgsReference(stream[0], stream[1], ReferencePoint.APC))
    framer = Framer()
    out: dict[tuple[GnssTime, Sat], State] = {}
    k = 0

    def query(t: GnssTime) -> None:
        """[design: the state at t_tx moved to t by the broadcast velocity and drift]"""
        for sat in sats[t]:
            rec = engine.query(sat, t + 0.075, _RECEIVER)
            if not rec.available or not rec.corrected:
                continue
            s = rec.satellite
            dt = t - rec.t_tx
            out[(t, sat)] = State(
                tuple(p + v * dt for p, v in zip(s.position, s.velocity)), s.clock + s.clock_drift * dt,
                tuple(p + v * dt for p, v in zip(s.broadcast.position, s.velocity)),
                s.broadcast.clock + s.clock_drift * dt, s.velocity, rec.clock_correction.correction,
                (rec.orbit.radial, rec.orbit.along, rec.orbit.cross), s.broadcast.iod)

    for unix, transport, data in pieces:
        received = gps_time(unix, leap)
        while k < len(epochs) and received - epochs[k] > _LATENCY:
            query(epochs[k])
            k += 1
        for m in framer.push(transport, data, received):
            engine.feed(m, received)
    return out


# ---- statistics ----------------------------------------------------------------------------


def _rms(values) -> float:
    """[design: root mean square; nan for none]"""
    values = list(values)
    return math.sqrt(sum(v * v for v in values) / len(values)) if values else float("nan")


@dataclass
class Row:
    """[design: one comparison] One satellite at one epoch against the product (m):
    orbit (radial, along, cross) corrected and broadcast, clock corrected and broadcast
    before the epoch's offset is removed, and the stream's corrections."""

    t: GnssTime
    sat: Sat
    orbit: tuple
    broadcast_orbit: tuple
    clock: float
    broadcast_clock: float
    state: State


def compare(product, states, offsets, leap: int) -> list[Row]:
    """[spec: IGS SSR v1.00 4.1 (frame), 4.2 (clock of satellite times, compared as
    offsets); design: the product's clock is of the centre of mass and the APC alike,
    its offset is not moved] the rows of the satellites the stream corrects."""
    rows = []
    for t, sats in sorted(product.items()):
        for sat, (xyz, clock) in sats.items():
            st = states.get((t, sat))
            if st is None or clock is None:
                continue
            ref = xyz
            if offsets is not None:
                ref = phase_centre(sat, t, xyz, offsets, leap)
                if ref is None:
                    continue
            rows.append(Row(t, sat, rac([p - q for p, q in zip(st.position, ref)], ref, st.velocity),
                            rac([p - q for p, q in zip(st.broadcast_position, ref)], ref, st.velocity),
                            (st.clock - clock) * C, (st.broadcast_clock - clock) * C, st))
    return rows


def offsets_removed(rows: list[Row]) -> dict[int, dict[str, float]]:
    """[design: one offset per epoch and system, the mean over the satellites the stream
    corrects, removed from the clock and from radial minus clock: a product's and a
    stream's clocks have each their own reference] per row index, the clock and radial
    minus clock residuals, corrected and broadcast."""
    groups = defaultdict(list)
    for i, r in enumerate(rows):
        groups[(r.t, r.sat.system)].append(i)
    out: dict[int, dict[str, float]] = {}
    for idx in groups.values():
        series = {"clock": [rows[i].clock for i in idx],
                  "broadcast clock": [rows[i].broadcast_clock for i in idx],
                  "clock, δC negated": [2.0 * rows[i].broadcast_clock - rows[i].clock for i in idx],
                  "radial - clock": [rows[i].orbit[0] - rows[i].clock for i in idx],
                  "broadcast radial - clock": [rows[i].broadcast_orbit[0] - rows[i].broadcast_clock for i in idx]}
        means = {k: statistics.fmean(v) for k, v in series.items()}
        for j, i in enumerate(idx):
            out[i] = {k: v[j] - means[k] for k, v in series.items()}
    return out


def report(name: str, rows: list[Row], per_satellite: str) -> None:
    """[design: the summary of the module docstring]"""
    print(f"\n== {name}: {len({r.t for r in rows})} epochs, {len(rows)} satellite epochs corrected")
    res = offsets_removed(rows)
    for system in sorted({r.sat.system for r in rows}, key=lambda s: s.value):
        idx = [i for i, r in enumerate(rows) if r.sat.system == system]
        sats = len({rows[i].sat for i in idx})
        for label, get in (("corrected", lambda r: r.orbit), ("broadcast", lambda r: r.broadcast_orbit)):
            comp = [_rms(get(rows[i])[k] for i in idx) for k in range(3)]
            mean_r = statistics.fmean(get(rows[i])[0] for i in idx)
            print(f"  {system.name:8s} {label:9s} n {len(idx):4d} ({sats} satellites)  orbit RMS radial {comp[0]:.3f} "
                  f"along {comp[1]:.3f} cross {comp[2]:.3f} m (mean radial {mean_r:+.3f})")
        for key in ("clock", "broadcast clock", "clock, δC negated", "radial - clock", "broadcast radial - clock"):
            v = [res[i][key] for i in idx]
            print(f"  {system.name:8s} {key:24s} offset per epoch removed: RMS {_rms(v):.3f} m, "
                  f"largest {max(abs(x) for x in v):.3f} m")
    if not per_satellite:
        return
    print("  per satellite over the epochs (m; drift in m per 10 min): clock corrected mean, spread, drift | "
          "broadcast mean | radial - clock corrected, broadcast | radial corrected | stream δC, radial; IOD")
    by_sat = defaultdict(list)
    for i, r in enumerate(rows):
        if _LETTERS[r.sat.system] in per_satellite:
            by_sat[r.sat].append(i)
    for sat, idx in sorted(by_sat.items(), key=lambda x: (x[0].system.value, x[0].number)):
        c = [res[i]["clock"] for i in idx]
        ts = [rows[i].t - rows[idx[0]].t for i in idx]
        drift = float("nan")
        if len(idx) > 1 and max(ts) > 0:
            mt, mc = statistics.fmean(ts), statistics.fmean(c)
            drift = 600.0 * sum((a - mt) * (b - mc) for a, b in zip(ts, c)) / sum((a - mt) ** 2 for a in ts)
        st = rows[idx[-1]].state
        print(f"    {_LETTERS[sat.system]}{sat.number:02d} n {len(idx)}  {statistics.fmean(c):+.3f} "
              f"{statistics.pstdev(c):.3f} {drift:+.3f} | {statistics.fmean(res[i]['broadcast clock'] for i in idx):+.3f}"
              f" | {statistics.fmean(res[i]['radial - clock'] for i in idx):+.3f} "
              f"{statistics.fmean(res[i]['broadcast radial - clock'] for i in idx):+.3f}"
              f" | {statistics.fmean(rows[i].orbit[0] for i in idx):+.3f}"
              f" | {st.clock_correction:+.3f} {st.orbit_correction[0]:+.3f}; {st.iod}")


def products_agree(products: dict[str, dict]) -> None:
    """[design: two products are a check on each other before either is a check on the
    stream; one clock offset per epoch and system removed, as in ``offsets_removed``]
    per pair of products and system, the position and clock differences at their
    common epochs and satellites."""
    names = list(products)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pos, clk = defaultdict(list), defaultdict(list)
            for t in sorted(set(products[a]) & set(products[b])):
                groups = defaultdict(list)
                for sat, (xyz, clock) in products[a][t].items():
                    other = products[b][t].get(sat)
                    if other is None:
                        continue
                    pos[sat.system].append(math.dist(xyz, other[0]))
                    if clock is not None and other[1] is not None:
                        groups[sat.system].append((clock - other[1]) * C)
                for system, v in groups.items():
                    m = statistics.fmean(v)
                    clk[system] += [x - m for x in v]
            for system in sorted(pos, key=lambda s: s.value):
                c = clk.get(system, [])
                print(f"{a} against {b}, {system.name}: position 3D RMS {_rms(pos[system]):.3f} m (n {len(pos[system])}); "
                      f"clock RMS {_rms(c):.3f} m, largest {max((abs(x) for x in c), default=float('nan')):.3f} m "
                      f"(n {len(c)})")


def main(argv: list[str]) -> int:
    """[design: command line of the module docstring]"""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--eph", type=Path, action="append", required=True, help="broadcast ephemeris capture")
    p.add_argument("--ssr", type=Path, action="append", required=True, help="IGS SSR capture (APC)")
    p.add_argument("--sp3", type=Path, action="append", required=True, help="SP3 product (GPS time)")
    p.add_argument("--atx", type=Path, help="ANTEX file of the products' satellite antennas")
    p.add_argument("--stream", default="0:0", help="the stream declared APC, PROVIDER:SOLUTION")
    p.add_argument("--leap", type=int, default=18, help="GPS − UTC (s) over the recording")
    p.add_argument("--satellites", default="", help="system letters to list per satellite, e.g. GE")
    a = p.parse_args(argv[1:])
    first = min(min(c[0] for c in chunks(f, "")) for f in a.ssr + a.eph)
    last = max(max(c[0] for c in chunks(f, "")) for f in a.ssr + a.eph)
    start, end = gps_time(first, a.leap) + 60.0, gps_time(last, a.leap) - _LATENCY
    products = {f.name: read_sp3(f, start, end) for f in a.sp3}
    epochs = sorted({t for prod in products.values() for t in prod})
    if not epochs:
        print(f"no product epoch within the recording ({start} to {end})")
        return 1
    sats: dict[GnssTime, set] = defaultdict(set)
    for prod in products.values():
        for t, s in prod.items():
            sats[t] |= set(s)
    provider, solution = (int(x) for x in a.stream.split(":"))
    states = replay(a.eph, a.ssr, epochs, sats, (provider, solution), a.leap)
    offsets = read_atx(a.atx) if a.atx else None
    print(f"recording {start} to {end}; satellite reference point of the product: "
          f"{'APC of IGS SSR Table 2 from ' + a.atx.name if offsets else 'centre of mass (no ANTEX given)'}")
    products_agree(products)
    for name, prod in products.items():
        report(name, compare(prod, states, offsets, a.leap), a.satellites)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
