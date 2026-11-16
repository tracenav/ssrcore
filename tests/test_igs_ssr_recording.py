"""The IGS SSR adapter on real recordings (local data, never distributed).

``../gnsswire/build/private/live`` holds recordings of public IGS real-time
streams, each ``name.bin`` with ``name.bin.idx``, a line "unix_time
byte_count" per chunk as it arrived. The tests feed them through the framer
and the engine in arrival order and check properties of the result; nothing
of the data is copied here, and each test skips with a named reason when
its recording is absent.

The recording ``igs_*_ssr_cas_apc.bin`` is an antenna-phase-centre stream
with the broadcast ephemerides of the same time in ``igs_*_eph_bcep.bin``;
``ga_SSRA01CAS1.bin`` is an APC stream without ephemerides and
``ga_SSRC01CAS1.bin`` a centre-of-mass stream that also carries the VTEC
(IM201). [design: recordings are not committed (DESIGN 5.1, Public tests),
so a test on local data skips when it is absent; the bounds are stated per
test]
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from gnsswire import rtcm3

from ssrcore.adapters import igs_ssr
from ssrcore.adapters.framing import Framer
from ssrcore.adapters.igs_ssr import IgsReference, IgsSsr, IgsStreamId, ReferencePoint
from ssrcore.engine import Engine
from ssrcore.model.ephemeris import EphemerisStore
from ssrcore.model.igs_vtec import VtecInput
from ssrcore.model.types import (
    Config,
    FeedOutcome,
    FeedReason,
    FrameFailure,
    GnssTime,
    LeapSeconds,
    MaxAge,
    Missing,
    System,
)
from ssrcore.toolkit import geodesy
from ssrcore.toolkit import ionosphere as iono
from ssrcore.toolkit import orbits

LIVE = Path(__file__).resolve().parents[2] / "gnsswire" / "build" / "private" / "live"
C = 299792458.0
# [design: GPS - UTC over the recordings (2026), no leap second since 2017]
LEAP = 18
_GPS_UNIX = 315964800          # 1980-01-06 00:00:00 UTC in Unix seconds
# [design: a receiver anywhere; the checks are of the satellite states]
RECEIVER = geodesy.geodetic_to_ecef((math.radians(-35.3), math.radians(149.1), 600.0))


def _recording(pattern: str) -> Path:
    found = sorted(LIVE.glob(pattern)) if LIVE.is_dir() else []
    if not found or not Path(str(found[0]) + ".idx").is_file():
        pytest.skip(f"no recording {LIVE / pattern} with its .idx (local data, never distributed)")
    return found[0]


def _gps(unix: float) -> GnssTime:
    s = unix - _GPS_UNIX + LEAP
    week = int(s // 604800)
    return GnssTime(week, s - week * 604800)


def _chunks(path: Path, transport: str) -> list[tuple[float, str, bytes]]:
    """(unix time, transport, bytes) per arrival; bytes after the last stamp go at its time."""
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


def _messages(*files: tuple[Path, str]):
    """Every framed message of the files in arrival order (stable for equal stamps), with its time."""
    pieces = sorted((c for path, transport in files for c in _chunks(path, transport)), key=lambda c: c[0])
    framer = Framer()
    for unix, transport, data in pieces:
        received = _gps(unix)
        for m in framer.push(transport, data, received):
            yield m, received


# ---- the APC stream with its ephemerides ---------------------------------------------------


@dataclass
class Replay:
    outcomes: Counter = field(default_factory=Counter)
    found: dict = field(default_factory=dict)            # (sat, IGS IOD) -> [found at each orbit message]
    steps: dict = field(default_factory=lambda: defaultdict(list))   # (system, issue changed) -> (clock m, position m, sat, t)
    largest: dict = field(default_factory=lambda: defaultdict(lambda: [0.0] * 4))  # system -> radial, along, cross, C0
    rejected: list = field(default_factory=list)
    engine: Engine | None = None


@pytest.fixture(scope="module")
def replay() -> Replay:
    ssr = _recording("igs_*_ssr_cas_apc.bin")
    eph = _recording("igs_*_eph_bcep.bin")
    r = Replay()
    e = r.engine = Engine(Config(max_age=MaxAge(orbit=120.0, clock=30.0, bias=300.0)))
    e.feed(LeapSeconds(LEAP))
    declared: set[IgsStreamId] = set()
    for m, received in _messages((ssr, "ssr"), (eph, "eph")):
        satellites = ()
        if isinstance(m, (rtcm3.IgsOrbit, rtcm3.IgsClock)):
            stream = IgsStreamId.of(m)
            if stream not in declared:      # the file is an APC stream (its name, "apc")
                e.feed(IgsReference(m.provider, m.solution, ReferencePoint.APC))
                declared.add(stream)
            system = igs_ssr.system_of(m)
            satellites = [(s, igs_ssr.satellite(system, s.satellite)) for s in m.satellites]
            satellites = [(s, sat) for s, sat in satellites if sat is not None]
        before = {sat: e.query(sat, received, RECEIVER) for _, sat in satellites}
        result = e.feed(m, received)
        key = (type(m).__name__, result.outcome.value, result.reason.value if result.reason else None)
        r.outcomes[key] += 1
        if result.outcome == FeedOutcome.REJECTED and type(m).__name__.startswith("Igs"):
            r.rejected.append((key, result.detail))
        if not satellites:
            continue
        epoch = igs_ssr.epoch_time(m, system, received)
        for s, sat in satellites:
            if isinstance(m, rtcm3.IgsOrbit) and result.outcome != FeedOutcome.DUPLICATE \
                    and not (sat.system == System.BEIDOU and orbits.beidou_geo(sat.number)):
                iod = igs_ssr.ephemeris_iod(sat, s.iod, e.ephemerides)
                hit = None if iod is None else \
                    e.ephemerides.find(sat, iod, igs_ssr.ephemeris_type(sat), epoch).stored
                r.found.setdefault((sat, s.iod), []).append(hit is not None)
                big = r.largest[sat.system.name]
                for i, v in enumerate((s.radial * 1e-4, s.along * 4e-4, s.cross * 4e-4)):
                    big[i] = max(big[i], abs(v))
            if isinstance(m, rtcm3.IgsClock):
                big = r.largest[sat.system.name]
                big[3] = max(big[3], abs(s.c0 * 1e-4))
            old, new = before[sat], e.query(sat, received, RECEIVER)
            if old.corrected and new.corrected:
                changed = old.orbit.iode != new.orbit.iode
                r.steps[(sat.system.name, changed)].append(
                    (abs(new.satellite.clock - old.satellite.clock) * C,
                     math.dist(new.satellite.position, old.satellite.position), sat, received))
    return r


def test_recording_every_igs_message_is_taken(replay: Replay):
    """Every 4076 message of the stream decodes and is stored or a repeat; nothing is refused."""
    igs = {k: v for k, v in replay.outcomes.items() if k[0].startswith("Igs")}
    assert igs and not replay.rejected, replay.rejected[:5]
    assert all(k[1] in (FeedOutcome.STORED.value, FeedOutcome.DUPLICATE.value) for k in igs), igs
    failures = [k for k in replay.outcomes if k[0] == FrameFailure.__name__]
    print("\noutcomes:", dict(sorted(replay.outcomes.items(), key=str)))
    print("left out by the adapter:", dict(replay.engine.igs.left_out))
    assert all(k[1] == FeedOutcome.COUNTED.value for k in failures)


def test_recording_orbit_corrections_find_their_ephemeris_once_received(replay: Replay):
    """Once an issue is in the store, every later correction naming it finds it; the issues
    never found are named (broadcast before the recording began, so never received)."""
    never = []
    for (sat, iod), hits in replay.found.items():
        if any(hits):
            first = hits.index(True)
            assert all(hits[first:]), f"{sat} IOD {iod} lost after it was found"
        else:
            never.append((sat, iod, len(hits)))
    total = sum(len(h) for h in replay.found.values())
    unmatched = sum(n for *_, n in never)
    print(f"\norbit corrections {total}, never matched {unmatched}:",
          ", ".join(f"{s.system.name[0]}{s.number:02d} IOD {i} ({n})" for s, i, n in sorted(never, key=str)))
    # [design: a stream may name an issue the recording began after; most must match]
    assert unmatched <= 0.1 * total
    # never received: at the end the store still holds no ephemeris of that issue
    e = replay.engine
    end = max(t for v in replay.steps.values() for *_, t in v)
    for sat, iod, _ in never:
        mapped = igs_ssr.ephemeris_iod(sat, iod, e.ephemerides)
        if mapped is not None:
            miss = e.ephemerides.find(sat, mapped, igs_ssr.ephemeris_type(sat), end).miss
            assert miss == Missing.IODE_MISMATCH, (sat, iod, miss)


def test_recording_the_corrected_clock_is_continuous_across_updates(replay: Replay):
    """At each orbit or clock update, the same query before and after: the corrected clock
    moves by the provider's change over one interval, and across a broadcast issue change
    the correction follows the broadcast clock (IG5 sign, IG4 pairing).
    [design: bounds well above the C0 step of 0.1 mm and the provider's few-cm changes
    per interval, and far below a wrong sign or pairing (0.2 m and more)]"""
    for (system, changed), v in sorted(replay.steps.items()):
        d = sorted(x[0] for x in v)
        worst = max(v, key=lambda x: x[0])
        print(f"\n{system} issue change {changed}: {len(v)} updates, clock step p50 {d[len(d) // 2]:.4f} "
              f"p99 {d[int(len(d) * 0.99)]:.4f} max {d[-1]:.4f} m ({worst[2]} at {worst[3].second:.0f})")
        assert d[-1] < (0.02 if changed else 0.10)
        assert d[int(len(d) * 0.99)] < 0.02


def test_recording_corrections_are_small_against_the_broadcast(replay: Replay):
    """[design: a correction of the broadcast is metres, the broadcast itself 2e7 m and up
    to 1e-3 s (3e5 m) of clock: orbit corrections under 10 m and clocks under 100 m]"""
    for system, (radial, along, cross, c0) in sorted(replay.largest.items()):
        print(f"\nlargest {system}: radial {radial:.3f} along {along:.3f} cross {cross:.3f} C0 {c0:.3f} m")
        assert max(radial, along, cross) < 10.0 and c0 < 100.0


def test_recording_the_corrected_position_changes_smoothly(replay: Replay):
    """[design: an update moves the corrected position by the provider's change of orbit
    over one interval, millimetres; 1 cm bounds it]"""
    for (system, changed), v in sorted(replay.steps.items()):
        print(f"\n{system} issue change {changed}: position step max {max(x[1] for x in v):.4f} m")
        assert max(x[1] for x in v) < 0.01


# ---- the streams without ephemerides -------------------------------------------------------


def test_recording_an_apc_stream_without_ephemerides_is_stored():
    """Every 4076 message decodes; the satellite corrections are stored, and a query has no
    ephemeris to apply them to."""
    path = _recording("ga_SSRA01CAS1.bin")
    e = Engine()
    e.feed(LeapSeconds(LEAP))
    outcomes: Counter = Counter()
    streams, sats, last = set(), set(), None
    for m, received in _messages((path, "ssr")):
        if isinstance(m, igs_ssr.IGS_SATELLITE_MESSAGES):
            stream = IgsStreamId.of(m)
            if stream not in streams:
                e.feed(IgsReference(m.provider, m.solution, ReferencePoint.APC))
                streams.add(stream)
            system = igs_ssr.system_of(m)
            sats |= {igs_ssr.satellite(system, s.satellite) for s in m.satellites} - {None}
        r = e.feed(m, received)
        outcomes[(type(m).__name__, r.outcome.value, r.reason.value if r.reason else None)] += 1
        last = received
    print("\noutcomes:", dict(sorted(outcomes.items(), key=str)))
    igs = {k: v for k, v in outcomes.items() if k[0].startswith("Igs")}
    assert igs
    assert all(k[1] in (FeedOutcome.STORED.value, FeedOutcome.DUPLICATE.value)
               or k[2] == FeedReason.UNSUPPORTED_SYSTEM.value for k in igs), igs
    assert sats and last is not None
    missing = {e.query(sat, last, RECEIVER).available for sat in sats}
    assert missing == {False}


def test_recording_a_com_stream_is_refused_and_its_vtec_is_plausible():
    """The satellite messages are refused (IG11, centre of mass); the IM201 messages adapt,
    every coefficient present, and the VTEC they give is that of the Earth's ionosphere.
    [design: a wrong normalization of the Legendre functions scales the terms of degree n
    by factors up to sqrt(2(2n+1)) and more, a wrong sign of the sun-fixed longitude moves
    the daytime maximum to night; a VTEC between 0 and 200 TECU everywhere, with its
    maximum on the day side, is what the Earth's ionosphere gives]"""
    path = _recording("ga_SSRC01CAS1.bin")
    adapter = IgsSsr()
    refused: Counter = Counter()
    vtecs: list[VtecInput] = []
    failures: Counter = Counter()
    for m, received in _messages((path, "ssr")):
        if isinstance(m, FrameFailure):
            failures[(m.message_type, m.subtype)] += 1
            continue
        if isinstance(m, igs_ssr.IGS_SATELLITE_MESSAGES):
            adapter.declare(IgsReference(m.provider, m.solution, ReferencePoint.COM))
        a = adapter.input(m, received, EphemerisStore()) if igs_ssr.is_igs(m) else None
        if isinstance(a, igs_ssr.Refused):
            refused[a.reason.split(" gives ")[-1][:40]] += 1
        elif isinstance(a, igs_ssr.Adapted) and isinstance(a.item, VtecInput):
            vtecs.append(a.item)
    print("\nrefused:", dict(refused), "IM201:", len(vtecs), "frame failures:", dict(failures))
    assert refused and all("centre-of-mass" in k or "SBAS" in k for k in refused)
    # [spec: IGS SSR v1.00 Table 5, the sub-types; a 4076 frame of another sub-type is counted]
    defined = {base + k for base in (20, 40, 60, 80, 100, 120) for k in range(1, 8)} | {201}
    assert all(t != 4076 or s not in defined for t, s in failures)
    assert vtecs
    v = vtecs[-1]
    assert v.layers and all(layer.complete for layer in v.layers)
    day = v.epoch.second % 86400
    # the sub-solar longitude at that time (GPS ≈ UT to the second here)
    subsolar = math.radians(180.0 - day / 240.0)
    grid = []
    for lat in range(-80, 81, 10):
        for lon in range(-180, 180, 15):
            total = 0.0
            for layer in v.layers:
                total += iono.vtec(layer.degree, layer.order, layer.cosine, layer.sine,
                                   math.radians(lat), math.radians(lon), day)
            grid.append((total, lat, lon))
    values = [g[0] for g in grid]
    peak = max(grid)
    lon_from_sun = (math.radians(peak[2]) - subsolar + math.pi) % (2 * math.pi) - math.pi
    print(f"VTEC over the grid: min {min(values):.1f} max {peak[0]:.1f} TECU at {peak[1]}, {peak[2]}; "
          f"{math.degrees(lon_from_sun):+.0f} deg from the sub-solar longitude")
    assert -5.0 < min(values) and max(values) < 200.0 and peak[0] > 5.0
    assert abs(lon_from_sun) < math.radians(90.0)
