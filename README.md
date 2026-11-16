# ssrcore

An engine for GNSS state-space corrections. It takes corrections (SPARTN or
IGS SSR) and broadcast ephemerides as typed messages, each with an explicit
time. For any satellite, time and receiver position it then gives:

- the satellite's state at transmit time;
- every correction term of the line of sight, each separate, with its
  validity, its quality and the source that gave it.

A virtual base station built on the same engine synthesises RTCM 3
observations for receivers that take RTK corrections.

Version 0.1.0 ([RELEASES.md](RELEASES.md)). Licence: MIT.

```text
                 correction service, receiver, network
                                   │ bytes
                                   ▼
                 gnsswire (SPARTN, RTCM 3, UBX, NMEA codecs)
                                   │ typed messages
                                   ▼
                    ┌────────── ssrcore ──────────┐
                    ▼                             ▼
        line-of-sight records        RTCM 3 (MSM, 1005, 1033, 1230)
```

The codecs are the `gnsswire` distributions, which turn bytes into frozen
dataclasses of raw integers and back. ssrcore does no bit parsing. It also
does no I/O: it opens no socket or file and reads no clock, so every run can
be replayed from recorded inputs.

## Correction sources

Per system, section 2.2.1 of [docs/DESIGN.md](docs/DESIGN.md) gives every
signal each source is taken for, what the record and the virtual base
station give for it, and the ephemerides and satellites accepted.

**SPARTN v2.0.3**:
- Messages: OCB (orbit, clock, biases), HPAC (troposphere and ionosphere)
  and GAD (areas), for GPS, GLONASS, Galileo, BeiDou and QZSS.
- Delivery: over IP, or over L-band inside UBX-RXM-PMP.
- Decryption is the caller's, with the key it holds (`gnsswire`'s
  `spartn.decrypt`).

**IGS SSR, RTCM 4076**: orbit, clock, combined, high-rate clock, code and
phase bias, URA and VTEC messages, through `adapters.igs_ssr` (rules IG1 to
IG12 of [docs/DESIGN.md](docs/DESIGN.md)).

Each adapter rule is tested with synthetic messages, and
`tools/sp3_compare.py` compares a recorded stream's corrected orbits and
clocks with SP3 products (DESIGN 5.2).

Before using IGS SSR:
- A stream's satellite messages are used only once the caller declares it
  phase-centre (`IgsReference`). A centre-of-mass stream is refused.
- BeiDou GEO orbit corrections are left out.
- The SBAS sub-types are refused.
- The engine applies a stream's clocks as sent. With a stream whose clocks
  are phase-defined, pseudoranges corrected with its clocks and code biases
  are off by a constant per satellite.

## Install

ssrcore needs Python 3.11 or later and the `gnsswire` codecs, which are
installed as its dependencies. It is distributed as source, so installing
it compiles the C99 kernels and needs a C compiler. The build fetches its
own requirements (setuptools and Cython).

- Linux: gcc or clang (for example the `build-essential` package), and the
  Python headers (`python3-dev`) for a distribution's Python.
- macOS: Apple's Command Line Tools (`xcode-select --install`).
- Windows: the Microsoft C++ Build Tools (MSVC). The build has not yet been
  run on Windows.

```sh
pip install ssrcore
```

For development, check out `ssrcore` and `gnsswire` side by side. `uv sync`
installs the codecs from `../gnsswire` (`[tool.uv.sources]` in
`pyproject.toml`) and builds the kernels in place:

```sh
uv sync
uv run pytest
```

To install from a local source distribution, with the codecs' wheels built
from `../gnsswire`:

```sh
uv build --sdist --no-sources --out-dir dist
uv build --directory ../gnsswire --all-packages --out-dir "$PWD/dist/gnsswire"
python -m venv env
env/bin/pip install --find-links dist/gnsswire dist/ssrcore-0.1.0.tar.gz
```

## Feed a stream

`Framer` turns pieces of byte streams into gnsswire messages, one buffer per
named transport. `Engine.feed` takes each message with the GPS time it was
received and says what it changed.

`examples/one_satellite.py` builds a short synthetic stream with gnsswire:
- an RTCM 1019 ephemeris of GPS 5, from the public broadcast of
  2024-01-01;
- one SPARTN set: GAD, HPAC and OCB.

The blocks below run in order, from the repository's root with `examples/`
on the import path.

```python
import sys
sys.path.insert(0, "examples")
from one_satellite import stream      # [(bytes, receipt time), ...]

from ssrcore import Engine, Framer, LeapSeconds

engine, framer = Engine(), Framer()
engine.feed(LeapSeconds(18))          # GPS − UTC, as a receiver's UBX-NAV-TIMELS gives it
for piece, received in stream():
    for message in framer.push("ip", piece, received):
        result = engine.feed(message, received)
        print(type(message).__name__, result.outcome.value)
```

## Query a line of sight

```python
import math
from ssrcore import GnssTime, Sat, Signal
from ssrcore.toolkit import geodesy

t_rx = GnssTime(2295, 87400.0)        # GPS week and second of week
receiver = geodesy.geodetic_to_ecef((math.radians(-51.4), math.radians(-162.7), 100.0))
los = engine.query(Sat.gps(5), t_rx, receiver, signals=(Signal.GPS_L1C, Signal.GPS_L2W))

assert los.available and los.corrected and los.missing == {}
print(f"elevation {math.degrees(los.geometry.elevation):.2f} deg, range {los.geometry.range:.3f} m")
print(f"satellite clock {los.satellite.clock * 1e6:.6f} us, troposphere {los.troposphere.slant:.4f} m")
for signal in los.signals:
    print(signal.code, f"code bias {los.code_bias[signal].value:+.3f} m,",
          f"ionosphere {los.ionosphere.code_delay[signal]:.4f} m")
```

A term the engine cannot give is absent, never zero. `los.missing` names
each absent term with its reason: no correction, too old, outside the area,
and so on. `uv run python examples/one_satellite.py` prints every term of
the record.

## Virtual base station

`VirtualBase` places a station near the rover and keeps it fixed while the
rover stays within 5 km. At each epoch it gives the station's messages, in
send order:
- 1005 (position), 1033 (antenna and receiver), and 1230 when GLONASS is
  observed;
- one MSM7 per constellation (MSM4 on request).

```python
from gnsswire import rtcm3
from ssrcore.model.types import VrsConfig
from ssrcore.vrs.station import VirtualBase

station = VirtualBase(engine, VrsConfig(rate=1.0))
messages = station.epoch(t_rx, receiver)          # rover position, ECEF
wire = b"".join(rtcm3.encode(m) for m in messages)
print([m.number for m in messages], len(wire), "bytes")
```

`examples/vrs_stream.py` runs five minutes of a synthetic service with
seven GPS satellites and one GLONASS satellite. It includes an issue change,
a new bias arc and a station move.

## Atmosphere providers

The troposphere and ionosphere come from a list of providers,
`Config.atmosphere`. For each query the engine asks them in order and takes
the first that answers. The record names that provider (`origin`) and gives
its `sigma`.
- `SpartnSource()` is the satellite stream's SPARTN atmosphere, and the
  default.
- `ModelSource()` is the UNB3m troposphere, an a-priori model.
- `CallerSource(name)` is an object the caller passes to `Engine`.
  `model.igs_vtec.IgsVtecIonosphere`, the IGS SSR VTEC, is one such object.

DESIGN 3.6 gives the interface. Here a receiver outside the SPARTN area gets its troposphere from the
caller's provider, and its ionosphere from nobody:

```python
from ssrcore import CallerSource, Config, Kind, ModelSource, Origin, SpartnSource, Term
from ssrcore.model.providers import Zenith


class FixedTroposphere:
    """A caller's provider: the same zenith delays everywhere (an illustration)."""

    name, kind, terms = "fixed", Kind.MODEL, frozenset({Term.TROPOSPHERE})

    def troposphere(self, where):
        return Zenith(hydrostatic=2.30, wet=0.10, origin=Origin(self.name, self.kind), sigma=0.05)

    def ionosphere(self, where, sat, look):
        raise NotImplementedError("not in terms, so never asked")


config = Config(atmosphere=(SpartnSource(), CallerSource("fixed"), ModelSource()))
engine2, framer2 = Engine(config, providers={"fixed": FixedTroposphere()}), Framer()
engine2.feed(LeapSeconds(18))
for piece, received in stream():
    for message in framer2.push("ip", piece, received):
        engine2.feed(message, received)

elsewhere = geodesy.geodetic_to_ecef((math.radians(-48.0), math.radians(-162.7), 100.0))
los2 = engine2.query(Sat.gps(5), t_rx, elsewhere, signals=(Signal.GPS_L1C,))
print(los2.troposphere.origin.provider, f"{los2.troposphere.slant:.4f} m")
print({term.value: reason.value for (term, _), reason in los2.missing.items()})
```

## What it does not do

- **Positioning.** No position, velocity or time solution. What the caller
  builds on the record is outside this package.
- **I/O.** No sockets, NTRIP, serial ports, files or clocks. The caller
  moves bytes and passes times.
- **Bit parsing or decryption.** These are gnsswire's. SPARTN keys stay with
  the caller.
- **Other formats.** No RTCM SSR messages other than IGS SSR 4076, no SPARTN
  BPAC (gnsswire decodes it; the engine refuses it), and no SBAS.
- **Centre-of-mass products.** These need satellite antenna offsets and
  attitude, which no input carries. The record is for the antenna phase
  centre.
- **Ephemerides for the RTK receiver.** The virtual base station sends
  none; the receiver uses its own broadcast data.
- **Single precision.** The C kernels compile unchanged for
  microcontrollers (ESP32-S3, ESP32, ESP32-C3, C6, Cortex-M4F) and stay in
  double precision there (DESIGN 4, Q17).

## Tests

```sh
uv run pytest
```

`uv run pytest` is the gate. It runs on public data only:
- worked examples from the specifications;
- public broadcast and product values;
- synthetic messages built with gnsswire;
- the examples above, this README's Python blocks included.

Tests that read local recordings of public streams, and cross-compilers
that are not installed, skip with the reason. `uv run pytest tests/test_portability.py -rP` prints the kernels'
sizes and operation counts on each target found.

## Design and documents

- [docs/DESIGN.md](docs/DESIGN.md) is the design. It gives the
  specification section behind every rule, and the open questions.
- [Known limitations](docs/DESIGN.md#29-known-limitations), section 2.9 of
  the design, lists what the code does not do as it stands.
- [docs/ROADMAP.md](docs/ROADMAP.md) gives the work intended and what is not planned.
- [RELEASES.md](RELEASES.md) gives the notes of each release.
- [refs/README.md](refs/README.md) lists the specifications cited. They
  are not committed.

| Path | Layer | Contents |
|---|---|---|
| `kernels/` | kernels | C99 sources and headers: stateless arithmetic, compilable unchanged for a microcontroller |
| `src/ssrcore/_kernels/` | kernels | Cython bindings, one extension module per kernel |
| `src/ssrcore/toolkit/` | public | stateless functions: time scales, geodesy, orbits, troposphere, ionosphere, tides, wind-up, relativity |
| `src/ssrcore/model/` | model | stores, selection, consistency rules, atmosphere providers |
| `src/ssrcore/adapters/` | model | gnsswire messages to typed, scaled inputs; the framer |
| `src/ssrcore/engine.py`, `record.py` | interface | feed and query; the line-of-sight record |
| `src/ssrcore/vrs/` | interface | virtual base station |
| `tools/` | | the SP3 and tide checks, the C query and portability report, the astro table generator |

## Licence

MIT ([LICENSE](LICENSE)). [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
lists the third-party values the source holds, whose they are and their
terms.
