# ssrcore design

The design of ssrcore: the public interface first, with its known
limitations gathered in section 2.9, then the internal
structure and the rules that keep data consistent, then numerics and
portability, testing, and the questions the specifications leave open.

Every rule here, and every function of the code, is one of two kinds;
there is no third (section 5.1, *Specification-driven*):

- **[spec: DOCUMENT section]**: the rule implements that public document,
  and the code that implements it cites the same place;
- **[design: reason]**: the specifications do not decide, and this is the
  choice made, with its reason.

Document names are short forms; [refs/README.md](../refs/README.md) lists
each. "SPARTN" is SPARTN v2.0.3 (November 2025), "RTCM" is RTCM 10403.3,
"IS-GPS-200N with IRN-IS-200N-003 and -004" is IS-GPS-200 revision N (1
August 2022) as the change notices IRN-IS-200N-003 (27 September 2023) and
IRN-IS-200N-004 (16 June 2026) amend it, "TM-23" is ESA GNSS Data Processing Vol. I
(TM-23/1), "TN36" is IERS Technical Note 36 chapter 7, "IGS SSR" is IGS
State Space Representation Format v1.00 (October 2020). Rules are named by
identifier (S1, C4, A9, …, section 3.3) and open questions by Q and number
(section 6); a question that has been settled keeps its identifier where its
answer is stated.

## 1. Purpose and scope

ssrcore sits on the codecs:

| Package | Owns | Does not do |
|---|---|---|
| gnsswire (`../gnsswire`) | bytes to frozen dataclasses of raw integers and back: SPARTN, RTCM 3, UBX, NMEA | scaling, physics, state |
| **ssrcore** | scaling, stores, consistency, satellite state, every correction term, virtual base station observations | bit parsing, I/O, clocks, positioning |

Its outputs are the line-of-sight record (section 2.4) and, from the same
engine, the RTCM 3 observations of a virtual base station (section 2.6).
What the caller builds on them is outside its scope.

Constellations: GPS, GLONASS, Galileo, BeiDou and QZSS, the five SPARTN
covers (SPARTN Table 6.1). SPARTN messages handled: OCB (SM 0-0 to 0-4),
HPAC (SM 1-0 to 1-4) and GAD (SM 2-0). BPAC (SM 3-0) is decoded by gnsswire
and refused by the engine until a `bpac` kernel (SPARTN 8.13) and its
adapter exist. Dynamic key and authentication messages (SM 4) are not the
engine's: decryption is done before messages reach it, with the key the
caller holds (gnsswire `spartn.decrypt`). IGS SSR (RTCM 4076) is taken
through its own adapter (section 3.3, rules IG1 to IG12).

Targets: Linux, macOS and Windows hosts; the kernels compile unchanged for a
microcontroller with a single-precision FPU or none (section 4).

## 2. Public interface

### 2.1 Shape

```python
from ssrcore import Engine, GnssTime, Sat, Signal

engine = Engine(config)                       # no I/O, no clock
engine.feed(message, received=t)              # any input of section 2.2, with its receipt time
los = engine.query(Sat.gps(5), t_rx, receiver_ecef, signals=(Signal.GPS_L1C, Signal.GPS_L2W))
if los.available:
    los.satellite.position, los.troposphere.slant, los.code_bias[Signal.GPS_L1C]
```

- `Engine(config, providers=None)` takes a frozen `Config` (section 2.5
  lists its fields) and, for each `CallerSource` of `config.atmosphere`,
  the caller's provider object (`providers={name: object}`, rule A9).
  Nothing else configures it; there is no global state. [design: replay
  needs the whole state to be the inputs and the configuration; a caller's
  provider is an input the engine cannot record, so a replay needs the same
  provider, a deterministic function of what its caller fed it]
- `feed(message, received=None)` takes one typed input and, when the
  caller has it, the GPS time it was received. Without a receipt time the
  engine's time base derives "now" from the data itself (rule T5), so a
  recording that holds only bytes replays exactly. It returns a
  `FeedResult`: the outcome (`stored`, `superseded`, `held`, `duplicate`,
  `counted`, `rejected`), the reason, and the satellites and areas the
  message touched. It never raises for bad data; it raises only for a
  programming error (a type it does not know). [design: a stream with
  errors must not stop the engine; the result makes every rejection
  visible; recordings often carry no receipt times]
- `query(sat, t_rx, receiver, signals)` is pure with respect to the stores:
  it reads them and changes nothing, so the same feed sequence and the same
  queries give the same results in any order of queries. [design:
  replayability and safe concurrent reading] It keeps caches of values
  that are functions of the stores and of its arguments alone (section
  2.8), so a cache changes the cost of a query and never its result.
- `engine.status()` returns, without a satellite, the active stream's
  quality and every stream's state (section 2.5), the time base (current
  time, leap seconds, messages held and dropped), the count of messages
  fed, the configured satellite stream and the providers' names.

**I/O at the edge.** The core takes typed messages and explicit times. It
opens no socket or file and reads no clock, so every run can be replayed
from recorded inputs; moving bytes and passing times is the caller's.
`tests/test_rules.py` checks it (5.1). [design: replay, and a core that runs
unchanged wherever its caller does] Times are explicit everywhere: the
engine never asks what time it is. "Now" for ages and expiry is the `t_rx`
of the query. [design: *I/O at the edge*]

### 2.2 Inputs

**Wire formats are gnsswire's.** No bit parsing is done here: the codecs
are the `gnsswire` distributions (section 1), the package's dependencies,
which turn bytes into frozen dataclasses of raw integers and back. For
development they come from a sibling checkout at `../gnsswire`
(`[tool.uv.sources]` in `pyproject.toml`). [design: one codec per wire
format, shared by every package that reads it, and the engine sees only
typed values]

The engine accepts gnsswire messages directly and converts them with the
adapters in `ssrcore.adapters` (raw integers to SI units; no bit parsing).
Each adapter is also public, so a caller can convert once and feed the typed
input.

| Input | gnsswire type | Adapter output | Basis |
|---|---|---|---|
| SPARTN OCB, SM 0-0 to 0-4 | `spartn.Message` with body `Ocb` | `OcbInput` | SPARTN 6.3, Table 6.2 |
| SPARTN HPAC, SM 1-0 to 1-4 | `spartn.Message` with body `Hpac` | `HpacInput` | SPARTN 6.3, Table 6.2 |
| SPARTN GAD, SM 2-0 | `spartn.Message` with body `Gad` | `GadInput` | SPARTN 6.3, Table 6.2 |
| IGS SSR orbit, clock, combined, high-rate clock, code and phase biases, URA (RTCM 4076 sub-types 21 to 27, 41 to 47, 61 to 67, 81 to 87, 101 to 107) | `rtcm3.IgsOrbit`, `IgsClock`, `IgsCombined`, `IgsHighRateClock`, `IgsCodeBias`, `IgsPhaseBias`, `IgsUra` | `OcbInput` of an `IgsStreamId`, through `adapters.igs_ssr.IgsSsr` (rules IG1 to IG12); the URA, and an orbit waiting for its clock, are kept by the adapter | IGS SSR 4, 8.1 to 8.3, Table 7; RTCM 3.5.13 |
| IGS SSR reference point of a stream | `IgsReference(provider, solution, point)` | itself | IG11 |
| IGS SSR VTEC (4076 sub-type 201) | `rtcm3.IgsVtec` | `VtecInput`, fed to each `model.igs_vtec.IgsVtecIonosphere` among the providers (A9); with none, rejected, `not_an_input` | IGS SSR 4.5.1, 8.4.1 |
| IGS SSR SBAS sub-types (121 to 127) | the gnsswire message | none | rejected, `unsupported_system` |
| GPS LNAV ephemeris | `rtcm3.Msg1019` | `KeplerEphemeris` | RTCM Table 3.5-21; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4 |
| GLONASS ephemeris | `rtcm3.Msg1020` | `GlonassEphemeris` | RTCM Table 3.5-22; GLONASS ICD 4.4 |
| BeiDou D1/D2 ephemeris | `rtcm3.Msg1042` | `KeplerEphemeris` | RTCM Table 3.5-113; BDS ICD 5.2.4 |
| Galileo F/NAV, I/NAV | `rtcm3.Msg1045`, `Msg1046` | `KeplerEphemeris` | RTCM Tables 3.5-110, 3.5-111; Galileo OS SIS ICD 2.2 5.1 |
| QZSS ephemeris | `rtcm3.Msg1044` | `KeplerEphemeris` | RTCM Table 3.5-112 (Q14) |
| navigation words from a receiver | `ubx.RxmSfrbx` (from the framing helper as `Carried(transport, message)`) | through `adapters.navdata`, per transport, to the RTCM message of its system (1019, 1020, 1042, 1046), then as that row; leap seconds as a `LeapSeconds`; a word that fails its code as a `FrameFailure` (Q15); a word of a signal for which gnsswire-nav lays out no table (every QZSS word, GPS L2C and L5, Galileo E5a F/NAV, SBAS, NavIC, and those F9 HPG 1.51 and the X20 add: Galileo E6, BeiDou B3I, BeiDou B1C and B2a by pilot and data component, and the X20's QZSS L1C/B) rejected, `unsupported_system` | u-blox F9 HPG 1.51 1.5.4 and 3.17.9, u-blox X20 HPG 2.11 1.5.4 and 3.16.8 (gnssId, sigId); u-blox M8 Receiver Description UBX-13003221 R28 chapter 10 and the ICDs (IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3, Galileo OS SIS ICD 2.2 4.3.5, BDS ICD 5.2.4, GLONASS ICD 4.4, 4.5) through gnsswire-nav |
| leap seconds | `ubx.NavTimels` or a `LeapSeconds` | `LeapSeconds` | u-blox F9 HPG 1.51 3.15.26, u-blox X20 HPG 2.11 3.14.22; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 for the effectivity time |
| UTC date and time | `nmea.Zda` | `UtcTime` | u-blox F9 HPG 1.51 2.7.21.1, u-blox X20 HPG 2.11 2.7.20.1 (ZDA), standing in for NMEA 0183, which is not public |
| rover position | `nmea.Gga` | `RoverPosition` | u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.7.5.1, 2.5.3, 2.5.5 (GGA); a GGA behind an NMEA TAG block, the secondary output (both 2.8.1.1), is read as any other; used by the VRS (2.6.1) and kept as `Engine.rover` for the caller |
| other systems' ephemerides (RTCM 1041 NavIC, 1043 SBAS) | none: gnsswire has no table for them, so the framing helper gives `FrameFailure` of kind `unsupported` | itself | counted, `frame_failure` (S7, section 3.5) |
| header of an OCB message whose subtype has no table (SM 0-5 to 0-15) | `OcbHeaderInput`, from the framing helper | itself | SPARTN 8.4 item 5; rule S5 (Q19) |
| SPARTN BPAC, SM 3-0 | `spartn.Message` with body `Bpac` | none yet | rejected, `not_an_input` (section 1) |
| station and observation messages (RTCM 1005, 1006, 1033, 1230, MSM1 to MSM7), RXM-RAWX, an RXM-PMP fed whole | the gnsswire message | none | rejected, `not_an_input` [design: a receiver port or a caster carries them beside the inputs, and the framing helper passes on every message with a table; refusing them visibly keeps `feed` total over what the framer yields] |
| undecodable frame | `FrameFailure(kind, protocol, message_type, subtype, received)` | itself | section 3.5 |

SPARTN frames also arrive inside UBX-RXM-PMP (L-band); unwrapping the user
data is the framing helper's (below), not the engine's. Navigation words
(`RxmSfrbx`) are joined into ephemerides per transport, so the framing
helper gives them with the transport they came on (`Carried`); a bare
`RxmSfrbx` fed to the engine belongs to the unnamed transport `""`. [design:
the words of two receivers must not meet, and only the framer knows the
transport] gnsswire-nav decodes the words and checks their codes;
`adapters.navdata` joins them, and its module documentation gives the
assembly rules (issue of data, stale parts, GLONASS rows) and the fields one
source has and the other lacks (Q15). [design: the RTCM message is the
hand-over, so one code path stores an ephemeris whatever carried it]

The ephemeris adapter (`adapters.rtcm3`) takes, besides the message, a GPS
time near its receipt (the receipt time, or the time base's current time,
rule T5), which resolves the week rollover and the week or day of a
reference time, and for GLONASS the leap seconds (rule T4); the offsets
between time scales are the toolkit's. A 1019 whose IODE differs from the 8
LSBs of its IODC is refused [spec: IS-GPS-200N with IRN-IS-200N-003 and -004
20.3.3.4.1, a cutover], as is a 1044 likewise [spec: RTCM DF456]; `feed`
reports such a message rejected, `inconsistent`.

The SPARTN adapters treat a block holding a field whose value is special or
outside its range as not sent, so that the other blocks of the message stay
usable, and leave out biases at spare mask positions [spec: SPARTN 6.3,
Table 6.2; design: a value that is not a number cannot be applied, and one
bad block need not discard the satellite's others].

**Framing helper.** Recorded and live inputs are often pieces of byte
streams: a frame may start in one piece and end in a later one, a piece may
hold several frames and bytes that belong to none, and several transports
(for example an IP and an L-band correction stream, and a receiver port) are
interleaved. `ssrcore.adapters.framing.Framer` holds one buffer per named
transport, runs gnsswire's `parse` over it, and yields typed messages and
`FrameFailure`s in arrival order. Each transport may carry every codec (told
apart by the first byte of a frame), as a receiver port does; the SPARTN
bytes of RXM-PMP user data are framed as an inner stream of the transport
they came on; `flush` ends a transport; UBX and NMEA messages without a
gnsswire table are not inputs and are skipped (section 3.5). It is at the
edge: bytes in, messages out, no I/O and no clock; the engine core never
sees bytes. [design: *Wire formats are gnsswire's* (above) and *I/O at the
edge* (2.1); a shared framer makes every caller frame the same way]

*Delay.* A frame is released as soon as its last byte has arrived and its
check passes; bytes are held only while a frame that starts at them is
incomplete. A false start (a start byte inside other data whose length field
claims more than has arrived) therefore holds every byte after it until
enough have arrived to fail the check: up to 65 543 bytes for UBX, 1 029 for
RTCM, about 1 100 for SPARTN. The caller bounds this with
`Framer(max_frame=…, max_wait=…, protocols=…)`: a start whose frame would
exceed `max_frame` bytes, or whose first byte was received more than
`max_wait` seconds before the piece being pushed, is taken as false and
scanning moves past it; `protocols` names, for a transport known to carry
only some protocols, those protocols, so that a start byte of another is
never a frame there. The defaults (None) give each protocol's own maximum
and every protocol on every transport. [design: the limits act only on the
bytes and on the receipt times the caller passes, so a replay of the same
pieces with the same times gives the same messages, and a replay without
receipt times is not affected by `max_wait`] Frames found beyond a pending
start are not released before it is resolved [design: UBX RXM-PMP payloads
contain whole SPARTN frames, so a scan past a pending UBX start would
release the inner SPARTN frames a second time and out of order].

SPARTN inputs carry their own time tag (TF009); the engine resolves a 16-bit
tag against the receipt time, or the time base's current time when there is
none (rules T2, T5).

#### 2.2.1 Per constellation

What the package takes and gives for each system, gathered from the rules
that state it (C9, E1 to E5, T1 to T4, IG1, IG4, IG7, IG8, IG12), the table
of 2.2 and section 2.6. One table per system, one row per signal the code
knows for it: every signal of a SPARTN bias mask, of an IGS SSR identifier
or of an MSM signal ID, named by its RINEX 3 code without the observation
type. `tests/test_constellations.py` checks every row against the code's
own tables (`model.types.SPARTN_SIGNALS`, `adapters.igs_ssr.signal`,
`vrs.msm.SIGNAL_ID`) and the frequencies against `Signal.frequency`.

- *Band, MHz*: the band and the carrier frequency the record uses;
  "none" where the code has no frequency for the band, so that no `Signal`
  can name the signal.
- *SPARTN*: the bias mask position (Table 6.2) at which an OCB bias block
  gives the signal's code bias (C code) and phase bias (L code); both are
  taken (C9). A bias at a spare position is left out (2.2).
- *IGS SSR*: the IDF024 identifier under which the signal's code bias
  (IDF025) and phase bias (IDF028) are both taken (IG6, IG7); "RTCM" where
  IGS SSR reserves the identifier and RTCM defines it; "left out" where the
  adapter drops and counts it. Every identifier not in a table is reserved
  by IGS SSR and not defined by RTCM, and is left out and counted (IG7).
- *Line-of-sight record*: the per-signal terms of a record that asks for
  the signal, or that has a bias for it when no signal is asked (2.3), each
  when its inputs are present (2.4): b_s (`code_bias[s]`); B_s
  (`phase_bias[s]`, with its fix flag, continuity and discontinuity count);
  I_s, the ionosphere delay at the signal's frequency
  (`ionosphere.code_delay[s]`, the phase delay its negative), with the
  frequency itself (`frequency[s]`); and λs·w, the wind-up in metres, which
  is not a field of the record: `phase_model(s)` forms it from
  `frequency[s]` and the satellite's `modelled.windup`, one value in
  cycles. "no bias": no input gives a bias for the signal, so b_s and B_s
  are always absent (`no_bias`, or the selection's reason when the
  satellite has no bias block, E4).
- *Virtual base station*: the MSM signal ID (RTCM Tables 3.5-91, 3.5-96,
  3.5-99, 3.5-105, 3.5-108) under which the signal is sent. It is sent when
  its satellite is observed (orbit and clock corrected, elevation at least
  `elevation_mask`, the troposphere a correction, the ionosphere a
  correction from the stream of the code biases, the tide present, a
  satellite ID in the MSM's range) and the record has the signal's code
  bias and frequency; its phase is sent when the record also has the phase
  bias and the wind-up and the bias is fixed or `float_phase` is set
  (2.6.3). "never sent": the ID exists, but no input gives the signal a
  bias, or its band has no frequency. "no MSM ID": RTCM 10403.3 defines
  none, so the signal is not sent.

For every system: a signal is never replaced by another code of its band
(C9, Q10); the broadcast group delays (TGD, BGD, Δτn) are stored with the
ephemeris and applied in no case, with code biases or without them (E4);
the HPAC of a system gives the ionosphere of that system's satellites
(A6), and any system's HPAC the troposphere (A5); the offsets between
time scales that the navigation messages broadcast (GGTO, τc, τGPS, A0GPS,
A1GPS, the UTC A0, A1) have no input and are not applied: GST, BDT and QZSS
time become GPS time by whole seconds (T1), GLONASS time by the leap
seconds alone (T4), and `toolkit.timescale` evaluates the offsets for a
caller that has them (2.9).

**GPS** (SPARTN SM 0-0 and 1-0; IGS SSR sub-types 21 to 27)

| Signal | Band, MHz | SPARTN | IGS SSR | Line-of-sight record | Virtual base station |
|---|---|---|---|---|---|
| `1C` | L1, 1575.42 | bit 0: C1C (SF027), L1C (SF025) | ID 0 | b_s, B_s, I_s, λs·w | ID 2 |
| `1L` | L1, 1575.42 | — | ID 4 | b_s, B_s, I_s, λs·w | ID 31 |
| `1P` | L1, 1575.42 | — | ID 1 | b_s, B_s, I_s, λs·w | ID 3 |
| `1S` | L1, 1575.42 | — | ID 3 | b_s, B_s, I_s, λs·w | ID 30 |
| `1W` | L1, 1575.42 | — | ID 2 | b_s, B_s, I_s, λs·w | ID 4 |
| `1X` | L1, 1575.42 | — | — | I_s, λs·w; no bias | ID 32: never sent, no bias |
| `2C` | L2, 1227.6 | — | ID 5 | b_s, B_s, I_s, λs·w | ID 8 |
| `2D` | L2, 1227.6 | — | ID 6 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `2L` | L2, 1227.6 | bit 2: C2L (SF027), L2L (SF025) | ID 8 | b_s, B_s, I_s, λs·w | ID 16 |
| `2P` | L2, 1227.6 | — | ID 10 | b_s, B_s, I_s, λs·w | ID 9 |
| `2S` | L2, 1227.6 | — | ID 7 | b_s, B_s, I_s, λs·w | ID 15 |
| `2W` | L2, 1227.6 | bit 1: C2W (SF027), L2W (SF025) | ID 11 | b_s, B_s, I_s, λs·w | ID 10 |
| `2X` | L2, 1227.6 | — | ID 9 (RTCM DF380) | b_s, B_s, I_s, λs·w | ID 17 |
| `5I` | L5, 1176.45 | — | ID 14 | b_s, B_s, I_s, λs·w | ID 22 |
| `5Q` | L5, 1176.45 | bit 3: C5Q (SF027), L5Q (SF025) | ID 15 | b_s, B_s, I_s, λs·w | ID 23 |
| `5X` | L5, 1176.45 | — | — | I_s, λs·w; no bias | ID 24: never sent, no bias |

- *Ephemerides*: RTCM 1019 (LNAV), or the L1 C/A LNAV subframes of
  UBX-RXM-SFRBX, assembled into a 1019 (2.2); L2C and L5 words are
  refused, `unsupported_system`. A 1019 whose IODE differs from the 8 LSBs
  of its IODC is refused, `inconsistent` (2.2). Valid over half the fit
  interval on either side of toe (2.5); unused when the SV health is not 0
  (E5).
- *Issue of data*: a SPARTN orbit names the IODE (SF018) and ephemeris type
  SF016 0 (L1 C/A), the 1019; another type matches nothing (E1, E2). IGS
  SSR: IDF012 is the IODE of the 1019 (IG8).
- *Satellites*: PRN 1 to 64 from SPARTN (SF011); 1 to 63, and 0 as 64,
  from IGS SSR (IG1); 1 to 63 from the 1019 (DF009, 0 refused), so a
  correction for 64 finds no ephemeris (2.9); PRN 1 to 63 in the MSM
  (2.6.3).
- *Time*: GPS time.
- *Limits*: IGS SSR identifiers 12, 13 and 16 to 31 are reserved and left
  out (IG7). 2D has no MSM signal ID and is not sent (2.6.3).

**GLONASS** (SPARTN SM 0-1 and 1-1; IGS SSR sub-types 41 to 47)

| Signal | Band, MHz | SPARTN | IGS SSR | Line-of-sight record | Virtual base station |
|---|---|---|---|---|---|
| `1C` | G1, 1602 + 0.5625·K | bit 0: C1C (SF028), L1C (SF026) | ID 0 | b_s, B_s, I_s, λs·w | ID 2 |
| `1P` | G1, 1602 + 0.5625·K | — | ID 1 | b_s, B_s, I_s, λs·w | ID 3 |
| `2C` | G2, 1246 + 0.4375·K | bit 1: C2C (SF028), L2C (SF026) | ID 2 | b_s, B_s, I_s, λs·w | ID 8 |
| `2P` | G2, 1246 + 0.4375·K | — | ID 3 | b_s, B_s, I_s, λs·w | ID 9 |
| `3I` | G3: none | — | ID 8: left out, no frequency | none: no frequency | — |
| `3Q` | G3: none | — | ID 9: left out, no frequency | none: no frequency | — |
| `4A` | G1a: none | — | ID 4: left out, no frequency | none: no frequency | — |
| `4B` | G1a: none | — | ID 5: left out, no frequency | none: no frequency | — |
| `6A` | G2a: none | — | ID 6: left out, no frequency | none: no frequency | — |
| `6B` | G2a: none | — | ID 7: left out, no frequency | none: no frequency | — |

- *Ephemerides*: RTCM 1020, or strings 1 to 5 of L1OF and L2OF in
  UBX-RXM-SFRBX, assembled into a 1020 (2.2, `adapters.navdata`). Held
  until the leap seconds are known (T4). A tb outside 1 to 95 or a channel
  outside −7 to 6 is refused. Valid within 15 min + P1/2 of tb (2.5);
  unused when the Bn MSB or ln is set (E5). The orbit is integrated by RK4
  (E3); the clock has no relativistic term (C2).
- *Issue of data*: a SPARTN orbit names tb, the 7 bits of SF019, and type
  SF017 0, the 1020 (E1, E2, Q3). IGS SSR: IDF012 with bit 7 clear is tb;
  with bit 7 set the orbit is left out (IG8, IG12).
- *Satellites*: slot 1 to 63 from SPARTN (SF012); 1 to 63 from IGS SSR,
  whose 64 has no slot and is left out (IG1); 1 to 63 from the 1020 (DF038,
  0 refused); slot 1 to 24 in the MSM. The
  frequency channel K is the stored ephemeris's: DF040 of the 1020, or for
  navigation words the receiver's freqId or else the almanac strings (an
  ephemeris without one is refused). Without a stored ephemeris a GLONASS
  signal has no frequency (`no_frequency`). The MSM carries K + 7
  (DF419), and the 1230 the L1 and L2 C/A biases zero, aligned (2.6.4).
- *Time*: GLONASS time, UTC(SU) + 3 h, converted with the time base's leap
  seconds; without them a query answers `no_leap_seconds` (T4). A 32-bit
  SPARTN tag counts the GLONASS clock from 2010 (T3, Q2); an IGS SSR epoch
  is GPS time (IG3).
- *Limits*: G3, G1a and G2a have no frequency here, so IGS SSR identifiers
  4 to 9 are left out (IG7); 10 to 31 are reserved and left out (IG7).
  τc and τGPS have no input and are not applied (2.9).

**Galileo** (SPARTN SM 0-2 and 1-2; IGS SSR sub-types 61 to 67)

| Signal | Band, MHz | SPARTN | IGS SSR | Line-of-sight record | Virtual base station |
|---|---|---|---|---|---|
| `1A` | E1, 1575.42 | — | ID 0 | b_s, B_s, I_s, λs·w | ID 3 |
| `1B` | E1, 1575.42 | — | ID 1 | b_s, B_s, I_s, λs·w | ID 4 |
| `1C` | E1, 1575.42 | bit 0: C1C (SF105), L1C (SF102) | ID 2 | b_s, B_s, I_s, λs·w | ID 2 |
| `1X` | E1, 1575.42 | — | — | I_s, λs·w; no bias | ID 5: never sent, no bias |
| `1Z` | E1, 1575.42 | — | — | I_s, λs·w; no bias | ID 6: never sent, no bias |
| `5I` | E5a, 1176.45 | — | ID 5 | b_s, B_s, I_s, λs·w | ID 22 |
| `5Q` | E5a, 1176.45 | bit 1: C5Q (SF105), L5Q (SF102) | ID 6 | b_s, B_s, I_s, λs·w | ID 23 |
| `5X` | E5a, 1176.45 | — | — | I_s, λs·w; no bias | ID 24: never sent, no bias |
| `6A` | E6, 1278.75 | — | ID 14 | b_s, B_s, I_s, λs·w | ID 9 |
| `6B` | E6, 1278.75 | — | ID 15 | b_s, B_s, I_s, λs·w | ID 10 |
| `6C` | E6, 1278.75 | bit 4: C6C (SF105), L6C (SF102) | ID 16 | b_s, B_s, I_s, λs·w | ID 8 |
| `6X` | E6, 1278.75 | — | — | I_s, λs·w; no bias | ID 11: never sent, no bias |
| `6Z` | E6, 1278.75 | — | — | I_s, λs·w; no bias | ID 12: never sent, no bias |
| `7I` | E5b, 1207.14 | — | ID 8 | b_s, B_s, I_s, λs·w | ID 14 |
| `7Q` | E5b, 1207.14 | bit 2: C7Q (SF105), L7Q (SF102) | ID 9 | b_s, B_s, I_s, λs·w | ID 15 |
| `7X` | E5b, 1207.14 | — | — | I_s, λs·w; no bias | ID 16: never sent, no bias |
| `8I` | E5, 1191.795 | — | ID 11 (RTCM DF382) | b_s, B_s, I_s, λs·w | ID 18 |
| `8Q` | E5, 1191.795 | bit 3: C8Q (SF105), L8Q (SF102) | ID 12 (RTCM DF382) | b_s, B_s, I_s, λs·w | ID 19 |
| `8X` | E5, 1191.795 | — | — | I_s, λs·w; no bias | ID 20: never sent, no bias |

- *Ephemerides*: RTCM 1045 (F/NAV, the health of E5a) and 1046 (I/NAV, the
  health of E1-B and E5b), or the I/NAV pages of E1 and E5b in
  UBX-RXM-SFRBX, assembled into a 1046; E5a F/NAV words, and the E6 words
  (sigId 8 to 10, u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 1.5.4), are refused,
  `unsupported_system`. Valid within 4 h of toe (2.5); unused when a
  signal's health is not 0 or its data validity is 1 (E5).
- *Issue of data*: a SPARTN orbit names the IODnav (SF099) and type SF096,
  0 for F/NAV (1045) or 1 for I/NAV (1046) (E1, E2). IGS SSR corrects I/NAV
  only: IDF012, the 8 LSBs of the IODnav, names the stored 1046 issue whose
  IODnav ends in them; without one the value is taken as sent and matches
  only by chance until the next correction (IG8; 2.9).
- *Satellites*: PRN 1 to 64 from SPARTN (SF093); 1 to 36 from IGS SSR
  (IG1); 1 to 63 from the 1045 and 1046 (DF252, 0 refused), so a correction
  for 64 finds no ephemeris (2.9); PRN 1 to 50 in the MSM.
- *Time*: GST, converted to GPS time by the week offset alone: the scales
  are taken to share whole seconds, and GGTO is not applied (T1; 2.9).
- *Limits*: IGS SSR identifiers 3, 4, 7, 10, 13 and 17 to 31 are reserved
  and left out (IG7). The E5a quasi-pilot of Galileo OS SIS ICD 2.2 has no
  row: no `Signal`, MSM signal ID or bias names it (2.9).

**BeiDou** (SPARTN SM 0-3 and 1-3; IGS SSR sub-types 101 to 107)

| Signal | Band, MHz | SPARTN | IGS SSR | Line-of-sight record | Virtual base station |
|---|---|---|---|---|---|
| `1A` | B1C, 1575.42 | — | ID 15 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `1D` | B1C, 1575.42 | — | ID 9 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `1P` | B1C, 1575.42 | bit 4: C1P (SF106), L1P (SF103) | ID 10 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `2I` | B1I, 1561.098 | bit 0: C2I (SF106), L2I (SF103) | ID 0 | b_s, B_s, I_s, λs·w | ID 2 |
| `2Q` | B1I, 1561.098 | — | ID 1 | b_s, B_s, I_s, λs·w | ID 3 |
| `2X` | B1I, 1561.098 | — | — | I_s, λs·w; no bias | ID 4: never sent, no bias |
| `5D` | B2a, 1176.45 | — | ID 12 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `5P` | B2a, 1176.45 | bit 1: C5P (SF106), L5P (SF103) | ID 13 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `6A` | B3, 1268.52 | — | ID 18 | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `6I` | B3, 1268.52 | bit 3: C6I (SF106), L6I (SF103) | ID 3 | b_s, B_s, I_s, λs·w | ID 8 |
| `6Q` | B3, 1268.52 | — | ID 4 | b_s, B_s, I_s, λs·w | ID 9 |
| `6X` | B3, 1268.52 | — | — | I_s, λs·w; no bias | ID 10: never sent, no bias |
| `7I` | B2b, 1207.14 | bit 2: C7I (SF106), L7I (SF103) | ID 6 | b_s, B_s, I_s, λs·w | ID 14 |
| `7P` | B2b, 1207.14 | bit 5: C7P (SF106), L7P (SF103) | — | b_s, B_s, I_s, λs·w | none: no MSM ID |
| `7Q` | B2b, 1207.14 | — | ID 7 | b_s, B_s, I_s, λs·w | ID 15 |
| `7X` | B2b, 1207.14 | — | — | I_s, λs·w; no bias | ID 16: never sent, no bias |
| `8P` | B2a+b, 1191.795 | bit 6: C8P (SF106), L8P (SF103) | — | b_s, B_s, I_s, λs·w | none: no MSM ID |

- *Ephemerides*: RTCM 1042 (D1 from MEO and IGSO satellites, D2 from GEO,
  both of the B1I ICD), or the D1 and D2 subframes of B1I and B2I in
  UBX-RXM-SFRBX, assembled into a 1042; B1C and B2a words (sigId 5 to 8,
  pilot and data), and the B3I words (sigId 4 and 10; u-blox F9 HPG 1.51 and
  u-blox X20 HPG 2.11 1.5.4), are refused, `unsupported_system`. Valid within 2 h of toe (2.5); unused when SatH1
  is set (E5). The GEO satellites are PRN 1 to 5 and 59 to 63 (BDS ICD
  Table 4-1).
- *Issue of data*: a SPARTN orbit names the IODC, mod(toc/720, 240)
  (SF100), and type SF097 0 (D1) for a MEO or IGSO satellite or 1 (D2) for
  a GEO, both the 1042; types 2 and 3 (B3I) are taken to have no RTCM
  message and match nothing (E1, E2; 2.9). IGS SSR: D2 for GEO and D1
  otherwise; IDF012, mod(toe/720, 240), names the stored issue with that toe; without one the
  value is taken as sent and matches only by chance until the next
  correction (IG8; 2.9).
- *Satellites*: PRN 1 to 64 from SPARTN (SF094); 1 to 63, and 0 as 64,
  from IGS SSR (IG1), whose GEO orbit corrections are left out (IG4, IG12),
  so that a GEO satellite of an IGS SSR stream has no orbit and clock pair
  and its record is the broadcast state (3.4); 1 to 63 from the 1042
  (DF488, 0 refused), so a correction for 64 finds no ephemeris (2.9);
  PRN 1 to 37 in the MSM, so GEO 59 to 63 are not observed.
- *Time*: BDT, GPS time − 14 s (BDS ICD 3.3), converted from SPARTN tags
  (T1) and from IGS SSR epochs (IG3).
- *Limits*: 1A, 1D, 1P, 5D, 5P, 6A, 7P and 8P have no MSM signal ID and
  are not sent (2.6.3). IGS SSR identifiers 2, 5, 8, 11, 14, 16, 17 and 19
  to 31 are reserved and left out (IG7). GEO orbit corrections from IGS
  SSR are left out (IG4; 2.9). BeiDou's A0GPS, A1GPS (BDT relative to GPS
  time, BDS ICD 5.2.4.19) have no input and are not applied (2.9).

**QZSS** (SPARTN SM 0-4 and 1-4; IGS SSR sub-types 81 to 87)

| Signal | Band, MHz | SPARTN | IGS SSR | Line-of-sight record | Virtual base station |
|---|---|---|---|---|---|
| `1C` | L1, 1575.42 | bit 0: C1C (SF107), L1C (SF104) | ID 0 | b_s, B_s, I_s, λs·w | ID 2 |
| `1L` | L1, 1575.42 | — | ID 2 | b_s, B_s, I_s, λs·w | ID 31 |
| `1S` | L1, 1575.42 | — | ID 1 | b_s, B_s, I_s, λs·w | ID 30 |
| `1X` | L1, 1575.42 | — | — | I_s, λs·w; no bias | ID 32: never sent, no bias |
| `2L` | L2, 1227.6 | bit 1: C2L (SF107), L2L (SF104) | ID 4 | b_s, B_s, I_s, λs·w | ID 16 |
| `2S` | L2, 1227.6 | — | ID 3 | b_s, B_s, I_s, λs·w | ID 15 |
| `2X` | L2, 1227.6 | — | — | I_s, λs·w; no bias | ID 17: never sent, no bias |
| `5I` | L5, 1176.45 | — | ID 6 | b_s, B_s, I_s, λs·w | ID 22 |
| `5Q` | L5, 1176.45 | bit 2: C5Q (SF107), L5Q (SF104) | ID 7 | b_s, B_s, I_s, λs·w | ID 23 |
| `5X` | L5, 1176.45 | — | — | I_s, λs·w; no bias | ID 24: never sent, no bias |
| `6E` | L6: none | — | ID 17: left out, no frequency | none: no frequency | — |
| `6L` | L6: none | — | ID 10: left out, no frequency | none: no frequency | ID 10: never sent, no frequency |
| `6S` | L6: none | — | ID 9: left out, no frequency | none: no frequency | ID 9: never sent, no frequency |
| `6X` | L6: none | — | — | none: no frequency | ID 11: never sent, no frequency |

- *Ephemerides*: RTCM 1044 only: QZSS words in UBX-RXM-SFRBX are refused,
  `unsupported_system`, since gnsswire-nav lays none out. A 1044 whose
  IODE differs from the 8 LSBs of its IODC is refused (DF456). Valid within
  1 h of toe (2.5, Q14); unused when one of the five defined health bits
  of DF454 is set (E5). Orbit and clock by IS-GPS-200N with IRN-IS-200N-003
  and -004 (Q14).
- *Issue of data*: a SPARTN orbit names the IODE (SF101) and type SF098 0,
  the 1044 (E1, E2). IGS SSR: IDF012 is the IODE (IG8).
- *Satellites*: numbered as SPARTN numbers them, the PRN less 192: 1 to 64
  from SPARTN (SF095); 1 to 10 (PRN 193 to 202) from IGS SSR (IG1); 1 to
  15 from the 1044 (DF429, 4 bits, of which RTCM defines 1 to 10, 11 to 15
  taken as satellites 11 to 15; 0 refused), so a SPARTN correction for 16
  to 64 finds no ephemeris (2.9); 1 to 10 in the MSM.
- *Time*: GPS time to the whole second (T1, Q14).
- *Limits*: L6 has no frequency here, so IGS SSR identifiers 9, 10 and 17
  are left out (IG7) and MSM IDs 9 to 11 are never sent (2.9);
  identifiers 5, 8, 11 to 16 and 18 to 31 are reserved and left out (IG7).

**SBAS and every other system** are not supported: `Sat` names only the
five systems above, so no query reaches another, and SPARTN carries none
(Table 6.1). Their data: the IGS SSR SBAS sub-types 121 to 127 are
rejected, `unsupported_system` (IG12); RTCM 1043 (SBAS) and 1041 (NavIC)
ephemerides, which gnsswire does not decode, reach the engine as an
`unsupported` `FrameFailure` and are counted (S7); an MSM of any system is
rejected, `not_an_input` (2.2); navigation words that gnsswire-nav does not
lay out (SBAS, NavIC and the signals named above) are rejected,
`unsupported_system`.

### 2.3 Query

```text
query(sat: Sat, t_rx: GnssTime, receiver: (x, y, z) ECEF metres,
      signals: tuple[Signal, ...] = ()) -> LineOfSight
```

- `sat` names a satellite by system and number (PRN; GLONASS slot).
- `t_rx` is the reception time in GPS time. The engine computes the
  transmit time by iterating the light time from the receiver position with
  the corrected orbit, until the change is below 1 ps (at most 5
  iterations; otherwise every term but the inputs is missing,
  `not_converged`). The clock does not enter: t_tx is the GPS time of
  transmission, and τ = t_rx − t_tx is geometric. Each step evaluates the
  broadcast orbit and its correction (C1) at t_tx = t_rx − τ, taken as a
  `GnssTime` (section 4), and sets τ = (ρ + E)/c, the distance in the frame
  of reception (section 2.4, row E) [spec: Ashby 2003 section 2]. [design:
  an observation is time-tagged at reception; 1 ps is 0.3 mm of range and
  the iteration converges quadratically from the first guess t_rx − 0.075 s;
  S is left out of the light time, since its 2 cm move t_tx by under 0.1 ns,
  under 1 µm of satellite motion]
- `receiver` is the approximate receiver position. It selects the
  atmosphere area, gives the elevation for mapping functions and tides, and
  the line of sight. Atmosphere terms need it within metres to tens of
  metres; geometric terms are computed from it exactly, for the position
  given.
- `signals` lists the signals asked for. Per-signal terms (code and phase
  biases, ionosphere delays, frequencies) are given for exactly those, and
  `phase_model(s)` forms the wind-up in metres for one of them; an empty
  tuple gives every signal for which the stream has a bias. [design: the
  record states what was asked and what is missing]

The query never extrapolates beyond the validity of section 2.5: a term that
is not valid is absent from the record with a reason (section 2.4), never
silently zero.

### 2.4 The line-of-sight record

`LineOfSight` is a frozen dataclass. Every term is separate, in SI units, and
carries the sign it has in the observation equations the record is built
for (metres; `s` a signal of wavelength `λs`):

```text
P_s   = ρ + E + S − c·(dt_sat + dt_rel) + T + I_s + b_s     + R_tide + (receiver terms)
Φ_s   = ρ + E + S − c·(dt_sat + dt_rel) + T − I_s + B_s + λs·w + R_tide + λs·N + (receiver terms)
```

| Symbol | Record field | Unit | Definition and basis |
|---|---|---|---|
| | `sat`, `t_rx`, `receiver`, `t_tx` | | inputs; `t_tx` GPS time of transmission |
| `ρ` | `geometry.range` | m | distance from receiver to the corrected antenna phase centre at `t_tx`, both in ECEF of their own epoch |
| | `geometry.range_rate` | m/s | ρ̇ = u·ṙ_s, the rate of ρ for a receiver at rest in ECEF [design: free from the state the query has; a moving receiver adds −u·v_r, and its velocity is not an input] |
| | `geometry.unit`, `azimuth`, `elevation` | –, rad | receiver to satellite in the ECEF frame of reception, R_z(ω_E τ)·X − r_r, the direction the signal arrives from [spec: Ashby 2003 section 2]; WGS-84 local frame [spec: TM-23 B.3]; azimuth in [0, 2π) clockwise from north. They serve every term that needs a direction (T, R_tide) [design: the unrotated position turns the line of sight by ω_E τ, about 5 µrad, which moves the hydrostatic delay at 5° elevation by about 1 mm] |
| `E` | `modelled.earth_rotation` | m | Sagnac term: E = c·τ − ρ, c·τ = \|R_z(ω_E τ)·r_s − r_r\| [spec: Ashby 2003 section 2; SPARTN 8.16]. The exact form, not Ashby's first-order Eq. (7) [design: the remainder of Eq. (7) reaches 0.10 mm for GPS, 0.13 mm for Galileo and 0.28 mm for GEO and IGSO above 5° elevation (`tests/kernels/test_relativity.py`)] |
| `S` | `modelled.shapiro` | m | Shapiro delay, the logarithmic term of Ashby Eq. (96), 1.3 to 1.9 cm for GPS [spec: Ashby 2003 section 10, "the usual expression for the Shapiro time delay"; SPARTN 8.16]; *The gravitational delay*, below |
| `dt_sat` | `satellite.clock` | s | broadcast clock corrected: t_s = t_b − δC/c [spec: SPARTN 8.7] |
| `dt_rel` | `modelled.relativity_clock` | s | −2 r·ṙ/c² of the broadcast position and velocity [spec: SPARTN 8.7, 8.8; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1]; zero for GLONASS [spec: SPARTN 8.7] |
| | `satellite.position`, `velocity` | m, m/s | X = X_b − δX at `t_tx` (C1, C11); velocity from the broadcast orbit. Without usable orbit and clock corrections (no stream, none selected, no ephemeris of their IODE) the state is the broadcast one of the nearest valid healthy issue, `orbit` and `clock_correction` absent with the reason [design: one rule for every way corrections can be missing (3.4)]; `corrected` says which |
| | `satellite.clock_drift` | s/s | broadcast clock drift [spec: the system's ICD]; for a clock polynomial (C11) less (C1 + 2·C2·(t_tx − t0))/c [spec: the derivative of RTCM eqs. (3.12-7), (3.12-8)]. It excludes the rate of dt_rel, which needs the acceleration [design: not free from the state the query has; the difference of two queries gives it] |
| | `satellite.broadcast` | | broadcast position, clock, issue of data, ephemeris type, toe, age |
| | `orbit` | m | δO radial, along, cross as sent, or δO(`t_tx`) for a correction with rates (C11); δX in ECEF; IODE; yaw (deg or absent); time tag; age; `rate` δȮ (m/s) and `reference` t0, None for SPARTN |
| | `clock_correction` | m | δC as sent, or δC(`t_tx`) for a polynomial (C11); IODE continuity (s); URE index and value; time tag; age; `rate` C1 (m/s), `acceleration` C2 (m/s²) and `reference` t0, None for SPARTN. URE index 7 (a lower bound, "> 1.0 m") and index 0 (unknown) give the index and no value [spec: SPARTN Table 6.2 SF024] |
| `b_s` | `code_bias[s]` | m | the bias the raw pseudorange contains, which SPARTN says to subtract [spec: SPARTN 8.9] |
| `B_s` | `phase_bias[s]` | m | the same for phase [spec: SPARTN 8.9]; with fix flag (SF023), continuity (SF015) and a discontinuity count (rule C7), which changes exactly when the continuity of that bias ends [design: the query keeps no state (2.1), and a count, unlike a flag, cannot be missed between two queries] |
| `T` | `troposphere` | m | slant total and its parts: hydrostatic zenith at the receiver height (*Hydrostatic height reduction*, below), hydrostatic and wet mapping, polynomial and grid wet residuals [spec: SPARTN 8.12]; quality SF042; area; age; from the first provider of rule A9 that answers: `origin` (provider name, kind `correction` or `model`, stream or None), `sigma` (m, 1σ of the total zenith delay as the provider states it, None when it states none), `wet_zenith` (T_P + T_G for SPARTN, the model's wet zenith delay). For another provider than SPARTN `polynomial`, `grid`, `quality`, `area_id` are None, and `age` is None for a model; the slant is formed the same way [spec: SPARTN 8.12; UNB3m eq. (8)] |
| `I_s` | `ionosphere` | TECU, m | slant TEC from polynomial and grid [spec: SPARTN 8.11]; per signal the code delay +40.3·10¹⁶ TEC/f² and the phase delay of opposite sign [spec: TM-23 eq. (5.40); SPARTN 8.11 prints the factor as 40.3 × 10⁻¹⁶, which with TEC in TECU and f in Hz is 40.3 × 10¹⁶]; quality SF055; area; age; `origin` and `sigma` (TECU) as for `T`; `with_biases` True when the ionosphere came from the stream of the record's code biases, False when not (the two are then consistent only up to a constant per satellite, rule A11), None without code biases. For another provider than SPARTN, `polynomial`, `grid`, `quality`, `area_id` are None and the delays come from its slant TEC |
| `R_tide` | `modelled.solid_tide` | m | −u·Δr for the site displacement Δr of the conventional tide-free model (*Tide system*, below) [spec: TN36 7.1.1, 7.1.1.2; SPARTN 8.16]; Δr also given, at `t_rx`. The Sun, the Moon and GMST take UT1 as UTC (0.9 s at most from UT1, 0.02 mm of tide); without leap seconds the term is absent, `no_leap_seconds` [design: GPS time for UT1 turns the Earth 18 s, moving Δr by up to 0.7 mm; T4 already refuses to guess the count] |
| `w` | `modelled.windup` | cycles | satellite and receiver antenna wind-up in (−½, ½] [spec: TM-23 5.5 after Wu et al. 1993; SPARTN 8.14, 8.16]; the satellite dipoles are x and y of the SPARTN 8.14 frame turned by the yaw SF021 (zero when absent, rule C8; Q23 for the sense), the receiver dipoles east and north of the local frame [spec: TM-23 5.5]; one value per record, in cycles, the same for every signal: λs·w in metres is not a field, and `phase_model(s)` forms it with λs = c / `frequency[s]`. Continuity across epochs needs the previous value, which the query does not keep: `toolkit.windup.unwrap(previous, current)` gives it, and the VRS uses it [design: the query keeps no state (2.1)] |
| `λs` | `frequency[s]` | Hz | carrier frequency of each requested signal, λs = c / f; for GLONASS from the satellite's frequency channel [spec: GLONASS ICD 3.3.1.1]; absent, `no_frequency`, when the channel is not known |
| | `quality` | | section 2.5 |
| | `missing` | | for each absent term, why, keyed by (term, signal), the signal `None` for a term of the satellite. `Term` and `Missing` are enumerations defined once in `model/types.py`, so that the stores give these reasons, and re-exported by `record.py` [design: the model may not import the record (3.1)]. `Missing`: `no_ephemeris`, `outside_validity`, `iode_mismatch`, `ephemeris_type`, `unhealthy`, `stale`, `outside_area`, `latitude_limit`, `dnu`, `siou_mismatch`, `aiou_mismatch`, `no_corrections`, `iode_change`, `no_bias`, `no_atmosphere`, `invalid_grid_node`, `grid_mismatch`, `no_leap_seconds`, `below_elevation`, `not_converged`, `no_frequency` |

`available` is true when position and clock exist; each other term may still
be absent. `corrected` is true when the orbit and clock corrections are
applied. Two convenience methods, `pseudorange_model(s)` and
`phase_model(s)`, return the sums above without receiver terms and without
`N`; they fail if a term they need is missing. [design: every term is also
given apart, so the sums are a convenience and never the only form]

**Tide system.** `R_tide` and every position the engine outputs are
*conventional tide-free*: Δr includes the permanent part, as TN36 7.1.1.2
prescribes for the conventional model [spec: TN36 7.1.1.2]. [design: ITRF
and WGS 84 coordinates, in which SPARTN 8.2 gives every orbit and station,
are tide-free; removing the permanent part would put the receiver in a
mean-tide frame that matches neither the satellites nor the corrections]
Coordinates in a mean-tide frame differ from tide-free ones by the permanent
deformation, `toolkit.tides.permanent_tide(site)` (TN36 7.1.1.2 eqs. (14a),
(14b); about −12 cm radial at the poles and +6 cm at the equator).

**The gravitational delay.** Ashby 2003 Eq. (96) gives the coordinate-time
delay along a path of length ρ as a term Φ0·ρ/c² plus the logarithm. S is
the logarithm, in the record and in both sums. The first term, −6.97 ×
10⁻¹⁰·ρ (−1.4 to −1.8 cm for GPS, −2.5 to −2.9 cm for GEO), is neither in
the record nor added by a sum. [design: it is the rate of the coordinate
time chosen for GPS, not a delay of the signal (Ashby: "due to the choice
of coordinate time rate"); SPARTN 8.16 asks for the Shapiro delay, the name
Ashby gives to the logarithm alone; and, being the same fraction of every
range, it is held almost wholly by an orbit and clock solution whose model
leaves it out] ρ is a coordinate distance, as are the positions SPARTN
corrects; Ashby's Eq. (97), the difference between proper and coordinate
distance, is not applied [spec: Ashby 2003 section 10, "l is a coordinate
distance"]. Both S and E use the WGS-84 GM and ω_E of section 4.

**Hydrostatic height reduction (Q7).** SF043 is the area's average
hydrostatic delay "reduced to a reference height of 0 m on the ellipsoid"
[spec: SPARTN 8.12], and SPARTN gives no formula to bring it to the
receiver's height (the Niell height term changes only the mapping). The
engine multiplies it by `toolkit.troposphere.hydrostatic_height_factor`,
d_zh(h)/d_zh(0) from UNB3m eqs. (4) and (6) with T0 and β of UNB3m Table 1
and eq. (3), heights from the ellipsoid: about 0.88 at 1 km. [spec: UNB3m;
design: a delay at 0 m applied unreduced at 1 km is 0.27 m too large]
`Config.troposphere_height_reduction` (default True) turns it off. The wet
residuals T_P, T_G are not reduced: SPARTN 8.12 gives them at the user's
location.

**Antenna phase centres.** SPARTN corrections give the satellite "nominal
antenna phase center" (SPARTN 8.8), so no satellite phase-centre offset is
applied. No receiver antenna offset or variation is in the record.

### 2.5 Validity and quality

**Validity in time.** Each stored item has the time it applies from (its
time tag, or toe/tb for an ephemeris) and an age at query time `t_rx`. An
item is used only while every rule below holds; otherwise it is absent with
the reason given. Of the SPARTN items of one kind (orbit, clock, bias,
troposphere, ionosphere of a satellite), a query at `t_rx` uses the newest
whose time tag is not after `t_rx`; an item tagged later is kept and is not
used for that query, so its age is never negative. [spec: SPARTN Table 7.1
TF009, the tag is the time the data refer to; design: a query of a past
time then gets the answer it would have had at that time, whatever arrived
since, which is what replay needs (2.1)]

| Item | Valid while | Basis |
|---|---|---|
| broadcast ephemeris, GPS | half the curve fit interval on either side of toe: 4 h when the fit flag is 0; with flag 1, half the interval Table 20-XII gives for the IODC (6, 8, 14 or 26 h) | [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1, 20.3.4.4, Tables 20-XII, 20-XIII] |
| broadcast ephemeris, QZSS | `\|t − toe\| ≤ 1 h` | [spec: RTCM DF457; design: both flags take the 2 h flag 0 guarantees (Q14), toe at the midpoint as for GPS] |
| broadcast ephemeris, Galileo | `\|t − toe\| ≤ 4 h` | [design: the Galileo OS SIS ICD 2.2 gives no fit interval (5.1.1, 5.1.9.2: a validity duration without a value); 4 h covers the nominal 3 h issue with margin] |
| broadcast ephemeris, BeiDou | `\|t − toe\| ≤ 2 h` | [design: the BDS ICD gives no fit interval; ephemerides are issued hourly] |
| broadcast ephemeris, GLONASS | `\|t − tb\| ≤ 15 min + P1/2` | [design: GLONASS ICD 4.4 gives tb and its interval P1, not a validity] |
| any SPARTN item | its SIOU still current, and age < 510 s | [spec: SPARTN Table 6.2 SF005] |
| GAD area, HPAC data | AIOU matches, and no outage of 450 s or more of the stream's GAD and HPAC messages up to the query | [spec: SPARTN Table 6.2 SF068] |
| orbit, clock, bias, atmosphere | age ≤ `config.max_age` per kind (defaults: orbit 120 s, clock 30 s, bias 300 s, atmosphere 120 s) (Q5) | [design: SPARTN gives no age limit and its corrections carry no rates, so an old value is a constant extrapolation; `max_age` bounds the extrapolation of a rate (C11) alike; the defaults lie within the 510 s and 450 s rules] |
| a satellite flagged DNU | never, from the time of the flag | [spec: SPARTN Table 6.2 SF013] |

**Validity in area.** Atmosphere data are used only when the receiver is
inside the area that defines them, and not above ±85° latitude [spec:
SPARTN 8.10]. A position on the area's edge is inside [design: the 8.10
figures show only positions strictly inside or outside; the bilinear
weights are defined on the edge, and adjacent areas then leave no gap
between them]. Interpolation never crosses from one area to another [spec:
SPARTN 8.13, stated for BPAC; design: applied to HPAC grids, since the grid
of one area is defined only on its nodes].

*Area geometry (Q8).* The reference point (SF032, SF033) is the north-west
corner, and the polynomial expansion point ("centroid") the centre of the
node grid, (ref_lat − (n_lat−1)·Δφ/2, ref_lon + (n_lon−1)·Δλ/2) [spec:
SPARTN 8.10 figure]. Longitudes are compared modulo 360° and the centroid
and δλ reduced to [−180°, 180°) [design: SPARTN 8.13 says only that the
±180° discontinuity "must be considered"; modulo 360 does it for every
reference longitude]. Grid residuals run from the reference point west to
east along a row, rows north to south, node (row, col) = row·n_lon + col
[spec: SPARTN 8.10 figure, "Grid node order"]; the bilinear weights are
those of the rectangular case of 8.13, W1 = xy (NE), W2 = (1−x)y (NW), W3 =
(1−x)(1−y) (SW), W4 = x(1−y) (SE) [spec: SPARTN 8.13; 8.11 and 8.12 say
only "bi-linear", and this is the only reading of 8.13 in which each weight
is 1 at its own node]. An invalid node is Q9.

**Quality indicators** (`LineOfSight.quality` and `engine.status()`), each
with the SPARTN field it comes from:

- stream: solution ID TF010, processor ID TF011, SIOU SF005, AIOU SF068,
  datum SF009 (ITRF or provider-defined), whether the last OCB set ended
  with an end-of-set flag or is incomplete, counts of frames that failed
  and why (S7), notices (S9), stream switches (S2) and datum changes (S6);
- satellite: URE SF024 (index and metres), fix flag SF023, continuity SF015
  of the satellite and of each phase bias, IODE continuity SF022, yaw
  present;
- atmosphere: troposphere quality SF042 (m, 1σ at zenith), ionosphere
  quality SF055 (TECU, 1σ), area ID, grid or polynomial only; per term the
  `origin` (who answered, a correction or a model) and `sigma` as the
  provider states it (SPARTN: the SF042, SF055 bound, None for index 0 and
  for the last step, which is a bound below [spec: SPARTN Table 6.2 SF042,
  SF055]; UNB3m: 4.9 cm [spec: UNB3m Table 7]). Which providers were passed
  over is not recorded: `origin` says who answered, and `missing` why no
  one did;
- ages of every term.

**Configuration** (`Config`): `max_age` per kind; `elevation_cutoff`, below
which the query returns geometry only (`below_elevation`, default 0°);
`troposphere_height_reduction` (default True); the lead `stream_switch` of
rule S2 (default 5 s); the `hold_limit` of rule T5 (default 256); `vrs`, the
virtual base station's parameters (section 2.6); `satellite_stream`, the
stream of orbit, clock and biases (default None: the active stream of S2;
rule S10); `atmosphere`, the sources of rule A9 in walk order, each
`SpartnSource(stream)` (None: the satellite stream, bound by its SIOU, A1),
`ModelSource("unb3m")` or `CallerSource(name)` (default
`(SpartnSource(),)`). Nothing else. [design: every default is listed here
so a replay can reproduce a run]

### 2.6 Virtual base station

The virtual base station (VRS) turns the engine into RTCM 3 observations for
a receiver that does RTK and accepts a network's non-physical reference
station [spec: RTCM 3.1.5]. It uses only the public query, as any caller
does.

```python
station = VirtualBase(engine, vrs_config)      # observations from vrs.observations.Observations
messages = station.epoch(t, rover_position)    # gnsswire.rtcm3 dataclasses, in send order
wire = b"".join(rtcm3.encode(m) for m in messages)
```

`VirtualBase(engine, config, source=None)` joins three parts: the placement
(2.6.1, `vrs.placement`), the schedule (2.6.2, `vrs.schedule`) and an
*observation source*, an object with `observe(t, station) -> [Msm]` that
gives the MSMs of epoch `t` at `station.position` (2.6.3). Without a
`source` it is `vrs.observations.Observations(engine, config)`. [design:
the station needs nothing of the engine itself, so placement and schedule
are tested without one] The station sets DF003 and DF393 of the MSMs and
orders them; when `station.new` is true the source restarts every arc
(2.6.3).

What an RTK receiver takes from a reference station, and how each is given:

| Need | Met by | Basis |
|---|---|---|
| observations of every satellite and signal, one epoch at a time | one MSM per constellation per epoch, MSM7 by default | [spec: RTCM 3.1.2 Table 3.1-2, 3.5.16] |
| the base position | 1005 (ARP ECEF, DF141 = 1, non-physical) | [spec: RTCM Table 3.5-6, DF141] |
| antenna and receiver | 1033 with antenna `ADVNULLANTENNA`, serial empty, setup 0, receiver descriptor `ssrcore` and firmware the package version | [spec: RTCM 3.1.6 for network data; Table 3.1-2] |
| GLONASS code-phase biases | 1230 with the L1 and L2 C/A biases zero and the alignment indicator set | [spec: RTCM Table 3.1-2; SPARTN 8.9, GLONASS biases are aligned to a theoretical zero-bias receiver] |
| a station near the rover | the station is placed at the rover (2.6.1) | [spec: RTCM 3.5.12.1] |
| continuity | lock time per signal reset only on a real discontinuity (2.6.3) | [spec: RTCM DF402] |

**Receiver descriptor.** RTCM 3.5.12.5 asks that the receiver descriptor
follow the IGS 20-character definition. This receiver is in no IGS table,
so the descriptor is `ssrcore` and the firmware field the package version.
[design: no IGS name describes a computed station, and a borrowed one would
claim a receiver's GLONASS biases, which the 1230 states instead]

#### 2.6.1 Placement and movement

- The station is placed at the first rover position the caller supplies
  (`rover_position`, ECEF, or a GGA through `adapters.nmea`), rounded to
  0.1 m, at the rover's ellipsoidal height. [spec: RTCM 3.5.12.1, the
  station at the rover's approximate position; 3.5.12.2 names the GGA it is
  sent in; design: the rover's height keeps the troposphere difference
  small, and 0.1 m rounding keeps the coordinates stable across GGA noise]
  The first placement keeps `station_id`.
- It stays fixed while the rover is within `move_distance` (default 5 km)
  horizontally, in the station's local horizontal plane: a change of height
  alone never moves it. [design: a reference station's position is fixed;
  5 km keeps the station near the rover while it moves seldom; no
  specification sets the distance]
- It moves when the rover has been beyond `move_distance` for
  `move_confirm` seconds (default 10 s), counted in the epochs that carried
  a rover position: ⌈`move_confirm`·`rate`⌉ consecutive such epochs, at
  least one. An epoch without a rover position neither confirms nor breaks
  the count. The epoch that confirms is the first at the new place, which
  is the rover position of that epoch; it sends 1005 and 1033 first, then
  the observations, with every lock time reset, and the station ID changes
  (2.6.4). [spec: RTCM 3.5.16.2, position messages before observations;
  design: the confirmation avoids moving on a single bad position, and a
  time rather than a count gives the same protection at every rate; an
  absent GGA says nothing about where the rover is, and resetting on it
  would delay a move indefinitely on a link that drops sentences]
- The caller can pin the station (`pinned_position`), which disables
  moves; a pinned position is held at the 0.1 mm of DF025 to DF027, so the
  observations are computed at exactly the ARP the receiver decodes.

#### 2.6.2 Messages and rates

Per epoch, in this order [spec: RTCM 3.1.3 item 3; 3.5.16.2]:

1. 1005, when due: at the first epoch, on every move, and so that no gap
   between two sendings exceeds `station_interval` (default 10 s);
2. 1033 and, if the epoch has a GLONASS MSM, 1230, at the same times as
   1005 [spec: Table 3.1-2];
3. MSMs: GPS, GLONASS, Galileo, BeiDou, QZSS, each with the multiple
   message bit DF393 set on all but the last of the epoch. Several MSMs of
   one system are sent, in the order the source gives, when the 64 cells of
   DF396 do not hold its observations (2.6.3) [spec: RTCM 3.5.16.3.4].

- Observation rate: `rate` (default 1 Hz; 0.5 to 2 Hz allowed), epochs on
  whole multiples of the interval in GPS time. [spec: RTCM 3.1.7, Tables
  3.1-3, 3.1-4]
- Station messages every 10 s at most. [design: RTCM sets no interval for
  1005, 1033 and 1230, and bounds by 15 s the completion of the auxiliary
  station data a rover needs at start-up (3.1.7); taken by analogy, 10 s
  keeps every rate within that bound]
- MSM type: MSM7 by default; MSM4 on request (`VrsConfig.msm`). [spec: RTCM
  Table 3.1-2, MSM6 minimum and MSM7 full service for high-precision RTK;
  design: MSM7 carries the phase range rate (2.6.3) at the high resolution,
  0.14 mm of phase range; MSM4 is kept for links and receivers that take
  only the standard-resolution full message without rates]
- 1005 rather than 1006: 1006 adds only DF028, the antenna height above a
  surveyed marker, which a computed station does not have. [spec: RTCM
  Table 3.1-2, "1005 or 1006"; 3.5.3, DF028]
- Ephemerides (1019, 1020, …) are not sent. [design: the receiver decodes
  its own broadcast data]
- Header fields: DF411 clock steering 1 (the synthetic clock is GNSS time
  exactly, 2.6.3), DF412 external clock 0, DF142 single oscillator 1,
  DF417/DF418 no smoothing, IODS DF409 0 (not used). [spec: RTCM 3.1.4;
  DF409]

#### 2.6.3 Observations

`vrs.observations.Observations` is the source. At epoch `t` it queries the
engine for every satellite with a stored broadcast ephemeris
(`engine.ephemerides.satellites()`), at `t_rx = t`, at the station's ARP
(`Station.position`), for every signal with a bias, and builds the MSMs of
the configured type. [design: a satellite without an ephemeris has no
record to give]

*The ranges.* The station's receiver clock and receiver terms are zero, so
each observation is the record's sum of section 2.4:

```text
P_s = pseudorange_model(s)
Φ_s = phase_model(s) + λs·(k + N_s)
```

with `k` the whole turns of the wind-up unwrapped across epochs and `N_s`
the integer of the arc. [spec: RTCM 3.1.4, epochs aligned to GNSS time,
DF411 = 1; design: a computed receiver has no clock error to model, and a
zero clock keeps DF397 within its 0–254 ms] The solid tide enters as the
record's R_tide at the ARP, the first-order form of displacing the station
(the remainder, |Δr|²/ρ, is under 10⁻⁸ m); troposphere and ionosphere are
those of the ARP (Q22). [design: a real station's observations contain the
tide and the atmosphere of its site, so the station's contain the same
model]

*Satellites.* A satellite is observed when its record is corrected (orbit
and clock), its elevation is at least `elevation_mask` (default 5°), its
troposphere, ionosphere and tide are present, and RTCM gives it a satellite
ID: GPS 1–63, GLONASS 1–24, Galileo 1–50, QZSS 1–10 (PRN 192 + the SPARTN
number), BeiDou 1–37. Otherwise it is left out whole. [spec: RTCM Tables
3.5-90, 3.5-95, 3.5-98, 3.5-104, 3.5-107; design: an observation without a
modelled term differs from the corrected range by the whole term, metres;
a satellite the station does not report is absent from its messages]

The station uses a troposphere only when its `origin.kind` is `correction`,
and an ionosphere only when it is a correction from the stream whose code
biases the record has (`with_biases` True, rule A11). Otherwise the
satellite is left out for that epoch, as for a missing term, and
`Observations.left` gives the reason (`troposphere_not_a_correction`,
`ionosphere_not_with_biases`; for a missing term, the record's). There is
no switch for it. [design: every observation the station sends is formed
from corrections of one consistent set: a model is not a correction of this
atmosphere (UNB3m's zenith error is 4.9 cm, 1σ), and an ionosphere from
another source than the code biases differs from one consistent with them
by a constant per satellite (A11); with the default configuration every
term is a correction of the satellite stream, so nothing is left out for
this]

*Signals.* A signal is sent when the record has its code bias and its
frequency and RTCM gives it a signal ID (Tables 3.5-91, 3.5-96, 3.5-99,
3.5-105, 3.5-108), under its own RINEX code only. GPS 2D and BeiDou 1A, 1D,
1P, 5D, 5P, 6A, 7P and 8P have no ID in RTCM 10403.3 and are not sent. No
other tracking mode of a band is offered in place of the one the stream has
a bias for. [spec: rule
C9; RTCM 3.5.16.3.2 aligns the phases of one band among its signals;
design: a relabelled signal would claim a bias the stream did not give] A
signal without a code bias is not sent at all. [design: its pseudorange
would be off by the bias, metres, and the phase's integer is chosen against
the signal's own pseudorange (DF401)]

*Phase.* The phase range of a sent signal is given when the record also has
its phase bias and the wind-up, and the bias is fixed (SF023) or
`float_phase` is set (default, Q21); otherwise the cell carries the
pseudorange alone, with the invalid pattern of DF401 or DF406 and lock time
indicator 0. [spec: DF401, DF406; design: a phase without its bias is off
by a fraction of a cycle that cannot be detected]

*Integer.* At the start of an arc, N_s is the integer that puts Φ_s within
half a cycle of P_s, and it is kept for the arc [spec: RTCM DF401]. Within
an arc Φ_s − `phase_model(s)` is λs times an integer, exactly, so the
station's phases differ between satellites and signals by the model and by
whole cycles only, to the MSM resolution (0.14 mm in MSM7, 0.56 mm in
MSM4).

*Wind-up* is unwrapped per satellite with `toolkit.windup.unwrap` from the
record's fraction, and started anew when every arc of the satellite ends.

*Arcs.* The arc of a signal ends, its lock time restarts at zero and a new
N_s is chosen, when:

1. the station moves (`Station.new`, 2.6.1);
2. the signal's phase was not sent at the station's previous epoch (the
   satellite absent or below the mask, a term missing, the phase withheld);
3. its phase bias discontinuity count changes (rule C7: an SF015 gap or
   SF015 = 0, a satellite discontinuity, DNU, a stream switch S2, a datum
   change S6);
4. its phase bias fix flag SF023 changes;
5. the SIOU of the record changes;
6. Φ_s leaves the range of DF401 or DF406 about the rough range;
7. the `origin` of the record's troposphere or ionosphere changes (A9).

[spec: RTCM DF402 and 3.5.16.3.3, a loss of continuity resets the lock
time; DF401 for item 6; SPARTN Table 6.2 SF015 for item 3; design: item 4,
a bias that turns from float to fixed may move by a fraction of a cycle;
item 5, SPARTN 6.2 forbids combining parts of different SIOU, so nothing
says the bias datum continues across issues (Q4); item 7, a change of
provider steps the modelled delay by centimetres to decimetres in one
epoch, a discontinuity of the phase; with the default configuration the
origin changes only with the satellite stream, which items 3 and 5 already
end arcs for]

An arc does *not* end at a correction update of the same SIOU, at a change
of the ephemeris issue (rule C4) or at an atmosphere update: the phase
follows the model with the same N_s. [design: corrected orbit and clock
describe the same satellite on either issue (SPARTN 8.3, 8.8), and a step
of the model is common to the code and every phase of the satellite, not a
loss of the phase's continuity (DF402); a reset there would signal a
discontinuity that the phase does not have]

*Fields of a cell* [spec: RTCM 3.5.16.3, Tables 3.5-78 to 3.5-88]:

- rough range DF397, DF398: the satellite's first sent pseudorange rounded
  to 2⁻¹⁰ ms; fine pseudorange and phase range DF400/DF405, DF401/DF406
  against it, rounded to their resolution, the restoration formulas of
  3.5.16.3.1 giving back P_s and Φ_s. A pseudorange that does not fit its
  field is not sent [spec: 3.5.16.3.4].
- lock time indicator DF402 (Table 3.5-74) or DF407 (Table 3.5-76) from the
  whole milliseconds since the arc began.
- half-cycle ambiguity DF420 = 0 [spec: DF420, set only for a phase of
  unresolved polarity; the synthesised phase has none].
- CNR DF403/DF408 = 0, "not computed, or not available" [spec: DF403;
  design: nothing is measured, and an invented CNR would pretend to be a
  measurement].
- phase range rate (MSM5, MSM7) DF399 + DF404: ρ̇ − c·ḋt_sat, the rate of the
  modelled range of the satellite (`geometry.range_rate`,
  `satellite.clock_drift`), the same for every signal; DF399 the whole
  metres per second nearest. [spec: DF399, DF404, the sign of the
  derivative of the phase range; design: the rates of the remaining terms
  (dt_rel, E, S, T, I_s, w, R_tide) are not in the record; together they
  stay under a few cm/s, and the test bounds the difference from the
  derivative of the sent phase range]
- extended satellite information (MSM5, MSM7): for GLONASS the frequency
  channel DF419 = K + 7, from the engine's GLONASS ephemeris; 0 for the
  other systems. [spec: 3.5.16.5.2, Table 3.4-6; 3.5.16.4.2]
- epoch time: GPS DF004 and QZSS DF428 the milliseconds of the GPS week,
  Galileo DF248 of the GST week, BeiDou DF427 of the BDT week (GPS − 14 s),
  GLONASS DF416 and DF034 the day of week and milliseconds of the day in
  GLONASS time, with the leap seconds of the time base (rule T4). [spec:
  3.5.16.4.2, 3.5.16.5.2, 3.5.16.6.2, 3.5.16.8.2, 3.5.15.8.2; GLONASS ICD
  3.3.3]
- satellites by ascending ID and signals by ascending signal ID, cells in
  mask order; when a system's satellites times signals exceed the 64 cells
  of DF396, the satellites are split, in ID order, over consecutive MSMs
  each of which holds what fits. [spec: Table 3.5-78, X ≤ 64; 3.5.16.3.4,
  splitting satellite by satellite]

*Phase alignment* (Q18). The phases are those of the SPARTN biases' codes,
not shifted to the reference signals of RTCM Table 3.5-73. [design: SPARTN
gives each bias per RINEX code (8.9), and whatever fraction of a cycle a
signal's convention holds is the same for every satellite of the signal, so
it cancels in any difference between satellites; a shift chosen here could
only be checked against a receiver]

The MSM arithmetic (rounding to the fields' resolutions, lock time
indicators) is written in `vrs/msm.py`, not as a kernel. [design: it is the
encoding of a wire format, a few roundings per cell beside a query that
costs a hundred times more, with no arithmetic a microcontroller port of the
kernels needs]

#### 2.6.4 Identity

- Station ID DF003: `station_id` (default 0) at the first placement,
  incremented by one modulo 4096 at every move and at no other time.
  [spec: RTCM DF003, IDs reassigned "only when absolutely necessary" and
  linking the observations to the station messages; design: a new station
  position is such a necessity, since observations of the new position
  under the old ID would be linked to the old one]
- 1005: DF141 = 1 (non-physical), DF142 = 1, DF021 = 0, DF364 = 0, DF022 to
  DF024 the systems of the epoch's MSMs. 1032 (physical reference station)
  is not sent: there is none. [spec: RTCM DF021, DF364 (3.5.16.2); 3.5.12.4]
- 1230: DF421 = 1 (aligned), DF422 mask L1 C/A and L2 C/A (1010), both
  biases zero, the P-code biases absent. [spec: SPARTN 8.9; Table 6.2 SF026,
  SF028, the GLONASS biases are of 1C and 2C; RTCM DF421 to DF426]

### 2.7 Toolkit

The stateless functions are public, in `ssrcore.toolkit`, one module per
kernel. Each is a typed Python function over a C kernel; none holds state.
Each module is usable without the engine.

| Module | Contents | Specification |
|---|---|---|
| `timescale` | GPS, GST, BDT, GLONASS time and UTC; week and second of week, week rollover; calendar and day of year; SPARTN 32-bit tag (seconds since 2010-01-01) and the one resolver of the 16-bit half-day tag (rule T2), which the time base uses | IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4, 20.3.3.3.1.1, 20.3.3.3.3.1, 20.3.3.5.2.4; Galileo OS SIS ICD 2.2 5.1.2, 5.1.7, 5.1.8; BDS ICD 3.3, 5.2.4.4, 5.2.4.18–21; RTCM DF427, DF428; GLONASS ICD 3.3.3, A.3.1.3; TM-23 A.1.4; SPARTN Table 7.1 TF009 |
| `geodesy` | WGS-84 ECEF ↔ geodetic, local ENU frame, line of sight, azimuth and elevation, range and its rate, rotation about z | TM-23 3.1.4 and Appendix B (B.1 to B.3); WGS-84 a from IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV and TM-23 Table 3.2, f from TM-23 Table 3.2 |
| `orbits` | Keplerian broadcast orbit and clock (GPS, Galileo, BeiDou MEO/IGSO and GEO, QZSS); GLONASS RK4 integration and clock. Velocity: GPS and QZSS by IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; Galileo and BeiDou the derivative of the specified position [design: their ICDs give position only]. `model.ephemeris.broadcast_state` feeds them from a stored ephemeris | IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3, 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.1, 5.1.4; BDS ICD 4.3 Table 4-1 (GEO PRNs), 5.2.4.9–12; GLONASS ICD 3.3.3 (clock), A.3.1.2 with the two corrections of RTCM 3.5.13 |
| `ssr` | SPARTN field scaling with the special values, ranges and note offsets of Table 6.2 (kernel `scale`); radial-along-cross to ECEF and clock correction; the orbit correction with rates and the clock polynomial at a time, and the polynomial's rate (kernel `correction`, rule C11) | SPARTN 6.3, Table 6.2, Table 8.0, 8.7, 8.8; RTCM 3.5.13.5, 3.5.13.6; IGS SSR 4.1, 4.2 |
| `troposphere` | Niell hydrostatic (with height correction) and wet mapping; SPARTN slant delay; height reduction of the hydrostatic delay (Q7) | SPARTN 8.12; TM-23 5.4.2.2.1 for Niell (1996); UNB3m eqs. (1)–(4), (6) and TM-23 5.4.2.1 eqs. (5.62), (5.64) |
| `unb3m` | UNB3m a-priori zenith delays: sea-level meteorology by latitude and season, hydrostatic and non-hydrostatic zenith delays at a height, the stated accuracy | UNB3m eqs. (1) to (11), Tables 1, 2, 7; TM-23 5.4.2.1 for g and D_min (Q25) |
| `ionosphere` | TEC to metres per frequency (the slant TEC of an area is `area`'s); the IGS SSR spherical-harmonic VTEC: line of sight on the 6370 km sphere, pierce point and mapping, VTEC at a pierce point | SPARTN 8.11; IGS SSR 4.5.1, 8.4.1 |
| `area` | GAD extent, inclusion, grid nodes, centroid, bilinear weights; the area polynomial (equation types 0–2) and the grid residual, for troposphere and ionosphere alike [design: 8.11 and 8.12 state one polynomial and one interpolation; one kernel evaluates both, and rule A7 holds for both] | SPARTN 8.10 and its figure, 8.11, 8.12, 8.13 (rectangular weights), Table 6.24, Table 6.2 SF041, SF054 |
| `astro` | Sun and Moon positions in ECEF from truncated series (Moon about 4.3″, 3.3 km; Sun 2″), ERA, GMST, nutation; `sky`, the one function that composes them for an epoch, whose result the tide and the wind-up take as an argument (Q11) [design: the truncation is far below what tides and wind-up need] | ELP 2000-82B (Chapront-Touzé and Chapront 1983, 1988); VSOP87D (Bretagnon and Francou 1988); TN36 ch. 5 (5.15, 5.32, 5.35 extract, 5.40, 5.43) |
| `tides` | solid Earth tide displacement, step 1 and step 2 of the conventional model; the permanent deformation | TN36 7.1.1, 7.1.1.2; ch. 5 eq. (5.43) and ch. 6 section 6.2 for the step 2 arguments; ch. 1 Table 1.1 constants |
| `windup` | satellite frame of SPARTN 8.14 (orbit-based, from the position and the inertial velocity) turned by the SF021 yaw; nominal Sun-steered frame of TM-23 5.6.3; wind-up of satellite and receiver antennas and its continuity | SPARTN 8.14 and its figure, Table 6.2 SF021; TM-23 5.5 (Wu et al. 1993), 5.6.3 |
| `relativity` | eccentricity clock term, from the state vector and from broadcast elements; Shapiro delay; Earth rotation during travel | SPARTN 8.7, 8.16; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4; BDS ICD 5.2.4.9; Ashby 2003 sections 2, 5, 10 |

*Sun and Moon (Q11).* The Moon is ELP 2000-82B, the Sun VSOP87D, both
truncated (the terms kept and their origin are in
`THIRD_PARTY_NOTICES.md`), rotated to the terrestrial frame of TN36
chapter 5. `tools/astro_tables.py` generates the tables of `astro.c` from
the published files.

Conventions every toolkit module follows [design: the modules are written
apart and used together; `tests/test_toolkit.py` checks them]:

- Vectors are `ssrcore.toolkit.Vector`, (x, y, z), ECEF unless the name
  says otherwise; angles are radians, except in `area`, whose SPARTN fields
  and polynomial offsets are degrees (section 4).
- A function of a line of sight takes the satellite first, then the
  receiver (`line_of_sight(satellite, receiver)`, `shapiro`, `sagnac`,
  `windup`), and so does its kernel.
- A kernel status other than success is raised as
  `ssrcore.toolkit.KernelError`, a `ValueError` carrying `kernel` and
  `status` (the kernel's own code). An argument the Python layer refuses
  before calling a kernel (a wrong number of values) is a plain
  `ValueError`. An outcome that is valid data is a value, not an exception:
  a field's special value (`ssr.field`), a position outside an area
  (`area.contains`), a grid with an invalid node (`area.interpolate`
  returns `None`, rule A7).
- Day of year is counted from January 0.0: 1.0 at 00:00 of 1 January, as
  TM-23 5.4.2.2.1 eq. (5.70) counts it for Niell and UNB3m. 
  `timescale.day_of_year` gives the whole day in that count (1 for 1
  January), to which the caller adds the fraction of the day. `astro` and
  `tides` take GPS week and seconds instead, and count TT Julian centuries
  themselves (TN36 ch. 5).

### 2.8 Cost of a query

The queries of one epoch usually share one `t_rx` and one receiver
position, over some tens of satellites. The query's work is of three kinds,
and each is done once where its inputs allow [design: a small processor;
every cache is keyed by everything its value depends on, so results are
those of a query without caches]:

- *per receiver*: geodetic coordinates and the local frame (the last
  receiver is kept);
- *per epoch* (`t_rx`, receiver, store version): the Sun and Moon, the
  tide displacement, the area choice (A4), the zenith troposphere and the
  day of year (the last epoch is kept); the atmosphere providers keep
  their per-epoch values in the epoch's dict, `Where.cache`, under keys of
  their own (the SPARTN provider: its area per stream and issue, its zenith
  per area);
- *per satellite*: the selection of 3.4 steps 1 to 3 is kept with the
  interval of `t_rx` over which it cannot change while no input is fed: from
  the query's `t_rx` to the earliest of the next stored block time of the
  satellite, the end of each chosen block's age limit and the end of the
  ephemeris validity; a later query in that interval does no selection
  work. The light time, the satellite state and the line-of-sight terms are
  computed for every query, since they depend on `t_rx` itself.

Every `feed` that may change a store increments a version, which retires
the caches. Kernels allocate nothing (3.1). `Engine.selections` counts the
selections made, so that a test can see that a query in a held interval
made none.

*Measured* on macOS arm64 (Apple silicon), on the synthetic scene of
`tests/test_engine.py` (one GPS and one GLONASS satellite, full corrections
and atmosphere, grid included), best of 40 runs under a load average near
50: CPython 3.11, 48.2 µs per query at one epoch and 65.0 µs when every
second query starts a new epoch, most of it Python object work (`GnssTime`
arithmetic, record construction, the provider walk); the kernel calls alone
in C (`tools/query/`, section 4), Apple clang `-O2`, 1.21 µs per query at one
epoch and 2.82 µs per new epoch (`uv run python tools/portability.py time
200000`).

### 2.9 Known limitations

What the code does not do as it stands, in one place. Each is also stated
where the reader meets it, in the sections named after it. Those meant to
be lifted are listed in [ROADMAP.md](ROADMAP.md).

*Inputs not taken*

- SPARTN BPAC (SM 3-0) is refused, `not_an_input`, until a `bpac` kernel
  exists. (1, 2.2)
- An IGS SSR stream of centre-of-mass orbits is refused (IG11). (3.3)
- IGS SSR BeiDou GEO orbit corrections are left out, since C1 does not use
  the IGS GEO frame, ṙ_b + Ω̇_e × X_b; such a satellite has the broadcast
  state (IG4). (2.2.1, 3.3)
- SPARTN BeiDou ephemeris types 2 and 3 (B3I, SF097) are taken to have no
  RTCM message, so corrections naming them match nothing (E2). (2.2.1)
- A band without a frequency here cannot be named: GLONASS G3, G1a, G2a and
  QZSS L6, whose IGS SSR biases are left out (IG7); the QZSS L6 MSM IDs 9
  to 11 that `vrs.msm.SIGNAL_ID` lists are never sent. Galileo's E5a
  quasi-pilot has no `Signal`, MSM signal ID or bias. (2.2.1)
- A caller's ionosphere is slant TEC in TECU; a delay in metres at one
  frequency is not accepted. (3.6)

*Quantities not applied*

- The offsets between time scales that the navigation messages broadcast:
  GGTO, GLONASS τc and τGPS, BeiDou A0GPS and A1GPS, the UTC A0 and A1. GST,
  BDT and QZSS time become GPS time by whole seconds (T1), GLONASS time by
  the leap seconds alone (T4). (2.2.1)
- IGS SSR: the dispersive-bias and MW-consistency flags, the phase bias's
  wide-lane integer indicator and the regional/global CRS flag, carried by
  gnsswire (IG12); the yaw is evaluated at the orbit's t0, not at the query
  time (IG10). (3.3)

*Terms absent in stated cases*

- Between a new OCB SIOU and the HPAC and GAD of that SIOU, troposphere and
  ionosphere are `siou_mismatch` and the virtual base station observes
  nothing; the old SIOU's complete combination is not used meanwhile (S3).
  (3.3, 3.4)
- The atmosphere of a SPARTN stream other than the satellite stream is
  taken under that stream's current SIOU issue, so a query of a time before
  its last SIOU change finds no area (A10). (3.3)
- An IGS SSR Galileo or BeiDou IOD that names no stored issue is taken as
  sent and matches only by chance until the next correction (IG8). (2.2.1)
- A `max_age` below an IGS SSR stream's update interval ends its parts'
  validity early: the adapter keeps the IGS lower bound, the store does not
  (C11). (3.3)

*Ranges that disagree between inputs*

- Satellite 64 of GPS, Galileo and BeiDou from SPARTN (SF011, SF093, SF094)
  and of GPS and BeiDou from IGS SSR (IG1), and SPARTN QZSS 16 to 64
  (SF095), have no ephemeris: the RTCM ephemerides reach 63, the 1044 15.
  Their corrections are stored and never used. (2.2.1)
- `adapters.rtcm3.from_1044` takes DF429 11 to 15 as QZSS satellites 11 to
  15, though RTCM defines 1 to 10. (2.2.1)

## 3. Internal structure

### 3.1 Layers

```text
interface   engine.py, record.py, vrs/          feed, query, records, RTCM synthesis
model       model/, adapters/                   stores, selection, consistency, scaling of gnsswire messages
toolkit     toolkit/                            typed public wrappers, no state
kernels     kernels/ (C99) + _kernels/ (Cython) arithmetic only
```

**Layers.** Kernels are typed, stateless arithmetic; the toolkit wraps
them; the model holds the stores, the selection and the consistency rules;
the interface feeds, queries and synthesises. A kernel takes numbers and
returns numbers: no Python objects, no allocation, explicit integer widths,
double precision unless its specification says otherwise (section 4), the
platform's math library. Kernels are written so that they compile for a
microcontroller unchanged (section 4). [design: arithmetic without state or
objects is checked against a document's numbers alone and runs on a small
processor; the state and the rules stay in the model, above it]

Dependencies point down only: a kernel includes only `<stdint.h>`,
`<math.h>` and other kernel headers; the toolkit imports only `_kernels`
and other toolkit modules; the model imports the toolkit and gnsswire; the
interface imports the model. [design: *Layers*] Two further limits
[design: `_kernels` is private, and wire formats stay at the edge]:
`ssrcore._kernels` is imported by the toolkit only, never by the model or
the interface; gnsswire is imported by neither the toolkit nor the kernels.
So the model cannot import `record.py`: the selection of 3.4
(`Corrections.select`, `AtmosphereStore.area`, `model/selection.py`)
returns model types, and the engine assembles the record from them; the
enumerations both need (`Term`, `Missing`) live in `model/types.py` and
`record.py` re-exports them. `tests/test_rules.py` checks all of these
mechanically (section 5.1).

**Kernel conventions** [design: *Layers*, and for the details the reason
each gives]:

- C99; one header `kernels/include/ssrcore/<name>.h` and one source
  `kernels/src/<name>.c` per kernel.
- Names [design: C has one namespace for every kernel linked into one
  image, and the kernels are written apart; a name then says which header
  declares it]: an exported function is `ssr_<kernel>_<what>`
  (`ssr_geodesy_ecef_to_geodetic`, `ssr_windup_cycles`), a status code
  `SSR_<KERNEL>_<CONDITION>`, a private helper `static` with any name. A
  constant is named for its quantity and, when the value is a system's,
  that system (`SSR_SPEED_OF_LIGHT`, `SSR_GPS_MU`), not for the kernel that
  holds it (section 4). The binding's Python name is the C name without
  `ssr_<kernel>_` (`geodesy.ecef_to_geodetic`).
- Inputs and outputs are numbers and caller-owned fixed-size arrays
  (`double out[3]`); no allocation, no static mutable state, no I/O, no
  `errno` reliance. Errors are returned as `int32_t` status codes, values
  through out-parameters: `SSR_<KERNEL>_OK` is 0 and means the outputs hold
  the result; every other condition has its own positive code, numbered from
  1 within the kernel, and leaves the outputs it names untouched [design: one
  test, `status != 0`, for every kernel].
- `double` throughout, `int32_t`/`uint32_t`/`int64_t` for integers: no
  `int`, `long`, `short`, `char`, `signed`, `unsigned` or `float` in kernel
  code, and no `static` variable that is not `const`. A kernel that needs
  another's result takes it as an argument (the Sun vector for tides or for
  the Sun-steered satellite frame, for example) rather than calling across
  kernels, so kernels can be written and tested independently.
- Every function carries `[spec: …]` or `[design: …]` in its comment.
- The Cython binding `src/ssrcore/_kernels/<name>.pyx` declares the header
  and exposes each function with a docstring repeating the citation; it adds
  no arithmetic.

**Kernels**, one per toolkit module, same names, except where a module wraps
two: `timescale`, `geodesy`, `kepler` and `glonass` (module `orbits`),
`scale` and `correction` (module `ssr`), `troposphere`, `unb3m`,
`ionosphere`, `area`, `astro`, `tides`, `windup`, `relativity`. [design:
the kernel of SPARTN 8.7 and 8.8 is `correction`, not `ssr`, so that its
names do not read `ssr_ssr_`] A `bpac` kernel (SPARTN 8.13) is planned
([ROADMAP.md](ROADMAP.md)).

### 3.2 Stores

All stores are plain Python objects owned by the engine, updated only by
`feed`, read by `query`.

| Store | Key | Holds | Module |
|---|---|---|---|
| time base | – | leap seconds (GPS − UTC) and when they change; last explicit time; held messages | `model/timebase.py` |
| ephemerides | satellite, message, issue of data | up to 3 issues per satellite and message (current, previous, next) with receipt time | `model/ephemeris.py` |
| stream | (solution ID, processor ID) | SIOU and AIOU issues (`Issue`: value and era), last end-of-set, datum, failure and notice counts (S7, S9); one instance shared by the corrections and atmosphere stores, each noting its messages there | `model/stream.py` |
| corrections | stream, satellite | orbit, clock and bias blocks per SIOU issue for the last 510 s, each with SIOU, time tag, receipt time; DNU; satellite and phase-bias arcs | `model/corrections.py` |
| atmosphere | stream, area ID | GAD definitions per issue (SIOU and AIOU `Issues`); per area the latest troposphere of any constellation (A5) and per satellite the latest ionosphere (A6), each per issue | `model/atmosphere.py` |
| atmosphere providers | – | the ordered list of rule A9, built once from `Config.atmosphere` (`build`): the SPARTN store presented as a provider, the UNB3m model, the caller's providers; holds nothing itself | `model/providers.py` |

[design: three ephemeris issues cover a provider that switches IODE up to
60 s after the receiver has a new issue, as RTCM 3.5.13 recommends for SSR,
and a receiver that has a new issue before the provider switches]

### 3.3 Consistency rules

Each rule has an identifier used in code comments and tests.

**Stream (S)**

| | Rule | Basis |
|---|---|---|
| S1 | A stream is identified by (solution ID TF010, processor ID TF011). Data of different streams are never combined. | [spec: SPARTN 6.2, Table 7.1 TF010, TF011: a change gives no guarantee of continuity] |
| S2 | Several streams may be stored; the query uses the *active* stream: the one whose latest OCB time tag is newest, kept until another stream is newer by more than `stream_switch` (default 5 s). A switch counts as a discontinuity for every satellite (C7); each count of the new stream is raised above the old stream's, so the counts the record gives only increase. | [design: redundant processors may interleave; switching only on a clear lead avoids flapping] |
| S3 | Within a stream, parts sent in separate messages are combined only when their SIOU (SF005) is equal. SIOU is solution-wide: it applies to OCB, HPAC and GAD together, so the OCB of one SIOU never combines with the atmosphere of another (Q4); when the OCB SIOU changes before the HPAC and GAD of the new SIOU arrive, troposphere and ionosphere are `siou_mismatch` until they do, and the complete combination of the old SIOU is not used meanwhile (2.9). | [spec: SPARTN 6.2, the SIOU "is a solution-wide indicator and is not specific to any message type"; Table 6.2 SF005] |
| S4 | After a gap of 510 s or more in a stream, every SIOU is treated as new even if equal. The gap is measured between time tags; a message 510 s or more older than the stream's newest is refused as stale, since its issue cannot be compared. | [spec: Table 6.2 SF005; design: receipt times may be absent (T5)] |
| S5 | An OCB set ends at the message with end-of-set SF010 = 1, also one of a subtype without a table, whose header alone is read (`OcbHeaderInput`; its time tag, in a time scale SPARTN has not defined, is not compared: its flag ends the set open in its stream); a set also ends when an OCB message with a newer time tag arrives. Data are usable per satellite as soon as received; the set state is a quality indicator. A failed OCB frame (3.5) marks the open set of every stream incomplete, or, between sets, the next set to open; a late message of an ended set changes no set state. | [spec: SPARTN Table 6.2 SF010, 8.4 item 5; design: each satellite block is complete in itself, and waiting for the end of set would add latency and lose a whole epoch when its last message is lost; a failed frame carries no readable stream key, and a set's frames are sent together] |
| S6 | The datum flag SF009 is reported; a change within a stream is a discontinuity for every satellite. | [spec: Table 6.2 SF009, the datum is constant per solution and processor ID; design: a change breaks that] |
| S7 | Frames that fail (section 3.5) change no store; they are counted per (protocol, reason, type, subtype). | [design: a frame that cannot be read carries nothing that can be trusted; the protocol keeps numbers of different protocols apart] |
| S8 | Per block, since one subtype's set may span several messages with the same tag: a block equal to the stored one of the same time tag (the same solution received over a second transport) is a duplicate and changes nothing; a different block of the same tag replaces it and is counted as a conflict; a block older than the newest of its kind and issue is refused as superseded. Applied to OCB; HPAC and GAD keep the latest received of an epoch (A5). | [design: one solution may arrive over two links; equal content is one message, and a conflict must be visible rather than resolved silently] |
| S9 | A message that is used although a reader may meet it in more than one way is counted as a *notice* per (protocol, notice, type, subtype) beside the failures of S7: `extension`, data after the table (3.5); `dnu_table_8_2`, an OCB message whose DNU satellites gnsswire read with the 7 bits Table 8.2 sizes, because only that layout fits the payload. The DNU block is used the same in both layouts: its other fields are not read (C6). | [spec: SPARTN 8.4 item 4; Table 6.2 SF013 and Table 8.2, which disagree on the DNU block; design: the data are used as the specification allows, and the reading taken is made visible] |
| S10 | Orbit, clock, code and phase biases come from one stream: `Config.satellite_stream` when named, otherwise the active stream of S2. One selection (one stream, one SIOU issue) gives all four; nothing configures the biases' stream apart. A named stream that goes quiet is not replaced: its blocks age out by `max_age` and the satellite falls back to the broadcast state (2.4); it returns with its own continuity (C7). | [spec: SPARTN 6.2, SF005 and 8.6: a stream's biases are defined by its own clocks; RTCM 3.5.13.6 and IGS SSR 4.2, "consistent transmission of clock and code bias parameters"; design: in Φ_s only −c·dt_sat + B_s is determined by a solution, and another solution's clock differs by a per-satellite quantity that is not a whole number of cycles; leaving `satellite_stream` None gives S2's failover] |

**Time (T)**

| | Rule | Basis |
|---|---|---|
| T1 | SPARTN time tags are in the system time of the message's constellation; GAD and multi-constellation messages are in GPS time. They are converted to GPS time on input. QZSS time is GPS time to the whole second. | [spec: Table 7.1 TF009; RTCM DF428, the QZSS week begins as GPS's and, unlike DF427 for BDT, no offset is given] |
| T2 | A 16-bit tag (seconds of the half day) is resolved to the half day that puts it nearest to the message's `received` time, the reference taken in the tag's time scale, by `toolkit.timescale.spartn_resolve_half_day`; a tag exactly 6 h away goes to the earlier half day. | [spec: TF009, "from another source"; design: the receipt time is the source, and nearest is unambiguous for latencies under 6 h; a message arrives after its tag] |
| T3 | The 32-bit tag counts seconds from 2010-01-01 00:00:00 of the same time scale. For GLONASS, the epoch is 2010-01-01 00:00:00 GLONASS time (2009-12-31 21:00 UTC) and the count is the GLONASS clock reading, 86 400 s to the day, so GLONASS half days begin at 21:00 and 09:00 UTC (Q2). | [spec: TF009; GLONASS ICD 3.3.3; design: SPARTN does not say how GLONASS tags count; elapsed SI seconds would differ by the leap seconds since 2010, and only the clock reading makes the 16-bit tag equal the 32-bit tag modulo 43 200] |
| T4 | GLONASS time is UTC(SU) + 3 h; the conversion to GPS time uses the leap seconds of the time base: the count in effect at the instant converted, the new count applying from the first second of the UTC day after the leap second. Without leap seconds, GLONASS inputs are held and GLONASS queries answer `no_leap_seconds`. A `LeapSeconds` comes from NAV-TIMELS only with a valid current count from a source other than the firmware default (0) or unknown (255), or from the broadcast navigation words: GPS subframe 4 page 18 of Data ID 2 (a page of another Data ID is not read), Galileo I/NAV word 6, BeiDou D1 subframe 5 page 10 or D2 page 102 (BDT − UTC plus 14 s). | [spec: GLONASS ICD 3.3.3; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4, 20.3.3.5.1.1 and 6.4.6.4 (IRN-IS-200N-003: Data IDs other than 2 may be broadcast, and the user determines the Data ID before using LNAV data); Galileo OS SIS ICD 2.2 5.1.7; BDS ICD 5.2.4.18, 3.3; u-blox F9 HPG 1.51 3.15.26.1 and u-blox X20 HPG 2.11 3.14.22.1; design: guessing the count would make a 1 s error silently] |
| T5 | The time base's current time is the maximum of every time observed: receipt times, 32-bit SPARTN tags, a ZDA date and time (only with leap seconds), and a GGA time of day placed in the UTC day nearest the current time (only with leap seconds and a current time). It never moves backwards. Until it is known, messages with a 16-bit tag are held (at most `hold_limit`, default 256, oldest dropped and counted) and released in arrival order when it becomes known. | [design: replay of byte recordings without receipt times; holding rather than guessing keeps T2 unambiguous] |

Readings of the UTC conversions of the toolkit: Galileo OS SIS ICD 2.2
5.1.7 Eq. 21 takes ΔtLS in case a, as IS-GPS-200N with IRN-IS-200N-003 and
-004 20.3.3.5.2.4 a does. Issue 2.1 printed ΔtLSF there (its Eq. 20), which
the toolkit already read as ΔtLS, since only so does case b give a
continuous UTC; issue 2.2 corrects it (its Document Change Record). BDS ICD
5.2.4.18 gives DN's range as 6 and windows from DN + 2/3 to DN + 5/4; DN is
read as the day of week 0 to 6, the leap second at the end of day DN.

**Ephemeris (E)**

| | Rule | Basis |
|---|---|---|
| E1 | A correction applies to the broadcast ephemeris whose issue of data equals the correction's IODE field: GPS IODE (SF018), GLONASS tb (SF019, 7 bits, Q3), Galileo IODnav (SF099), BeiDou IODC = mod(toc/720, 240) for D1/D2 (SF100), QZSS IODE (SF101). | [spec: SPARTN 8.3; Table 6.2 SF018 to SF101; RTCM DF392 for the GLONASS IOD, since SPARTN does not say which value SF019 is] |
| E2 | The ephemeris must also be of the type the OCB header names (SF016, SF017, SF096, SF097, SF098): Galileo F/NAV from 1045, I/NAV from 1046; BeiDou D1 (MEO and IGSO) and D2 (GEO) from 1042; a type with no RTCM message matches nothing. | [spec: Table 6.2 SF096, SF097; RTCM Tables 3.5-110, 3.5-111; BDS ICD 5.1.1, Table 4-1] |
| E3 | Orbit and clock by the ICD algorithms SPARTN names; GLONASS by the re-calculation of GLONASS ICD A.3.1.2 with RK4 and a step of at most 30 s. | [spec: SPARTN 8.3 Table 8.1; RTCM 3.5.13] |
| E4 | Broadcast group delays (GPS and QZSS TGD; Galileo BGD E5a/E1 and E5b/E1; BeiDou TGD1, TGD2; GLONASS Δτn) are applied in no case: not with code biases of either source, not to a signal without a code bias, not to the broadcast state. They are stored with the ephemeris, apart from its clock (`KeplerEphemeris.group_delays`, `GlonassEphemeris.delta_tau`), and are in no record. A signal without a code bias has no b_s: `code_bias[s]` is absent with its reason (`no_bias`, or the selection's when there is no bias block, 3.4 step 2) and `pseudorange_model(s)` refuses it. A caller that wants a group delay reads it from the issue the record used: of `engine.ephemerides.issues(sat)`, the one whose `iod` and `message` are those of `satellite.broadcast`. (Q6) | [spec: SPARTN 8.6, the parameters are free of satellite signal biases other than those of its code and phase bias fields, and 8.9, the biases are subtracted from the raw measurements: a SPARTN code bias is the signal's whole bias against the stream's clock, and a group delay would count part of it twice; IGS SSR 4.3, the biases are absolute, against a bias datum the provider chooses, and 4.2, the provider makes clocks and code biases consistent, while a group delay is relative to the broadcast clock's datum (for GPS the L1/L2 P(Y) ionosphere-free combination, IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.2); design: without a code bias the record gives no b_s rather than one from another datum, and the broadcast clock is given as broadcast, the dual-frequency clock of IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.2] |
| E5 | An ephemeris with an unhealthy flag is stored but not used: GPS SV health not 0; GLONASS Bn MSB or ln set; Galileo signal health not 0 or data validity 1 for any signal; BeiDou SatH1 set; QZSS one of the five defined bits of DF454 set (NAV data, L1 C/A, L2C, L5, L1C; the LSb, reserved for LEX, is not read). | [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.4; GLONASS ICD 4.4; Galileo OS SIS ICD 2.2 5.1.9.3; BDS ICD 5.2.4.6; RTCM DF454] |

**Corrections (C)**

| | Rule | Basis |
|---|---|---|
| C1 | Orbit: X = X_b − δX, δX the radial, along, cross correction in the frame of the broadcast state at `t_tx`. | [spec: SPARTN 8.8] |
| C2 | Clock: t_s = t_b − δC/c, plus the relativistic term for every system but GLONASS. | [spec: SPARTN 8.7] |
| C3 | SPARTN corrections are constants: no rate, no extrapolation; ages are limited by the table of 2.5. A format with rates is evaluated by C11. | [spec: SPARTN 6.5.1; design for the limits] |
| C4 | A clock block combines with the orbit block of the same satellite and SIOU. If the clock block's time tag is later than the orbit block's by more than its IODE continuity SF022, the IODE may have changed: the clock is held until an orbit block with a time tag after the change arrives (`iode_change`). | [spec: SPARTN 6.2, Table 6.2 SF022; Q24] |
| C5 | Code and phase biases combine with orbit and clock of the same SIOU. Raw observations contain the biases, which are subtracted from them. | [spec: SPARTN 6.2; 8.9] |
| C6 | A satellite block with DNU SF013 = 1 removes every stored correction of that satellite in the stream and ends its arc; the satellite is unavailable until a block without DNU and with a later time tag arrives. | [spec: Table 6.2 SF013, "immediately stop using"] |
| C7 | Discontinuity: when a new block arrives later than the continuity time of its SF015 after the previous block of the same satellite (or phase bias), or SF015 is 0, the satellite's (or signal's) arc ends: its discontinuity count increments. A satellite discontinuity also ends each of its phase-bias arcs. The VRS ends the signal's arc (2.6.3). | [spec: SPARTN 6.2, Table 6.2 SF015; design: one signal's count then changes at every discontinuity that touches it] |
| C8 | Yaw SF021 is used when the header's yaw flag SF008 is set and the value is not 0x3F; otherwise yaw is zero. | [spec: Table 6.2 SF008, SF021; SPARTN 8.14] |
| C9 | Signals map to RINEX 3 observation codes by the bias masks of Table 6.2 (GPS C1C, C2W, C2L, C5Q; GLONASS C1C, C2C; Galileo C1C, C5Q, C7Q, C8Q, C6C; BeiDou C2I, C5P, C7I, C6I, C1P, C7P, C8P; QZSS C1C, C2L, C5Q, and the L codes alike). A signal asked for with another code of the same band (C2X for C2L) is not substituted: the bias is missing (`no_bias`). | [spec: Table 6.2 SF025–SF028, SF102–SF107; design: substituting tracking modes would hide a possible bias difference (Q10)] |
| C10 | A new bias block replaces the satellite's previous one wholly; a signal it no longer carries is not kept from an older block. | [design: SF014 says only whether a bias block is present; a block is the provider's whole statement for its epoch] |
| C11 | An orbit correction with rates δȮ and a clock correction with a polynomial (C1, C2) are evaluated at a time t from their reference time t0 (`reference`, GPS time): δO(t) = δO + δȮ·(t − t0), δC(t) = C0 + C1·(t − t0) + C2·(t − t0)², t − t0 the exact `GnssTime` difference; δO at each light-time step's `t_tx`, δC at the final `t_tx`. Without rates (SPARTN) nothing is evaluated. The item applies from its time tag (the epoch time, not t0) while its age is within `max_age` of its kind (2.5); t − t0 may be negative. The update interval ends no validity and is not an input: the adapter of such a format keeps it, with the format's other rules, to itself. | [spec: RTCM 3.5.13.5 eq. (3.12-6), 3.5.13.6 eq. (3.12-8), DF365 and DF376 for t0, 3.5.13.11, "The validity interval of SSR parameters can be longer than the SSR Update Interval"; IGS SSR 4, 4.1, 4.2; design: t is the time of transmission because the broadcast orbit and clock they correct are evaluated there (2.3), 0.07 s from t_rx, under 0.1 mm for any rate in range; `max_age` bounds the extrapolation of a rate as it bounds a constant] |

**IGS SSR and RTCM.** C11 is read from RTCM 3.5.13 and from IGS SSR. They
agree on the equations and signs, the units and resolutions, t0 (the epoch
time plus half the update interval, the epoch itself for interval 0), the
update intervals and the high-rate clock added to the polynomial. Where they
differ, for an IGS SSR source the IGS document holds, read as follows:

- *Clock sign.* IGS SSR 4.2 (as RTCM 3.5.13.6) writes t_satellite =
  t_broadcast − δC/c of satellite *times*, the system time a satellite clock
  reading stands for, t = t_sv − Δt_sv (IS-GPS-200N with IRN-IS-200N-003 and
  -004 20.3.3.3.3.1); of the clock offset the engine holds, that is Δt_sv =
  Δt_b + δC/c. SPARTN 8.7 writes the same formula of the offset itself
  (`ssr.corrected_clock`), so the IGS adapter gives the engine −C0, −C1, −C2
  and the high-rate clock negated (IG5). [design: under this reading the
  corrected clock is continuous across a change of broadcast issue]
- *Epoch time.* IDF003 counts seconds of the week in one time scale for
  every system, GPS time (GLONASS's UTC leap seconds removed), where RTCM
  gives GLONASS the seconds of the GLONASS day (DF386, converted by T4). For
  BeiDou the note "−14 s offset from BDS" is read as BDT sent and converted
  [design: the reading under which the BeiDou and GPS epochs of one instant
  of a stream agree]. An IGS epoch and t0 are then GPS time, the week
  resolved from a time near receipt as for an ephemeris (2.2).
- *Validity.* IGS SSR 4: "The validity interval of SSR parameters is at
  least the SSR Update Interval"; RTCM 3.5.13.11: it "can be longer". Both
  say an interval ends no validity (C11); the IGS lower bound is kept by the
  adapter, not by the store, so a `max_age` below a stream's update
  interval ends its parts' validity early (2.9).
- *GEO velocity.* IGS SSR 4.1 defines the orbit frame of a GEO satellite
  with ṙ_b + Ω̇_e × X_b; RTCM 3.5.13.5 and SPARTN 8.8 do not. C1 uses ṙ_b, so
  an IGS source's BeiDou GEO orbits are left out (IG4).
- *Reference point.* IGS: the antenna phase centre of one reference
  frequency per system (Table 2) or the centre of mass; RTCM DF374 = 0: the
  phase centre of the ionosphere-free combination. The record holds an
  antenna phase centre (2.4); a centre-of-mass stream is refused (IG11).
- *Clock datum.* The record's clock is the stream's as sent: the engine
  does not adjust a stream's clock datum. IGS SSR 4.2 asks that clocks and
  code biases be consistent; for a stream whose clocks are phase-defined
  (the ionosphere-free phase datum in the clock, the code biases keeping the
  code datum), pseudoranges corrected with its clocks and code biases are
  off by a constant per satellite.

**IGS SSR adapter (IG)**: `adapters.igs_ssr` turns each 4076 message into
the inputs the engine takes, so that every rule above applies unchanged.

| | Rule | Basis |
|---|---|---|
| IG1 | A stream is (SSR Solution ID, SSR Provider ID), an `IgsStreamId`, never equal to a SPARTN `StreamId`; satellites by IDF011 (GPS, GLONASS, BDS 1 to 63 and 0 = 64; Galileo 1 to 36; QZSS 1 to 10 = PRN 193 to 202); GLONASS 64 has no name and is left out. | [spec: IGS SSR Table 7 IDF008, IDF009, IDF011; design: the pair plays TF010, TF011, so S1 to S10, A10, A11 apply as they are] |
| IG2 | The OCB's issue (SIOU) is the IOD SSR; end of set is never claimed (IDF005 ends one sub-type's sequence, not the epoch), so a set ends by a newer epoch (S5); a message equal to the previous one of its sub-type and stream is a duplicate (S8). | [spec: IGS SSR Table 7 IDF005, IDF007; RTCM 3.5.13.11; design: S3, S5, S8] |
| IG3 | The epoch IDF003 is GPS time (BeiDou's BDT converted, *Epoch time* above), the week nearest receipt; t0 = epoch + ½ update interval, the epoch for interval 0; the OCB is tagged on the GPS scale (its system GPS, T1), the satellites naming theirs. | [spec: IGS SSR 4, IDF003, IDF004, IDF013, IDF019] |
| IG4 | Orbit: δO, δȮ at t0 with the engine's sign (X = X_b − δX). An orbit of the IOD last sent for the satellite goes at once; one of a new IOD (or the first) is held and sent in the block of the first clock of its epoch or later of the same IOD SSR, taking that tag, dropped (counted) if none comes within its age. A BeiDou GEO orbit is left out. | [spec: IGS SSR 4.1, 4.2 (the clock is of the broadcast clock of the orbit's IOD), 4 (one IOD per orbit update interval); design: the store pairs the newest orbit and clock and C4 rejects only a clock later than its orbit, so a new IOD sent alone would meet the previous clock, made against the old broadcast clock; *GEO velocity* above] |
| IG5 | Clock: −C0, −C1, −C2 at t0 (*Clock sign*); a high-rate clock is added to the polynomial of the same stream, satellite and IOD SSR whose interval holds its epoch, else left out; the IODE continuity is the time since the start of the orbit update interval holding the epoch, so C4 pairs the clock only with an orbit of that interval or later. | [spec: IGS SSR 4, 4.2, IDF004, IDF019 to IDF022; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; design: SF022's statement, read for IGS's consistency by interval] |
| IG6 | Biases: b_s = −IGS code bias, B_s = −IGS phase bias (IGS biases are added to the measurement); the fix flag is IDF029; a phase bias's continuity is the time since its discontinuity counter IDF031 changed (C7); a block holds the code and phase biases of the same IOD SSR within their age, and a code message leaves the phase out when a newer phase message is due, so a counter change ends the arc (C10 replaces a block wholly). | [spec: IGS SSR 4.3, Table 7 IDF025, IDF028, IDF029, IDF031; design: SF015's meaning, the block whole] |
| IG7 | Signals by IDF024, and RTCM DF380, DF382 where IGS reserves an identifier RTCM defines (GPS 9 = 2X, Galileo 11 = 8I, 12 = 8Q); every identifier of every system that neither defines, and every identifier of a band without a frequency (GLONASS G3, G1a, G2a; QZSS L6), is left out and counted. 2.2.1 lists them per system. | [spec: IGS SSR Table 7 IDF024; RTCM 3.5.13 DF380, DF382] |
| IG8 | The ephemeris an orbit corrects: Galileo I/NAV, BeiDou D2 for GEO and D1 otherwise, the one message of every other system (E2). The IOD: GPS, QZSS IODE; GLONASS tb with bit 7 clear (set: left out); Galileo the stored I/NAV issue whose IODnav has those 8 LSBs; BeiDou the stored issue with mod(toe/720, 240); without one the value is taken as sent and matched by a later correction. | [spec: IGS SSR 4.1 Table 1, IDF012; BDS ICD 5.1.1] |
| IG9 | URA IDF034 (≤ 3^class (1 + value/4) − 1 mm) to the coarsest SPARTN SF024 index whose bound is not below it, 0 unknown, 7 above; taken by the clock corrections of the same IOD SSR that follow within their age. | [spec: IGS SSR IDF034; SPARTN SF024; design: the record's URE is SPARTN's index (2.5)] |
| IG10 | Yaw and yaw rate of the latest phase bias message, evaluated at the orbit's t0, in [0°, 360°); the OCB's yaw flag is set when a block has one. | [spec: IGS SSR 4.3, IDF026, IDF027; SPARTN 8.14 (C8); design: one yaw per orbit block] |
| IG11 | A stream's satellite messages are used only once the caller has declared it APC (`IgsReference`), an input like any other; undeclared or CoM, they are refused (`not_an_input`); IM201 is taken either way. | [spec: IGS SSR 4.1, Table 2; no field says which point; design: the record holds an APC (2.4), and CoM needs satellite antenna offsets no input carries] |
| IG12 | Left out, counted in `IgsSsr.left_out`: the satellites and signals of IG1 and IG7, BeiDou GEO orbits, GLONASS IODs with bit 7, high-rate clocks without their polynomial, held orbits without a clock. An SBAS message (sub-types 121 to 127) is refused whole, `unsupported_system`, and not counted there. The dispersive-bias and MW-consistency flags, the phase bias's wide-lane integer indicator and the regional/global CRS flag are carried by gnsswire and not applied. | [design: no engine input has a place for them (2.9); SBAS is a system `Sat` does not name (2.2.1), so its messages are refused like other unsupported systems' (2.2)] |

**Atmosphere (A)**

| | Rule | Basis |
|---|---|---|
| A1 | HPAC data combine with the GAD of the same area ID and the same AIOU (SF068), and the same SIOU (S3). | [spec: SPARTN 6.2, Table 6.2 SF068] |
| A2 | After a gap of 450 s or more between a stream's GAD and HPAC messages, every AIOU is treated as new; a GAD or HPAC message 450 s or more older than the stream's newest is refused as stale. Kept by the stream store (`Streams.observe`). | [spec: Table 6.2 SF068; design: gaps between time tags, as S4] |
| A3 | Area ID is unique within a processor ID, so areas are keyed by stream and ID. | [spec: Table 6.2 SF031] |
| A4 | The receiver must be inside the area, and within ±85° latitude. When several areas contain it, the one whose centroid is nearest is used, by δφ² + (δλ cos φ)², ties to the lowest area ID, and within an area the issue with the newest data; an area without usable data is passed over. The same area serves troposphere and ionosphere of that query. | [spec: SPARTN 8.10; design: the nearest centroid gives the most interior polynomial expansion; one area per query keeps the terms consistent; the local plane suffices for areas at most 22.4° wide] |
| A5 | When several troposphere models are received for one epoch (one per constellation's HPAC), the latest received is used; a model of an older epoch does not replace a newer one, and an HPAC without a troposphere block leaves the stored one. | [spec: SPARTN 8.4 item 2] |
| A6 | Ionosphere is per satellite, from the HPAC of the satellite's constellation; a satellite keeps its latest block until the age limit, even when a later message omits it. | [spec: SPARTN 6.5.2, Table 6.19; design: each satellite block is complete in itself] |
| A7 | Grid residuals are interpolated bilinearly from the four nodes surrounding the receiver (*Area geometry*, 2.5); with no grid block, the residual is zero. If any of the four nodes holds the invalid value, the grid residual is absent and the polynomial alone is used, with the reason recorded (`invalid_grid_node`); a residual count other than the GAD's node count gives `grid_mismatch`. | [spec: SPARTN 8.11, 8.12, 8.13; design for invalid nodes: no extrapolation from three nodes (Q9)] |
| A8 | Equation types 3 to 7 (SF041, SF054) and SF040 = 3 are undefined; gnsswire refuses them, so they never reach the store. | [spec: Table 6.2 SF040, SF041, SF054] |
| A9 | Troposphere and ionosphere each come from an ordered list of providers (`Config.atmosphere`, `model/providers.py`): per term, the first provider that has the term and answers wins; a provider whose `terms` lack the term is not asked. Without an answer the term is absent with the reason of the first provider that has it, `no_atmosphere` when none has it. The two terms may come from different providers; within SPARTN both come from one area (A4). The record names the provider in each term's `origin`. The default list is one SPARTN provider bound to the satellite stream. | [design: the caller ranks the providers; the first one that has a term is the one whose absence it most needs to know, and with the default list the value and the reason are those of A1 to A8 exactly] |
| A10 | The atmosphere of a SPARTN stream B other than the satellite stream A is used under B's own current SIOU issue (A1 and A2 within B), so a query of a time before B's last SIOU change finds no area of B (2.9), and between A and B by time alone: each B item applies from its tag while its age is within `max_age.atmosphere` (2.5). B's datum is not applied. | [spec: SPARTN 6.2, Table 7.1 TF010, TF011: an SIOU is a field of one solution, so A's and B's say nothing of each other (S1); design: time is the only quantity two streams share, both converted to GPS time by the one time base (T1 to T5); a datum difference of metres moves T by under 3 µm and I by under 7·10⁻⁴ TECU at the largest gradients SPARTN can send (SF049, SF061)] |
| A11 | An ionosphere from another stream than the code biases (or from a caller) is consistent with them only up to a constant per satellite, α_s·δ with α_s = f₁²/f_s². The record says which case holds (`ionosphere.with_biases`). | [spec: SPARTN 8.6, a provider's slant ionosphere is defined in the datum of its own satellite code biases; 2.4, the signs of I_s and b_s in P_s and Φ_s; design: the record states the case] |

*IGS SSR VTEC provider.* `model.igs_vtec.IgsVtecIonosphere` presents IM201
behind the provider interface (A9), made and named by the caller
(`CallerSource`). Per query and layer: the receiver and line of sight on
the document's 6370 km sphere, the pierce point and 1/sin(E + ψ), VTEC of
the expansion at the pierce point in the sun-fixed longitude, negative VTEC
as 0, the layers summed [spec: IGS SSR 4.5.1, 8.4.1, IDF035 to IDF041].
Steps the document leaves open, each decided [design]: the Legendre
normalization is the geodesists' 4π one without the Condon–Shortley phase
(the term "fully normalized" names no other); t of the sun-fixed longitude
is the query's GPS seconds of the day; the quality IDF041 is a 1σ of VTEC,
mapped by the largest layer factor; a coefficient at −163.84 TECU makes the
layer unusable. With `stream` None it takes the VTEC of the query's
satellite stream under the same IOD SSR (A1 read for IGS, `siou_mismatch`
otherwise); a named stream's VTEC is taken by time alone (A10). The age
limit is the larger of `max_age` and the message's update interval [spec:
IGS SSR 4, validity at least the interval]. The ionosphere is that of the
code biases' stream only when it is that stream's VTEC (A11).

### 3.4 Selection for a query

1. The satellite stream (S10): `Config.satellite_stream`, or the active
   stream (S2). If none: only broadcast state, `missing` says
   `no_corrections`.
2. Latest SIOU `s` of the stream for which the satellite has both orbit and
   clock (S3, C4), within age limits (2.5). Biases are taken from the same
   `s` (C5); if absent for `s`, they are missing, not taken from another
   SIOU. Atmosphere is A9's walk, whose default provider takes the same `s`
   (A1): the selection's SIOU `Issue` is passed to the atmosphere store's
   `area`. Without a selection that provider answers with the selection's
   reason.
3. The ephemeris with the orbit block's IODE and type (E1, E2), valid at
   `t_rx` (2.5); the type is the OCB header's, kept with the orbit block.
   Without it, or without a selection, the broadcast state of 2.4 (row
   `satellite.position`). An orbit kernel that refuses its input (a
   `KernelError`) makes the satellite absent, `not_converged`.
4. Transmit time by light-time iteration (2.3); broadcast state at `t_tx`;
   corrections C1, C2 at `t_tx`, with C11.
5. Geometry, modelled terms; troposphere and ionosphere from the walk (A9).
6. Record assembled with ages, quality and reasons.

[design: step 2 picks the newest complete combination rather than the newest
message, so that a SIOU change never mixes parts]

### 3.5 Frames that end early or fail

| Case | What reaches the engine | Engine action | Basis |
|---|---|---|---|
| bytes stop inside a frame | nothing (`NeedMore` stays with the caller) | none | [spec: gnsswire `parse`] |
| CRC fails | nothing (`scan` skips it) | none | [spec: SPARTN Table 7.1 TF018] |
| encrypted, no key | `FrameFailure("encrypted")` if the caller reports it | counted (S7) | [design: visible, like any failure] |
| payload shorter than its table, undefined value | `FrameFailure("decode")` | counted; if it was a SPARTN OCB message, the open set of every stream (or, between sets, the next) is marked incomplete (S5) | [spec: SPARTN 8.4 item 5; design: one lost message must not discard the satellites already received] |
| data after the table (later fields, non-zero fill) | the message, its `extension` set (gnsswire's tolerant `decode`) | used; counted as a notice (S9) (Q20) | [spec: SPARTN 8.4 item 4, decoders "should not use payload length to indicate where they stop decoding"; RTCM 3.5.17] |
| OCB message of a subtype without a table (SM 0-5 to 0-15) | `OcbHeaderInput` (gnsswire `decode_ocb_header`, which reads the header of every subtype) | counted as unsupported (S7); its SF010 ends the open set of its stream (S5); the set is not marked incomplete (Q19) | [spec: SPARTN 8.4 item 5; design: S5] |
| OCB message whose header cannot be read, or another unknown type or subtype | `FrameFailure("unsupported")` | counted; an OCB one marks the set incomplete (S5) | [design: as for a failed frame] |
| a whole message is lost (gap in a set) | nothing | the set has no end-of-set; S5 closes it by the next time tag | [spec: SF010; design: S5] |
| UBX or NMEA message for which gnsswire has no table (NAV-PVT, GSV, …) | nothing; the framing helper counts it in `Framer.skipped` per (protocol, class × 256 + ID or NMEA address) | none | [design: a receiver port sends many messages that are not inputs (2.2); they are neither failures nor data, and counting them where they are seen keeps them visible without the engine knowing them] |
| navigation word (RXM-SFRBX) that fails its parity, CRC, BCH or Hamming code | `FrameFailure("decode", "ubx", 0x0213)` from `adapters.navdata` | counted (S7); the word joins no assembly | [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-XIV, Galileo OS SIS ICD 2.2 5.1.9.4, BDS ICD 5.1.3, GLONASS ICD 4.7, through gnsswire-nav; design: the UBX frame passed its checksum and its content cannot be read] |
| UBX or NMEA message with a table that its payload does not fit | `FrameFailure("decode", protocol)`, with class × 256 + ID as the type for UBX | counted (S7) | [design: as for SPARTN] |

Every `FrameFailure` names its protocol (`spartn`, `rtcm3`, `ubx`, `nmea`),
and S7 counts per (protocol, kind, type, subtype); only a SPARTN OCB failure
touches the set state of S5. [design: an RTCM 1005 and a SPARTN type 0 or a
UBX class 0x01 would otherwise share a key]

### 3.6 Atmosphere providers

Rules A9 to A11 and S10 (3.3) are the rules; this section is the interface
they are built on (`model/providers.py`). [design: one seam for every source
of atmosphere, so that a second correction format, a climatology or a
caller's own map is one more entry in a list, and no SPARTN number changes
with the default list]

**Interface.** A provider is any object with these members
(`AtmosphereProvider`, a runtime-checkable `Protocol`):

| Member | Type | Meaning |
|---|---|---|
| `name` | `str` | unique within the list; recorded as the term's origin |
| `kind` | `Kind` | `CORRECTION` (a measurement of this atmosphere) or `MODEL` (a climatology) |
| `terms` | `frozenset[Term]` | the subset of {`Term.TROPOSPHERE`, `Term.IONOSPHERE`} the provider can ever give |
| `troposphere(where)` | `Zenith \| Absent` | zenith delays at the receiver |
| `ionosphere(where, sat, look)` | `SlantTec \| Absent` | slant TEC of one satellite |

- `Where`: `t` (`t_rx`), `latitude`, `longitude` (rad) and `latitude_deg`,
  `longitude_deg` (the engine's degrees, so the SPARTN store sees exactly
  its own values, section 4), `height` (m, ellipsoidal), `doy` (2.7),
  `binding` (`Binding(stream, issue)`, the satellite corrections' stream
  and SIOU issue of 3.4 step 2, None without a selection), `unbound` (the
  selection's `Missing` reason when `binding` is None) and `cache` (the
  epoch's dict of 2.8, in which a provider keys its own entries).
- `Look`: `elevation`, `azimuth` (rad) as the record's geometry has them,
  for a provider evaluated along the path (a global map at a pierce
  point); SPARTN's polynomial is per satellite at the receiver and ignores
  it (SPARTN 8.11).
- `Zenith`: `hydrostatic`, `wet` (m, zenith, at the receiver height),
  `origin`, `time`, `age` (None for a model), `sigma` (m, 1σ of the total
  zenith delay, None when not stated or bounded below only), `spartn` (the
  store's `TroposphereZenith` for the SPARTN provider). `SlantTec`: `sat`,
  `tec` (TECU), `origin`, `time`, `age`, `sigma` (TECU), `spartn`. The
  engine turns TECU into metres per signal (SPARTN 8.11).
- `Origin(provider, kind, stream)`, `stream` the SPARTN or IGS stream a
  value came from, None for any other provider. `Absent(reason)` as the
  stores give it.
- *Never and not now.* A term a provider will never have is absent from
  its `terms`, and the walk neither asks it nor takes its reason. A term it
  has but not for this query is `Absent(reason)`: not yet (`no_atmosphere`,
  `stale`, `siou_mismatch`, `aiou_mismatch`) or not here (`outside_area`,
  `latitude_limit`).
- The walk (A9): `walk_troposphere(providers, where)` and
  `walk_ionosphere(providers, where, sat, look)`. `check(providers)`
  refuses a non-provider, a `kind` that is not a `Kind`, a term that is not
  atmosphere and a name used twice. `build(config, store, callers)` makes
  the list from `Config.atmosphere` once, at `Engine(config, providers)`;
  every `CallerSource` must have its object, carrying the same `name`, and
  every object its `CallerSource`.

**The SPARTN provider**, `SpartnAtmosphere(store, stream=None,
height_reduction=True, name=None)`, presents the atmosphere store without
changing it. With `stream` None, or naming the satellite stream, it takes
`where.binding` (A1), and without one answers `Absent(where.unbound)`; a
named other stream is taken under its own current SIOU issue (A10). The
area choice (A4), the A1 and A2 checks, the age limit and A5 to A8 stay in
the store. The area is cached in `where.cache` per (store, stream, issue)
and the zenith per area, so one area serves both terms of an epoch (A4).
The hydrostatic delay is reduced to the receiver height (Q7) here, so that
every `Zenith` is at the receiver height; the reduction is computed before
the elevation test, so a receiver above UNB3m's lapse-rate profile (some 45
to 60 km) is refused for every satellite. `sigma` is SF042 (troposphere)
and SF055 (ionosphere), None for index 0 and for the last step, which is a
bound below [spec: SPARTN Table 6.2 SF042, SF055].

**The a-priori troposphere**, `Unb3mTroposphere(undulation=None,
name="unb3m")`, made by `ModelSource("unb3m")`: `kind` MODEL, `terms`
{TROPOSPHERE}, the zenith delays of `toolkit.unb3m` at the receiver height,
`sigma` 0.049 m, `time` and `age` None; a height beyond the profile is
`Absent(not_converged)`. UNB3m's height is orthometric: `undulation`, a
function of latitude and longitude (rad) giving the geoid height N (m),
turns the ellipsoidal height into it; without one the ellipsoidal height is
used [design: no geoid model is part of ssrcore; |N| up to about 100 m
moves the zenith delay by about 3 cm, inside the model's 4.9 cm]. [spec:
UNB3m eqs. (1) to (11), Tables 1, 2, 7: total zenith delay σ 4.9 cm against
ray tracing at North American stations; Q25 for the southern hemisphere]

**What the record says of each term's source.** Per term the `origin` (who
answered, a correction or a model) and `sigma`, and for the ionosphere
`with_biases` (A11); a term's `polynomial`, `grid`, `quality` and `area_id`
are SPARTN's and None for another provider. A change of a term's `origin`
between two epochs is a step of the modelled delay by the difference of the
two providers' values; the virtual base station ends the satellite's arcs at
one (2.6.3).

**Left out.**

- Other formats' atmosphere: the SPARTN BPAC grid and RTCM network
  residuals are each one more provider, to be built with their adapter (the
  IGS SSR VTEC is one, 3.3).
- A caller's ionosphere as a delay in metres at one frequency: a provider
  gives slant TEC in TECU (2.9).
- Centre-of-mass products: they need satellite antenna offsets and
  attitude, which no input carries for every constellation, while SPARTN
  8.8 and RTCM DF374 = 0 give the antenna phase centre the record holds. A
  different reference point would be a different record row, not a field.
- A merging of satellite corrections from several sources (the best per
  satellite, averages, the orbit of one and the clock of another): S10
  forbids it, since it would mix clock datums between satellites. Failover
  between whole streams stays S2.

## 4. Numerics and portability

**The 0.1 mm target.** No numerical choice made here (a time
representation, an iteration's tolerance, a truncated series, a first-order
form, a compiler flag) may move a term of the record by 0.1 mm of range or
more; each such choice states its effect where it is made. [design: an
order of magnitude below the millimetres of the corrections themselves, so
that the arithmetic adds nothing a user of the record could see]

**Double precision.** Kernels compute in IEEE 754 double precision unless
the specification they implement says otherwise; none does (Q17, below).
[design: binary32 holds neither an orbit radius nor a second of week to the
0.1 mm target (Q17)]

- Units SI, angles in radians inside kernels; SPARTN's degrees appear only
  where a SPARTN formula uses them: the area computations of 8.10 to 8.13,
  whose fields (SF032 to SF037) and polynomial offsets δφ, δλ (8.11, 8.12)
  are in degrees, so positions given to the `area` kernel and its
  polynomial are in degrees too [spec: SPARTN 8.10, 8.11, Table 6.2;
  design: converting the 0.1° area fields to radians would only round
  them, and an edge test must compare exactly].
- Constants: one definition of each constant of each system, in one kernel
  header, cited there; every other kernel includes that header, and the
  toolkit and model re-export the binding's value rather than retype it.
  [design: one definition cannot differ between kernels]

  | constant | value | header | basis |
  |---|---|---|---|
  | `SSR_SPEED_OF_LIGHT` | 299 792 458 m/s | `correction.h` | [spec: SPARTN Table 8.0] |
  | `SSR_WGS84_A`, `SSR_WGS84_F` | 6 378 137 m, 1/298.257223563 | `geodesy.h` | [spec: TM-23 Table 3.2] |
  | `SSR_WGS84_GM` | 3.986004418 × 10¹⁴ m³/s² | `geodesy.h` | [spec: TM-23 Table 3.2; Ashby 2003 Eq. (13)] |
  | `SSR_WGS84_OMEGA_E` | 7.2921151467 × 10⁻⁵ rad/s | `geodesy.h` | [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; Ashby 2003 section 2] |
  | `SSR_PI` | the double nearest π | `geodesy.h` | [design: the scale of the angles libm's `sin`, `cos` and `atan2` take and return; also the π of `astro` and `tides`, whose series fix no value of their own] |

  The Earth-fixed computations that need GM or ω_E (relativity, Sagnac, the
  inertial velocity of the SPARTN satellite frame) use these for every
  system [design: SPARTN 8.2 gives every orbit in one Earth-fixed frame;
  the ICDs' ω_E differ by 2 × 10⁻⁸ relative, under 1 µm in any of these
  terms]. The orbit kernels (`kepler`, `glonass`) define their system's μ,
  ω_E and the π of the ICD curve fit (3.1415926535898, IS-GPS-200N with
  IRN-IS-200N-003 and -004 20.3.3.4.3.2) in their own headers, named for the
  system (`SSR_GPS_MU`), because the ICD algorithm is specified with those
  values [spec: each ICD].
- Broadcast angles stay in semicircles (and angular rates in semicircles
  per second) in `KeplerEphemeris`, as the ICDs broadcast them [spec:
  IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-III]; the orbit kernel
  turns them into radians with the π its ICD fixes for that purpose [spec:
  IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2].
- Time. The spacing of doubles just below 604 800 is 2⁻³³ s (1.16 × 10⁻¹⁰
  s), 35 mm of light travel, so a single double second of week cannot hold
  a time whose difference is multiplied by c to 0.1 mm. Two representations
  follow from where a difference is multiplied by what:
  - **Engine times** (`model.types.GnssTime`): full week, whole second of
    week (an integer) and the fraction of the second (a double in [0, 1)).
    `t1 − t2` differences the weeks and whole seconds exactly in integers
    and the fractions in a double, so it is exact to about 10⁻¹⁶ s for any
    two times; `t + dt` splits `dt` into whole and fractional seconds. Every
    difference the engine multiplies by c is taken here: the light time
    t_rx − t_tx, the ages of corrections and the clock arguments t − toc.
    [design: the light time of 0.07 s to a satellite near the end of a week
    is then kept to 10⁻¹⁶ s, 0.03 µm of range; `tests/test_types.py` shows
    the resolution achieved and that a single double loses more than 1 mm]
  - **Kernel times**: integer week (`int32_t`) and seconds of week
    (`double`), the time of evaluation and the reference time (toe, toc, tb)
    both in full weeks of the same scale (GPS time for GPS, QZSS and, after
    T4, GLONASS; GST for Galileo; BDT for BeiDou), so t − toe is the actual
    difference without the half-week rule of IS-GPS-200N with
    IRN-IS-200N-003 and -004 20.3.3.4.3. The kernels multiply time
    differences only by rates: an orbit by velocities of at most 4 km/s, a
    clock by drifts of at most 10⁻⁹. At a spacing of 2⁻³⁴ s on the
    difference (the kernels difference weeks and seconds before they add),
    that is 0.2 µm of position. `GnssTime.seconds` gives the kernels their
    double. [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1,
    20.3.3.4.3; design: 3.1, *Layers*, kernels take numbers only, and the
    ICD algorithms are in seconds of week]
- Floating-point reproducibility on the hosts: C99, no contraction of
  multiply-add, no fast-math (`setup.py`); SSE2 on x86-64. The remaining
  difference is the last-ulp behaviour of `sin`, `cos`, `atan2` and `exp`
  between math libraries, orders of magnitude below 0.1 mm. [design: bit
  equality across libms cannot be had without shipping a math library]

**Microcontrollers.** The kernels compile unchanged for a microcontroller:

- `tests/test_portability.py` compiles every kernel with `-std=c99
  -pedantic -Wall -Wextra -Wconversion -Wdouble-promotion -Werror
  -ffp-contract=off -fno-fast-math` for the host, Cortex-M4F
  (`arm-none-eabi-gcc -mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard`),
  ESP32-S3, ESP32, ESP32-C3 and ESP32-C6 (Espressif GCC). A target whose
  compiler is absent is skipped by name.
- It checks that a kernel opens only `stdint.h`, `math.h` and its own
  headers, holds no mutable data, and needs only libm, the software double
  routines and `memcpy` (a structure copy).
- `tools/query/` is one query's kernel calls in plain C, in the engine's
  order, from the fixed inputs of `tests/test_engine.py`. Its results are
  the engine's bit for bit. A counting build (C++, `count.hpp`) counts its
  double operations by kind and kernel.
- `tools/portability.py` reports the sizes, the routines, the host time and
  the counts, and builds the ESP32-S3 archive (`uv run python
  tools/portability.py report | count | archive DIR`; `uv run pytest
  tests/test_portability.py -rP` prints the same). [design: the cost on a
  device is given as operations, not a guessed time]

*Double precision on every target (Q17).* Every kernel stays in double
precision. [design: *Double precision*, above: double unless the
specification says otherwise; spec: none of SPARTN 8, TM-23, IS-GPS-200N
with IRN-IS-200N-003 and -004, the GLONASS ICD or the IERS Conventions
allows less] On a single-precision FPU (ESP32-S3, Cortex-M4F) or none
(ESP32-C3, C6), double arithmetic runs as the compiler's software routines.
From the magnitudes (no kernel was run in binary32), binary32's relative
step of 1.2×10⁻⁷ is 2 m at an orbit radius, 0.5 m at a receiver's, 1/16 s at
a second of week near 604 800, so kepler, glonass, correction, geodesy,
timescale, astro (whose time argument would resolve 94 s), relativity's
clock term (its products of 10¹¹ m²/s cancel to 10⁹) and area's inclusion
cannot use it. What could (troposphere, unb3m, area's polynomial and grid,
ionosphere, windup, tides given the Sun and Moon in double, scale) is 415 of
a query's 1 590 double operations and 1 670 of the 8 565 of a new epoch, so
binary32 would shorten a query by a quarter at most, at the price of results
that are no longer the host's bit for bit.

*Measured* with Espressif GCC 14.2 (esp-14.2.0_20241119): compiled and
counted only; nothing has run on a device. Code and constant bytes per
kernel for the ESP32-S3 (`xtensa-esp32s3-elf-gcc`, the strict flags,
`-mlongcalls`). Mutable data is 0 for every kernel on every target.

| kernel | `-O2` code | `-O2` const | `-Os` code | `-Os` const |
|---|---:|---:|---:|---:|
| area | 2 544 | 0 | 2 280 | 0 |
| astro | 7 608 | 8 560 | 6 472 | 8 800 |
| correction | 1 564 | 0 | 1 556 | 0 |
| geodesy | 3 488 | 0 | 3 362 | 0 |
| glonass | 1 840 | 0 | 1 768 | 0 |
| ionosphere | 3 641 | 0 | 3 328 | 0 |
| kepler | 4 662 | 0 | 4 462 | 0 |
| relativity | 1 712 | 0 | 1 664 | 0 |
| scale | 398 | 2 656 | 390 | 2 656 |
| tides | 8 483 | 640 | 7 676 | 640 |
| timescale | 5 320 | 56 | 3 414 | 56 |
| troposphere | 2 392 | 520 | 2 360 | 520 |
| unb3m | 1 896 | 240 | 1 880 | 240 |
| windup | 4 343 | 0 | 2 403 | 0 |
| **total** | **49 891** | **12 672** | **43 015** | **12 912** |

Totals for the other targets, code + const at `-O2` and then `-Os`:
- ESP32 (LX6): 49 905 + 12 672 and 43 000 + 12 912.
- ESP32-C3 and C6 (RV32IMC, RV32IMAC, the same figures): 61 210 + 14 728
  and 52 206 + 14 840.

Every target needs the same routines from outside the kernels:
- Math library (14): asin, atan2, ceil, cos, exp, floor, fmod, hypot,
  log, pow, round, sin, sqrt, tan.
- Software double (15): `__adddf3 __subdf3 __muldf3 __divdf3 __eqdf2
  __gedf2 __gtdf2 __ledf2 __ltdf2 __unorddf2 __fixdfsi __fixunsdfsi
  __floatsidf __floatunsidf __floatdidf`.
- Other: `__divdi3` and `__moddi3` (timescale's 64-bit division), and
  `memcpy` (structure copies in correction and glonass).

The ESP32-S3 query archive (`libssrcore_query.a`: the kernels with
`tools/query/query.c` and `scene.c`, `-O2`; the objects before the link
relaxes their long calls):
- Objects: 54 495 code and 12 952 constant bytes.
- With newlib's libm and libgcc linked in (`ld -r`): 59 338 and 13 666.
  Of that, the math and run-time routines alone are 20 726 and 714.
- `__errno` and `memcpy` are left to the firmware's C library.

Double-precision operations, counted by the C++ counting build on the host.
A conversion of an integer literal is counted, so `conv` is an upper bound.

| part | add | mul | div | sqrt | cmp | neg | conv | sin | cos | atan2 | other |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| per query, epoch cached | 542 | 667 | 113 | 37 | 73 | 12 | 45 | 47 | 29 | 6 | fabs 7, floor 11, log 1 |
| per new epoch | 2 574 | 3 322 | 101 | 27 | 37 | 31 | 1 889 | 310 | 243 | 12 | fmod 7, floor 8, pow 2, fabs 2 |
| per receiver | 35 | 66 | 4 | 6 | 2 | 3 | 0 | 9 | 9 | 4 | fabs 2, isfinite 1 |
| per ephemeris | 0 | 7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | |

Where they go:
- Per query: kepler 514 operations, relativity 273, correction 223,
  windup 164, troposphere 150, geodesy 102, area 89, the sequence itself
  63, ionosphere 12.
- Per epoch: astro 6 877 of 8 565, tides 1 455, area 140, troposphere 75.

A query makes 3 light-time iterations.

## 5. Testing

`uv run pytest` is the gate; it runs the public tests of 5.1. Releases are also checked against a private regression suite on
recorded data, which is not part of the repository.

### 5.1 Public tests

**Public tests** use only:

- worked examples from the specifications (for example the RTCM 1005
  example of RTCM 4.2; the conversion examples of Table 6.2 ranges);
- synthetic messages built with gnsswire's `encode` from chosen values, so a
  test states the physical input and checks the store, the record or the
  RTCM output;
- public data (broadcast ephemerides and precise orbits from public
  archives, listed with their terms in `THIRD_PARTY_NOTICES.md`), and
  published numerical values (for example the Niell coefficients of TM-23
  5.4.2.2.1).

Kernels are tested through their bindings: each
`tests/kernels/test_<name>.py` checks every function against the document's
numbers. Consistency rules are tested by rule identifier
(`test_rule_C4_...`). The VRS is tested by encoding, decoding and
re-deriving ranges. Tests that read local recordings of public streams
(`tests/test_igs_ssr_recording.py`) and cross-compilers that are absent
skip with the reason. No capture of a paid correction service and no key
is committed, and recordings of public streams are read from a local path,
not committed. [design: anyone can run the gate and reproduce what it
checks, and nothing whose terms forbid it is published]

**Specifications are not committed.** They are staged in `refs/`, which git
ignores but for `refs/README.md`, the list of them with the short names the
citations use. A test that needs a staged document skips without it.
[design: several are copyrighted]

The third-party tables of `astro.c` are listed with their origin and terms
in `THIRD_PARTY_NOTICES.md`. Their generator is `tools/astro_tables.py`,
and `tests/test_astro_tables.py` regenerates them from `refs/` when the
source files are staged, skipping with the reason otherwise.

*Live-service facts.* Behaviours the specification allows and a service may
show are allowed for, not made rules; each has a synthetic test of the store
that meets it:

- orbit blocks of zeros, with the whole line-of-sight correction (computed for
  the user's position) in the clock field (`test_corrections`,
  `test_live_all_zero_orbit_…`);
- areas with constant terms only (`test_atmosphere`,
  `test_live_area_with_constant_terms_only`);
- corrections that move to a new ephemeris issue some time after it first
  appears, while some constellation changes issue every few minutes: the
  ephemeris store holds the issue in use and the newer ones (3.2) and
  chooses by issue number (`test_ephemeris`,
  `test_live_issue_every_ten_minutes_…`);
- a new connection that replays recent messages: S8 and the older-block
  rule leave the stores unchanged (`test_corrections`,
  `test_live_reconnection_replays_the_last_30_s`; A5 for the atmosphere).

**Specification-driven.** Every function either implements a public
document, cited where it does so as `[spec: DOCUMENT section]`, in the
clearest and most efficient form that document allows, or is a design
choice, marked `[design: reason]`. There is no third kind: a `REVIEW`
marker, a basis still to be found, fails the gate. [design: every line can
be traced to the document it implements or to the reason it was chosen]

A rules test (`tests/test_rules.py`) checks mechanically that every C
function and every public Python function has a `[spec:` or `[design:`
marker and that no `REVIEW` marker, a basis still to be found, remains;
that kernels include no header outside the allowed list and contain no
forbidden construct (3.1); that no module imports upward across the layers;
and that the core opens no file or socket and reads no clock.

### 5.2 Cross-checks

Two tools check the engine against outside data that cannot be committed;
each takes its inputs by path.

- *Tide.* `tools/tide_check.py` builds the complete ELP 2000-82B and
  VSOP87D series with `tools/astro_tables.py` (no term truncated) and
  evaluates the solid Earth tide apart from the `tides` kernel, by an
  equinox-based rotation where the engine takes a CIO-based one, over a
  grid of sites and times, and on the published test cases of TN36's
  software. It is rerun after any change to `astro.c` or the tide kernel.
  [design: an evaluation that shares neither the truncated series nor the
  kernel checks both]
- *Orbit and clock products.* `tools/sp3_compare.py` replays a recorded IGS
  SSR stream and its broadcast ephemerides and compares the record's
  corrected and broadcast satellite states with SP3 products, the products'
  centre of mass moved to the phase centre of IGS SSR Table 2 by an ANTEX
  file. [design: a tool, not a test, because its inputs are recordings and
  products that are not committed]

## 6. Open questions

Points the specifications leave open, each with what the code does until
data or a document settles it.

- **Q6, broadcast group delays.** E4 applies no broadcast group delay in
  any case. With SPARTN code biases this is what SPARTN 8.6 implies without
  saying it; only a statement from SPARTN's authors settles it. With IGS
  SSR code biases the bias datum is the provider's (IGS SSR 4.2, 4.3).
  Without a code bias the record has no b_s, and the group delay is the
  caller's to read from the ephemeris (E4).
- **Q9, grid nodes holding the invalid value.** Rule A7 drops the grid
  residual when any of the four surrounding nodes is invalid, also one of
  zero weight. Open: whether a node of zero weight should count.
- **Q10, tracking modes.** SPARTN biases are per RINEX code; whether C2L
  and C2X (or L2L and L2X) share a bias is not stated. Meanwhile a code is
  never substituted for another (C9).
- **Q12, the Niell and Wu originals.** SPARTN cites Wu et al. (1993) and
  Niell (1996); the kernels implement TM-23 5.5 and 5.4.2.2.1, which
  reproduce their formulae. Open: the southern-hemisphere phase of Niell,
  which TM-23 does not state.
- **Q14, QZSS.** SPARTN names IS-QZSS-PNT-003. Meanwhile QZSS uses the
  constants of IS-GPS-200N with IRN-IS-200N-003 and -004, QZSS time is GPS
  time (T1), and the fit interval behind DF457 = 1 ("greater than 2 hours")
  is taken as 2 h (2.5).
- **Q18, phase alignment in MSM.** The station's phases are not shifted to
  the reference signals of RTCM Table 3.5-73 (2.6.3). Open: whether a
  receiver that pairs the station's signal with another tracking mode of the
  band sees the fraction.
- **Q21, float phase biases.** SPARTN does not say whether a phase bias
  whose SF023 flag is float may be used. The record gives every phase bias
  with its flag; the VRS emits the phase for both unless
  `VrsConfig.float_phase` is False (default True).
- **Q22, atmosphere position for the VRS.** SPARTN 8.10 to 8.12 evaluate
  the atmosphere at "the user" position; for a virtual station that may be
  the station or the rover near it. The station position is used, so a
  station outside the area has no atmosphere even when the rover is inside
  it.
- **Q23, sense of the yaw rotation.** SPARTN SF021 and 8.14 give the yaw
  angle but not the sense in which it turns the 8.14 frame about e_z. The
  code turns it right-handed, x' = cos ψ e_x + sin ψ e_y [design: the usual
  sense of a rotation angle]; the other sense is the sign of one sine in
  `ssr_windup_frame_spartn`.
- **Q24, IODE continuity when the orbit is newer than the clock.** SPARTN
  6.2 speaks of SF022 only looking back from the clock block. The code
  combines the pair whenever the clock is not later than the orbit by more
  than SF022 (C4); the other reading would also hold the clock when a newer
  orbit names another IODE than the one the clock was last combined with.
- **Q25, UNB3m south of the equator, and its height datum.** UNB3m was
  validated in North America only and its tables give no southern season;
  the provider shifts the season by half a year south of the equator (D_min
  211), as TM-23 eq. (5.62) does for UNB3 and the troposphere kernel's
  temperature row already does. Its height is orthometric; without the
  caller's undulation the ellipsoidal height is used (3.6). Settled only by
  southern-hemisphere data and a geoid model.
