# Roadmap

ssrcore takes SPARTN or IGS SSR corrections and broadcast ephemerides and
gives, for any satellite, time and position, the satellite's state and every
correction term of the line of sight, and the RTCM 3 observations of a
virtual base station. The aim is to take more of what those inputs define
and to lift the limitations of
[DESIGN.md section 2.9](DESIGN.md#29-known-limitations) listed below.
An item leaves this page when the code implements it.

## Next

- **SPARTN BPAC (SM 3-0)**: a `bpac` kernel (SPARTN 8.13) and its adapter,
  as an atmosphere provider; refused today.
- **IGS SSR BeiDou GEO orbit corrections**: C1 in the IGS GEO frame,
  ṙ_b + Ω̇_e × X_b (IG4); left out today.
- **IGS SSR yaw rate**: the yaw at the query time rather than at the
  orbit's t0 (IG10).
- **IGS SSR flags**: the dispersive-bias and MW-consistency flags and the
  wide-lane integer indicator carried with the corrections (IG12).
- **IGS SSR validity**: a `max_age` below a stream's update interval no
  longer ends its parts' validity early (C11).
- **IGS SSR issue received after its correction**: a Galileo or BeiDou IOD
  matched to its issue when that issue arrives (IG8).
- **Another SPARTN stream's atmosphere at a past time**: the newest issue
  with an area at the query time, not only the current one (A10).
- **Satellite numbers the inputs disagree on**: QZSS 1044 IDs 11 to 15,
  which RTCM does not define; satellite 64 (and SPARTN QZSS 16 to 64),
  whose corrections no ephemeris can match; QZSS L6 MSM IDs, which no
  signal can use.

## Later

- **Time offsets from navigation messages**: GPS − UTC, GGTO, BDT − UTC,
  BeiDou A0GPS, A1GPS and GLONASS τc, τGPS as an input, and whether the
  engine applies them; BeiDou's waits on the codecs.
- **SPARTN availability across an SIOU change**: the old SIOU's complete
  combination while each part is within its age (S3).
- **IGS SSR centre-of-mass streams**: need an input of satellite antenna
  offsets (IG11).
- **Galileo E5a quasi-pilot** (OS SIS ICD 2.2): waits on a signal ID and a
  bias identifier that name it.
- **A caller's ionosphere in metres** at one frequency, beside slant TEC.
- **SPARTN BeiDou ephemeris types 2 and 3 (B3I)**: whether they may use
  RTCM 1042 waits on BDS-SIS-ICD-B3I.
- **QZSS by IS-QZSS-PNT-003** in place of the IS-GPS-200 constants and
  fit interval (Q14).
- **IGS SSR identifiers** that v1.00 reserves, once a later version defines
  them.

## Not planned

- RTCM SSR messages 1057 to 1068 and 1240 to 1270: SPARTN and IGS SSR
  (RTCM 4076) are the state-space inputs.
- SBAS, NavIC and other systems: SPARTN covers GPS, GLONASS, Galileo,
  BeiDou and QZSS (Table 6.1).
- SPARTN dynamic key and authentication messages (SM 4): decryption is
  done before messages reach the engine.
- Positioning: no position, velocity or time solution.
- The GPS Integrity Support Message (CNAV Message Type 40, IS-GPS-200N with
  IRN-IS-200N-004 section 30.3.3.10): its parameters serve ARAIM, the
  integrity of a position solution.
- I/O: no sockets, NTRIP, serial ports, files or clocks, so that every run
  replays from its inputs.
- Bit parsing and decryption: these are the gnsswire codecs'.
- Ephemerides from the virtual base station: an RTK receiver decodes its
  own broadcast data.
- Single-precision kernels: the 0.1 mm target needs double precision.
