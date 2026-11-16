Version 0.1.0 (2026-11-16)
==========================
* First release: an engine for GNSS state-space corrections
* Systems: GPS, GLONASS, Galileo, BeiDou and QZSS
* SPARTN v2.0.3 OCB, HPAC and GAD, over IP or inside UBX-RXM-PMP; decryption is the caller's
  * Galileo E5 and E6 biases (C8Q, L8Q, C6C, L6C) from the bias mask positions v2.0.3 defines
* IGS SSR (RTCM 4076): orbit, clock, high-rate clock, code and phase biases, URA, VTEC
  * Taken only from a stream declared antenna phase centre
  * BeiDou GEO orbit corrections and SBAS are left out
* Broadcast ephemerides from RTCM 1019, 1020, 1042, 1044, 1045 and 1046, or from navigation words in UBX-RXM-SFRBX
  * GPS leap seconds only from a subframe 4 page of Data ID 2 (IS-GPS-200N with IRN-IS-200N-003 and -004)
* Out: for any satellite, time and position, the line-of-sight record
  * The satellite's state and every correction term apart, with validity, quality and source
  * A term it cannot give is absent with its reason, never zero
* Out: a virtual base station's RTCM 3 messages: MSM7 (or MSM4), 1005, 1033 and 1230
* Atmosphere from an ordered list of providers: SPARTN, the UNB3m model, IGS SSR VTEC or the caller's own
* No I/O and no clock: every run replays from its recorded inputs
* Every function implements a cited public specification or is a marked design choice
* C99 kernels in double precision that compile unchanged for microcontrollers
* Codecs: the gnsswire distributions bits, spartn, rtcm3, ubx, nmea and nav (gnsswire 0.2.0)
* Python 3.11 or later; source distribution only, so installing needs a C compiler
* MIT licence
