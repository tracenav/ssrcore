# refs: staged specifications

The public documents ssrcore is written from. They are staged here locally
for reading and are not committed: several are copyrighted
(`.gitignore` keeps everything in `refs/` but this file out of git). Each PDF
has a `.txt` beside it, extracted for searching; when the text loses a
figure or a formula, the PDF is authoritative.

Short names are those used in `[spec: …]` citations and in
[docs/DESIGN.md](../docs/DESIGN.md).

| File | Short name | Document | Used for |
|---|---|---|---|
| `SPARTN_v2.0.3.pdf` | SPARTN | SPARTN Interface Control Document, version 2.0.3, November 2025 | field scaling (6.3), continuity controls (6.2), message contents (6.5), time tags (7), orbit, clock, bias, atmosphere and area computations (8); every consistency rule S, C, A of DESIGN 3.3 |
| `RTCM_10403.3.pdf` | RTCM | RTCM 10403.3, Differential GNSS Services, Version 3 | ephemeris message contents (1019, 1020, 1042, 1045, 1046); SSR consistency and GLONASS orbit notes (3.5.13); virtual base station: service levels (Table 3.1-2), scheduling (3.1.3, 3.1.7), clock steering (3.1.4), non-physical stations (3.1.5, 3.5.12), MSM (3.5.16), 1005, 1033, 1230 |
| `IGS_SSR_v1.00.pdf` | IGS SSR | IGS State Space Representation (SSR) Format, version 1.00, 2020-10-05 | orbit corrections with rates and the clock polynomial at a reference time (4, 4.1, 4.2; IDF003, IDF004, IDF013 to IDF022), checked against RTCM 3.5.13 for rule C11; the IGS SSR adapter, every message of 4 and 8 with Table 7 (rules IG1 to IG12, DESIGN 3.3), and the VTEC of 4.5.1 and 8.4.1 for `model.igs_vtec` |
| `IS-GPS-200N.pdf` | IS-GPS-200N | Interface Specification IS-GPS-200, revision N, 1 August 2022, without its change notices; cited as "IS-GPS-200N with IRN-IS-200N-003 and -004", the revision as the two notices below amend it | GPS LNAV orbit and clock algorithms (20.3.3.3.3.1, 20.3.3.4.3), fit intervals, GPS time and UTC, Data ID (20.3.3.5.1.1); also the algorithm QZSS LNAV follows until IS-QZSS-PNT is staged |
| `IRN-IS-200N-003.pdf` | IRN-IS-200N-003 | Change notice IRN-IS-200N-003 to IS-GPS-200N, 27 September 2023 (RFC-00502, 2023 changes to the public documents) | LNAV Data ID: other values than 2 no longer invalid (20.3.3.5.1.1, 40.3.3.5.1.1), the Data ID determined before LNAV data are used (6.4.6.4), which the GPS leap seconds follow; its other changes (PRN 33, 6.4.3; cutovers, 20.3.4.1; alarm and marginal conditions, 6.4.6.2.2, 6.4.6.3; CNAV schedules and GGTO ID, 30.3.3, 30.3.4.1) touch nothing ssrcore uses |
| `IRN-IS-200N-004.pdf` | IRN-IS-200N-004 | Change notice IRN-IS-200N-004 to IS-GPS-200N, 16 June 2026 (RFC-000519, civil Integrity Support Message formats) | checked only: CNAV Message Type 40 ISM packet (30.3.3.10) and the nominal and integrity-assured NED URA of CNAV (30.3.3.1.1, 30.3.3.2.4); ssrcore reads no CNAV |
| `Galileo_OS_SIS_ICD_v2.2.pdf` | Galileo OS SIS ICD 2.2 | Galileo Open Service Signal-in-Space ICD, issue 2.2, November 2025 (doi:10.2878/7420173) | Galileo carrier frequencies (Table 2); I/NAV words, through gnsswire-nav (4.3.1, 4.3.3, 4.3.5); orbit and clock (5.1.1 to 5.1.5), GST, GST-UTC and GPS-Galileo offsets (5.1.2, 5.1.7, 5.1.8); issue of data, health and data validity, CRC (5.1.9.2 to 5.1.9.4) |
| `BDS_B1I_ICD_v3.0.pdf` | BDS ICD | BeiDou Signal-in-Space ICD, B1I, version 3.0 | BeiDou D1/D2 orbit and clock including GEO satellites (5.2.4.9 to 5.2.4.12), BDT offsets (5.2.4.18 to 5.2.4.21) |
| `GLONASS_ICD_5.1_en.pdf` | GLONASS ICD | GLONASS Interface Control Document, edition 5.1 | GLONASS time (3.3.3), ephemeris parameters (4.4), orbit integration (A.3.1.2) |
| `IERS_TN36_ch7.pdf` | TN36 | IERS Conventions (2010), IERS Technical Note 36, chapter 7 | solid Earth tide displacement (7.1.1), permanent tide (7.1.1.2) |
| `Q11_sun_moon_tides/elp82b_README`, `elp82b.ps` (notice), `ELP1` … `ELP36`, `elp82b_1`, `elp82b_2` (the authors' Fortran, not read) | ELP 2000-82B | M. Chapront-Touzé and J. Chapront, ELP 2000-82B lunar theory, from https://ftp.imcce.fr/pub/ephem/moon/elp82b/ | Moon position (`astro`, Q11) |
| `Q11_sun_moon_tides/elp82b_ps_strings.txt` | ELP 2000-82B | text strings extracted locally from `elp82b.ps` | searching the notice |
| `Q11_sun_moon_tides/geopos_README`, `MOON.Y05` | IMCCE geopos | IMCCE geocentric Moon positions, from https://ftp.imcce.fr/pub/ephem/moon/geopos/ | test values for the Moon |
| `Q11_sun_moon_tides/VSOP87D.ear`, `vsop87.doc`, `vsop87.chk` | VSOP87 | P. Bretagnon and G. Francou, VSOP87 version D, from https://ftp.imcce.fr/pub/ephem/planets/vsop87/ | Sun position (`astro`, Q11) and its check values |
| `Q11_sun_moon_tides/icc1.pdf`, `icc5.pdf`, `icc6.pdf` (+ `.txt`) | TN36 ch. 1/5/6 | IERS Conventions (2010) chapters 1, 5 and 6, from https://iers-conventions.obspm.fr/content/chapter1/icc1.pdf, …/chapter5/icc5.pdf, …/chapter6/icc6.pdf | constants (Table 1.1), terrestrial frame, ERA, GMST and nutation (ch. 5), tide arguments (ch. 6 section 6.2) |
| `Q11_sun_moon_tides/DEHANTTIDEINEL.F` | DEHANTTIDEINEL | IERS Conventions software for TN36 7.1.1, from https://iers-conventions.obspm.fr/content/chapter7/software/dehanttideinel/DEHANTTIDEINEL.F | only its header test cases were used as test values; none of its code |
| `Ashby_2003_LRR.pdf` | Ashby 2003 | N. Ashby, Relativity in the Global Positioning System, Living Reviews in Relativity 6 (2003) | Sagnac effect (2), satellite clock relativity (5), Shapiro delay (10), as SPARTN 8.16 directs |
| `ESA_GNSS_TM-23_Vol_I.pdf` | TM-23 | J. Sanz Subirana et al., GNSS Data Processing, Vol. I, ESA TM-23/1 | geodesy (3.1.4, Appendix B), Niell mapping functions (5.4.2.2.1), carrier phase wind-up after Wu et al. 1993 (5.5), relativistic terms (5.1.2, 5.2.1); stands in for Niell 1996 and Wu 1993, which are not staged |
| `UNB3m_Leandro_2006.pdf` | UNB3m | R. Leandro, M. Santos, R. B. Langley, UNB Neutral Atmosphere Models: Development and Performance (2006) | height dependence of the hydrostatic delay (DESIGN Q7); the UNB3 table and constants also from TM-23 5.4.2.1 |
| `ublox_F9_HPG_1.51_InterfaceDescription.pdf` | u-blox F9 HPG 1.51 | u-blox F9 HPG 1.51 Interface Description (UBXDOC-963802114-13124 R01, UBX protocol version 27.50) | the meanings of the inputs from an F9 receiver: GNSS and signal identifiers (1.5.2, 1.5.4), NAV-TIMELS (3.15.26), RXM-PMP (3.17.3), RXM-SFRBX (3.17.9), NMEA GGA, ZDA and their fields (2.5.3, 2.5.5, 2.7.5, 2.7.21) and secondary output behind a TAG block (2.8); UBX and NMEA framing (3.2, 3.4, 2.1); the wire format itself is gnsswire's |
| `ublox_X20_HPG_2.11_InterfaceDescription.pdf` | u-blox X20 HPG 2.11 | u-blox X20 HPG 2.11 Interface Description (UBXDOC-304424225-21617 R01, UBX protocol version 50.11) | the same meanings for an X20 receiver: GNSS and signal identifiers (1.5.2, 1.5.4), NAV-TIMELS (3.14.22), RXM-PMP (3.16.3), RXM-SFRBX (3.16.8), NMEA GGA, ZDA and their fields (2.5.3, 2.5.5, 2.7.5, 2.7.20) and secondary output behind a TAG block (2.8); UBX and NMEA framing (3.2, 3.4, 2.1); the wire format itself is gnsswire's |
| `u-blox-X20-HPG-2.11.json` | u-blox X20 HPG 2.11 | the same protocol in machine-readable form, from https://github.com/u-blox/u-blox-X20-interface-description-json, path `HPG/ZED-X20P/u-blox-X20-HPG-2.11.json`; that repository's README permits use, copying, modification and distribution without fee | cross-check of the NAV-TIMELS, RXM-PMP and RXM-SFRBX fields against the PDF |
| `ublox_D9_PMP_1.04_InterfaceDescription.pdf` | u-blox D9 PMP 1.04 | u-blox D9 PMP 1.04 Interface Description (UBX-21040023) | RXM-PMP version 0, the L-band carrier of SPARTN |

Not staged and needed: IS-QZSS-PNT-003 (QZSS, named by
SPARTN 8.3), and if obtainable Niell (1996) and Wu et al. (1993), which
SPARTN 8.16 names and TM-23 reproduces. USNO season times used by the
`astro` tests were read from https://aa.usno.navy.mil/api/seasons?year=2005
and are not staged. SPARTN 8.3 names older ICD issues than
those staged (IS-GPS-200H, Galileo OS SIS ICD issue 1.3, BeiDou ICD SIS version 2); the
algorithms used here are taken as unchanged between those issues and the
staged ones. Not yet confirmed: IS-GPS-200N 20.3.3.3.3.1 (F and the
relativistic clock term), which neither staged change notice amends,
against IS-GPS-200H, which is not staged.
