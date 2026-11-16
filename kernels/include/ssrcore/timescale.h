/* Time scales: GPS, Galileo (GST), BeiDou (BDT), GLONASS and UTC; weeks,
 * calendar days, SPARTN time tags.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * Representation [design: DESIGN 4, an integer week and a double second of
 * week keep 2^-34 s (58 ps) or better for any epoch, enough for arithmetic
 * that multiplies a time difference by a rate, not by c; and one
 * representation for every scale makes each conversion a plain offset]: a time of scale X is
 * (week, sow), the weeks and seconds elapsed since 1980-01-06 00:00:00 read
 * on the clock of scale X, counting 86 400 s to the day. For GPS time that is
 * the GPS week and time of week (IS-GPS-200N with IRN-IS-200N-003 and -004
 * 3.3.4, 6.2.4). The native week numbers of GST and BDT are converted by
 * ssr_timescale_gst_to_gps and ssr_timescale_bdt_to_gps. A calendar day (UTC,
 * or GLONASS time) is a day count since 1980-01-06 and a second of day; the
 * second of day reaches 86 400 only inside a leap second.
 *
 * Every function returning int32_t returns SSR_TIMESCALE_OK or an error code;
 * values are written through the pointers only on SSR_TIMESCALE_OK.
 */
#ifndef SSRCORE_TIMESCALE_H
#define SSRCORE_TIMESCALE_H

#include <stdint.h>

/* Status codes (DESIGN 3.1). */
#define SSR_TIMESCALE_OK 0
/* [design: the documents give each argument a domain (week widths, calendar
 * range, DN, tag range); a value outside it has no result] */
#define SSR_TIMESCALE_RANGE 1

/* [spec: SPARTN Table 6.1, the OCB and HPAC subtypes 0 to 4] A satellite
 * system, numbered as SPARTN numbers the subtypes of its messages. It names
 * the time scale of a SPARTN message (TF009: the constellation's system time,
 * GPS for messages of no or several constellations) and the ICD whose
 * constants an orbit kernel uses. [design: one enumeration for every kernel,
 * named for the quantity (DESIGN 4); the kepler kernel includes this header
 * for it] */
#define SSR_SYSTEM_GPS 0
#define SSR_SYSTEM_GLONASS 1
#define SSR_SYSTEM_GALILEO 2
#define SSR_SYSTEM_BEIDOU 3
#define SSR_SYSTEM_QZSS 4

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] One week, 604 800 s. */
#define SSR_WEEK_SECONDS 604800.0
/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] One day, 86 400 s. */
#define SSR_DAY_SECONDS 86400.0

/* ---- weeks and seconds of week -------------------------------------- */

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, end of week
 * crossover] (week, sow) with sow brought into [0, 604800) by whole weeks. */
void ssr_timescale_normalize(int32_t week, double sow, int32_t *week_out, double *sow_out);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] (week, sow) +
 * dt seconds, normalized. */
void ssr_timescale_add(int32_t week, double sow, double dt, int32_t *week_out, double *sow_out);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, 20.3.3.5.2.4:
 * 604800 (WN - WNt) + t - tot] a - b in seconds. */
double ssr_timescale_diff(int32_t week_a, double sow_a, int32_t week_b, double sow_b);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1, 20.3.3.5.2.4;
 * Galileo OS SIS ICD 2.2 5.1.2, 5.1.7, 5.1.8; BDS ICD 5.2.4.4, 5.2.4.18] The
 * full week whose `bits` low bits are `truncated`, nearest `reference` (a full
 * week of the same numbering). */
int32_t ssr_timescale_week_resolve(int32_t truncated, int32_t bits, int32_t reference, int32_t *week);

/* ---- calendar ---------------------------------------------------------- */

/* [spec: TM-23 A.1.4 (A.4), (A.5); 3.1.1 for 1980-01-06 = JD 2 444 244.5]
 * Days from 1980-01-06 to the civil date. Valid 1900-03-01 to 2100-02-28. */
int32_t ssr_timescale_days_from_date(int32_t year, int32_t month, int32_t day, int32_t *days);

/* [spec: TM-23 A.1.4 (A.6), (A.7)] Civil date of a day count from
 * 1980-01-06. Valid 1900-03-01 to 2100-02-28. */
int32_t ssr_timescale_date_from_days(int32_t days, int32_t *year, int32_t *month, int32_t *day);

/* [design: ordinal day, 1 January = 1, as the day count from 1 January of
 * the same year by ssr_timescale_days_from_date] Day of year of a day count. */
int32_t ssr_timescale_day_of_year(int32_t days, int32_t *doy);

/* [spec: GLONASS ICD A.3.1.3; 4.5 for N4, from 1996 = 1] Day count from
 * 1980-01-06 of the GLONASS date given by the four-year interval N4 and the
 * day NT within it. */
int32_t ssr_timescale_glonass_day(int32_t n4, int32_t nt, int32_t *days);

/* [spec: GLONASS ICD A.3.1.3, inverted] N4 and NT of a day count. */
int32_t ssr_timescale_glonass_n4_nt(int32_t days, int32_t *n4, int32_t *nt);

/* ---- scales ------------------------------------------------------------ */

/* [spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8 Eq. 24] GPS time of a GST
 * (WN, TOW), nominal scales: week + 1024, same seconds. */
void ssr_timescale_gst_to_gps(int32_t wn, double tow, int32_t *week, double *sow);
/* [spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8 Eq. 24] GST (WN, TOW) of a GPS
 * time, the inverse of ssr_timescale_gst_to_gps. */
void ssr_timescale_gps_to_gst(int32_t week, double sow, int32_t *wn, double *tow);

/* [spec: Galileo OS SIS ICD 2.2 5.1.8 Eq. 24] GGTO, t_Galileo - t_GPS =
 * A0G + A1G (TOW - t0G + 604800 (WN - WN0G)), weeks untruncated. */
double ssr_timescale_ggto(int32_t wn, double tow, double a0g, double a1g, int32_t wn0g, double t0g);

/* [spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] GPS time of a BDT (WN, SOW),
 * nominal scales: week + 1356, seconds + 14. */
void ssr_timescale_bdt_to_gps(int32_t wn, double sow_bdt, int32_t *week, double *sow);
/* [spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] BDT (WN, SOW) of a GPS time, the
 * inverse of ssr_timescale_bdt_to_gps. */
void ssr_timescale_gps_to_bdt(int32_t week, double sow, int32_t *wn, double *sow_bdt);

/* [spec: BDS ICD 5.2.4.19 to 5.2.4.21] Delta t = A0 + A1 tE, the broadcast
 * offset of BDT to GPS time, Galileo time or GLONASS time; tE the BDT SOW. */
double ssr_timescale_bdt_offset(double a0, double a1, double sow_bdt);

/* [spec: GLONASS ICD 3.3.3; IS-GPS-200N with IRN-IS-200N-003 and -004
 * 20.3.3.5.2.4] GPS time of a GLONASS time (day count from 1980-01-06, second
 * of day): UTC(SU) = GLONASS - 3 h, GPS = UTC + leap_seconds (GPS - UTC, delta
 * t_LS). */
void ssr_timescale_glonass_to_gps(int32_t days, double sod, int32_t leap_seconds,
                        int32_t *week, double *sow);
/* [spec: GLONASS ICD 3.3.3; IS-GPS-200N with IRN-IS-200N-003 and -004
 * 20.3.3.5.2.4] GLONASS time (day count, second of day) of a GPS time, the
 * inverse of ssr_timescale_glonass_to_gps. */
void ssr_timescale_gps_to_glonass(int32_t week, double sow, int32_t leap_seconds,
                        int32_t *days, double *sod);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a, solved for
 * tE with A0 = A1 = 0] GPS time of a UTC (day count, second of day) given GPS -
 * UTC. A second of day of 86 400 or more is the leap second 23:59:60, given
 * with the count before it. */
void ssr_timescale_utc_to_gps(int32_t days, double sod, int32_t leap_seconds,
                    int32_t *week, double *sow);

/* UTC from a satellite system's time and its broadcast UTC parameters, all
 * weeks untruncated (ssr_timescale_week_resolve), dt_ls and dt_lsf in seconds. The
 * result is a UTC day count from 1980-01-06 and second of day; `which` is the
 * case of the document used, 1 (a), 2 (b, the leap second window) or 3 (c). */

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 cases a, b, c]
 * GPS time; DN 1 to 7, day one the first day of week WN_LSF; window six hours
 * either side of the end of day DN; A1 (tE - tot + 604800 (WN - WNt)). */
int32_t ssr_timescale_gps_utc(int32_t week, double sow, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t wn_t, double t_ot,
                    int32_t *days, double *sod, int32_t *which);

/* [spec: Galileo OS SIS ICD 2.2 5.1.7 Eq. 21 to 23] GST (WN, TOW) and GST
 * weeks WN_LSF, WN0t; otherwise as ssr_timescale_gps_utc. */
int32_t ssr_timescale_gst_utc(int32_t wn, double tow, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t wn_0t, double t_0t,
                    int32_t *days, double *sod, int32_t *which);

/* [spec: BDS ICD 5.2.4.18 cases 1 to 3] BDT (WN, SOW) and BDT week WN_LSF;
 * DN 0 to 6, the day of week; window from DN + 2/3 to DN + 5/4 days;
 * delta t_UTC = dt + A0UTC + A1UTC tE with tE the SOW. */
int32_t ssr_timescale_bdt_utc(int32_t wn, double sow_bdt, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t *days, double *sod, int32_t *which);

/* ---- SPARTN time tags -------------------------------------------------- */

/* [spec: SPARTN Table 7.1 TF009] GPS time of a 32-bit time tag, the seconds
 * since 2010-01-01 00:00:00 of the time scale of system `scale` (an
 * SSR_SYSTEM_ code). leap_seconds (GPS - UTC) is used
 * for GLONASS only. */
int32_t ssr_timescale_spartn_tag32_to_gps(uint32_t tag, int32_t scale, int32_t leap_seconds,
                                int32_t *week, double *sow);

/* [spec: SPARTN Table 7.1 TF009] The 32-bit tag, in the scale, of a GPS time,
 * whole seconds rounded down. */
int32_t ssr_timescale_spartn_tag32_from_gps(int32_t week, double sow, int32_t scale,
                                  int32_t leap_seconds, uint32_t *tag);

/* [spec: SPARTN Table 7.1 TF009, "the half-day must be resolved with the
 * assistance of an unambiguous time tag"] The 32-bit tag of a 16-bit
 * half-day tag (seconds of the half day, 0 to 43 199), the half day being the
 * one that puts it nearest the 32-bit reference tag. */
int32_t ssr_timescale_spartn_resolve_tag16(uint32_t tag16, uint32_t reference, uint32_t *tag);

#endif
