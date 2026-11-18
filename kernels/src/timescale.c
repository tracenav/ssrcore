/* Time scales, weeks, calendar days, SPARTN time tags. See ssrcore/timescale.h. */

#include <math.h>
#include <stdint.h>

#include "ssrcore/timescale.h"

/* [spec: Galileo OS SIS ICD 2.2 5.1.2] GST week 0 starts 13 s before
 * 1999-08-22 00:00:00 UTC; GPS - UTC was then 13 s (IERS Bulletin C), so it
 * starts with GPS week 1024, and 5.1.8 Eq. 24 leaves only the broadcast
 * fraction A0G, A1G between the scales. RTCM DF289 repeats the epoch. */
#define GST_WEEK0 1024
/* [spec: BDS ICD 3.3, 5.2.4.4] BDT week 0 starts 2006-01-01 00:00:00 UTC,
 * a Sunday 1356 weeks after 1980-01-06; [spec: RTCM DF427] BDT time of week
 * is 14 s less than GPS time of week at the same epoch. */
#define BDT_WEEK0 1356
#define BDT_GPS_SECONDS 14.0
/* [spec: GLONASS ICD 3.3.3] T_GLONASS = T_UTC(SU) + 03 h 00 min. */
#define GLONASS_UTC_SECONDS 10800.0
/* [spec: SPARTN Table 7.1 TF009] 2010-01-01 00:00:00, 10 953 days after
 * 1980-01-06: week 1564, 5 days. */
#define SPARTN_EPOCH_WEEK 1564
#define SPARTN_EPOCH_SOW 432000
/* [spec: SPARTN Table 7.1 TF009] the 16-bit tag ranges 0 to 43 199 s. */
#define HALF_DAY 43200
/* [spec: TM-23 3.1.1] 6 January 1980 0h = JD 2 444 244.5; with the 1 720 981.5
 * of TM-23 (A.4) the day count is the sum of (A.4) less 723 263. */
#define JD_GPS_EPOCH_OFFSET 723263
/* JD + 0.5 of 1980-01-06 0h, the a of TM-23 (A.6) for day count 0. */
#define JD_GPS_EPOCH_A 2444245

/* [design: integer division rounded towards minus infinity, as int[] of the
 * documents is used here only on positive values but week and day counts
 * before an epoch must still fall in the right week or day] */
static int32_t floor_div(int32_t a, int32_t b)
{
    int32_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) {
        q -= 1;
    }
    return q;
}

/* [design: shared by every function that splits a count of seconds into a
 * whole number of periods and a remainder in [0, period)] */
static void split(double seconds, double period, double *whole, double *rest)
{
    double k = floor(seconds / period);
    double r = seconds - k * period;
    /* the quotient may round across an integer; one step restores the range */
    if (r < 0.0) {
        k -= 1.0;
        r += period;
    } else if (r >= period) {
        k += 1.0;
        r -= period;
    }
    *whole = k;
    *rest = r;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, end of week crossover] */
void ssr_timescale_normalize(int32_t week, double sow, int32_t *week_out, double *sow_out)
{
    double k, r;
    split(sow, SSR_WEEK_SECONDS, &k, &r);
    *week_out = week + (int32_t)k;
    *sow_out = r;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] */
void ssr_timescale_add(int32_t week, double sow, double dt, int32_t *week_out, double *sow_out)
{
    ssr_timescale_normalize(week, sow + dt, week_out, sow_out);
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4, tE - tot +
 * 604800 (WN - WNt)] the week difference is exact in double, and the seconds
 * are subtracted before they are added to it (DESIGN 4). */
double ssr_timescale_diff(int32_t week_a, double sow_a, int32_t week_b, double sow_b)
{
    return (double)(week_a - week_b) * SSR_WEEK_SECONDS + (sow_a - sow_b);
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1 (WN modulo
 * 1024), 20.3.3.5.2.4 (WNt, WN_LSF modulo 256, within 127 of WN); Galileo OS
 * SIS ICD 2.2 5.1.2 (WN modulo 4096), 5.1.7, 5.1.8 (WN0G within 31 of WN); BDS
 * ICD 5.2.4.4 (WN modulo 8192), 5.2.4.18]
 * [design: the full week in [reference - M/2, reference + M/2), M = 2^bits,
 * the one interval in which the documents' bounds on the difference place
 * it] */
int32_t ssr_timescale_week_resolve(int32_t truncated, int32_t bits, int32_t reference, int32_t *week)
{
    int32_t modulus, offset;
    if (bits < 1 || bits > 30) {
        return SSR_TIMESCALE_RANGE;
    }
    modulus = (int32_t)1 << bits;
    if (truncated < 0 || truncated >= modulus) {
        return SSR_TIMESCALE_RANGE;
    }
    offset = (truncated - reference) % modulus;   /* in (-M, M) */
    if (offset < 0) {
        offset += modulus;                          /* in [0, M) */
    }
    if (offset >= modulus / 2) {
        offset -= modulus;                          /* in [-M/2, M/2) */
    }
    *week = reference + offset;
    return SSR_TIMESCALE_OK;
}

/* [spec: TM-23 A.1.4: (A.4), (A.5) hold from 1 March 1900 to 28 February
 * 2100, where the Julian leap-year rule they use agrees with the Gregorian
 * calendar] [design: dates outside it are refused] */
static int32_t date_in_range(int32_t year, int32_t month, int32_t day)
{
    int32_t key = year * 10000 + month * 100 + day;
    return key >= 19000301 && key <= 21000228;
}

/* [spec: TM-23 A.1.4 (A.4), (A.5)] JD = int[365.25 y] + int[30.6001 (m+1)] +
 * D + 1 720 981.5, written in exact integer arithmetic: int[365.25 y] =
 * (1461 y) div 4 and int[30.6001 (m+1)] = (306001 (m+1)) div 10000 for the
 * positive y and m of the valid range. [design: the date must exist; it is
 * checked by converting back] */
int32_t ssr_timescale_days_from_date(int32_t year, int32_t month, int32_t day, int32_t *days)
{
    int32_t y, m, d, yy, mm, dd;
    if (month < 1 || month > 12 || day < 1 || day > 31
        || !date_in_range(year, month, day)) {
        return SSR_TIMESCALE_RANGE;
    }
    y = (month <= 2) ? year - 1 : year;
    m = (month <= 2) ? month + 12 : month;
    d = (1461 * y) / 4 + (306001 * (m + 1)) / 10000 + day - JD_GPS_EPOCH_OFFSET;
    if (ssr_timescale_date_from_days(d, &yy, &mm, &dd) != SSR_TIMESCALE_OK
        || yy != year || mm != month || dd != day) {
        return SSR_TIMESCALE_RANGE;   /* 31 April, 29 February of a common year */
    }
    *days = d;
    return SSR_TIMESCALE_OK;
}

/* [spec: TM-23 A.1.4 (A.6), (A.7)] a = int[JD + 0.5], b = a + 1537,
 * c = int[(b - 122.1)/365.25], d = int[365.25 c], e = int[(b - d)/30.6001],
 * D = b - d - int[30.6001 e], M = e - 1 - 12 int[e/14],
 * Y = c - 4715 - int[(7 + M)/10]; in exact integer arithmetic,
 * (b - 122.1)/365.25 = (100 b - 12210)/36525 and (b - d)/30.6001 =
 * 10000 (b - d)/306001. A day count has frac[JD + 0.5] = 0. */
int32_t ssr_timescale_date_from_days(int32_t days, int32_t *year, int32_t *month, int32_t *day)
{
    int32_t a, b, c, d, e, dd, mm, yy;
    /* 1900-03-01 and 2100-02-28 are days -29 165 and 43 883 */
    if (days < -29165 || days > 43883) {
        return SSR_TIMESCALE_RANGE;
    }
    a = days + JD_GPS_EPOCH_A;
    b = a + 1537;
    c = (100 * b - 12210) / 36525;
    d = (1461 * c) / 4;
    e = (10000 * (b - d)) / 306001;
    dd = b - d - (306001 * e) / 10000;
    mm = e - 1 - 12 * (e / 14);
    yy = c - 4715 - (7 + mm) / 10;
    *year = yy;
    *month = mm;
    *day = dd;
    return SSR_TIMESCALE_OK;
}

/* [design: ordinal day, 1 January = 1; 1 January of 1900 and of 2100 lie
 * outside the range of TM-23 (A.4), so those years are refused] */
int32_t ssr_timescale_day_of_year(int32_t days, int32_t *doy)
{
    int32_t y, m, d, jan1;
    if (ssr_timescale_date_from_days(days, &y, &m, &d) != SSR_TIMESCALE_OK
        || ssr_timescale_days_from_date(y, 1, 1, &jan1) != SSR_TIMESCALE_OK) {
        return SSR_TIMESCALE_RANGE;
    }
    *doy = days - jan1 + 1;
    return SSR_TIMESCALE_OK;
}

/* [spec: GLONASS ICD A.3.1.3] first NT of the years J = 1 to 4 of the
 * four-year interval: 1, 367, 732, 1097. */
static const int32_t GLONASS_YEAR_START[4] = {1, 367, 732, 1097};

/* [spec: GLONASS ICD A.3.1.3] J from NT, Y = 1996 + 4 (N4 - 1) + (J - 1);
 * the day and month, which the ICD takes from a table of the civil calendar,
 * follow from the day of year NT - NT(J) + 1. [spec: GLONASS ICD 4.5] N4 is
 * 5 bits, 1 to 31; NT is 1 to 1461. */
int32_t ssr_timescale_glonass_day(int32_t n4, int32_t nt, int32_t *days)
{
    int32_t j, year, jan1;
    if (n4 < 1 || n4 > 31 || nt < 1 || nt > 1461) {
        return SSR_TIMESCALE_RANGE;
    }
    j = (nt >= 1097) ? 4 : (nt >= 732) ? 3 : (nt >= 367) ? 2 : 1;
    year = 1996 + 4 * (n4 - 1) + (j - 1);
    if (ssr_timescale_days_from_date(year, 1, 1, &jan1) != SSR_TIMESCALE_OK) {
        return SSR_TIMESCALE_RANGE;
    }
    *days = jan1 + (nt - GLONASS_YEAR_START[j - 1]);
    return SSR_TIMESCALE_OK;
}

/* [spec: GLONASS ICD A.3.1.3, the same relations solved for N4 and NT] */
int32_t ssr_timescale_glonass_n4_nt(int32_t days, int32_t *n4, int32_t *nt)
{
    int32_t y, m, d, jan1, k;
    if (ssr_timescale_date_from_days(days, &y, &m, &d) != SSR_TIMESCALE_OK || y < 1996
        || ssr_timescale_days_from_date(y, 1, 1, &jan1) != SSR_TIMESCALE_OK) {
        return SSR_TIMESCALE_RANGE;
    }
    k = y - 1996;
    *n4 = k / 4 + 1;
    *nt = GLONASS_YEAR_START[k % 4] + (days - jan1);
    return SSR_TIMESCALE_OK;
}

/* [spec: Galileo OS SIS ICD 2.2 5.1.2] */
void ssr_timescale_gst_to_gps(int32_t wn, double tow, int32_t *week, double *sow)
{
    *week = wn + GST_WEEK0;
    *sow = tow;
}

/* [spec: Galileo OS SIS ICD 2.2 5.1.2] */
void ssr_timescale_gps_to_gst(int32_t week, double sow, int32_t *wn, double *tow)
{
    *wn = week - GST_WEEK0;
    *tow = sow;
}

/* [spec: Galileo OS SIS ICD 2.2 5.1.8 Eq. 24] */
double ssr_timescale_ggto(int32_t wn, double tow, double a0g, double a1g, int32_t wn0g, double t0g)
{
    return a0g + a1g * ssr_timescale_diff(wn, tow, wn0g, t0g);
}

/* [spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] */
void ssr_timescale_bdt_to_gps(int32_t wn, double sow_bdt, int32_t *week, double *sow)
{
    ssr_timescale_normalize(wn + BDT_WEEK0, sow_bdt + BDT_GPS_SECONDS, week, sow);
}

/* [spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] */
void ssr_timescale_gps_to_bdt(int32_t week, double sow, int32_t *wn, double *sow_bdt)
{
    int32_t w;
    ssr_timescale_normalize(week, sow - BDT_GPS_SECONDS, &w, sow_bdt);
    *wn = w - BDT_WEEK0;
}

/* [spec: BDS ICD 5.2.4.19 to 5.2.4.21] */
double ssr_timescale_bdt_offset(double a0, double a1, double sow_bdt)
{
    return a0 + a1 * sow_bdt;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a with
 * A0 = A1 = 0: tE = tUTC + delta t_LS] */
void ssr_timescale_utc_to_gps(int32_t days, double sod, int32_t leap_seconds,
                    int32_t *week, double *sow)
{
    int32_t w = floor_div(days, 7);
    double s = (double)(days - 7 * w) * SSR_DAY_SECONDS + sod + (double)leap_seconds;
    ssr_timescale_normalize(w, s, week, sow);
}

/* [spec: GLONASS ICD 3.3.3: UTC(SU) = GLONASS - 3 h; then as ssr_timescale_utc_to_gps] */
void ssr_timescale_glonass_to_gps(int32_t days, double sod, int32_t leap_seconds,
                        int32_t *week, double *sow)
{
    ssr_timescale_utc_to_gps(days, sod - GLONASS_UTC_SECONDS, leap_seconds, week, sow);
}

/* [spec: GLONASS ICD 3.3.3; IS-GPS-200N with IRN-IS-200N-003 and -004
 * 20.3.3.5.2.4 a with A0 = A1 = 0] [design: the second of day is in [0, 86400);
 * a GLONASS leap second has no GPS time of its own to come from] */
void ssr_timescale_gps_to_glonass(int32_t week, double sow, int32_t leap_seconds,
                        int32_t *days, double *sod)
{
    double k, r;
    split(sow - (double)leap_seconds + GLONASS_UTC_SECONDS, SSR_DAY_SECONDS, &k, &r);
    *days = 7 * week + (int32_t)k;
    *sod = r;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4; Galileo OS SIS
 * ICD 2.2 5.1.7; BDS ICD 5.2.4.18] Case a (and 1): the leap second is not in
 * the past and t is outside the window, or dt_ls = dt_lsf: tUTC = (tE - delta
 * t_UTC) mod 86400. Case b (2): t inside the window: W = (tE - delta t_UTC -
 * 43200) mod 86400 + 43200, tUTC = W mod (86400 + dt_lsf - dt_ls), delta t_UTC
 * with dt_ls. Case c (3): past and outside the window: case a with dt_lsf.
 * [design: the day count carries with the modulo, as 20.3.3.5.2.4 b asks
 * ("carries or borrows into any year/week/day counts"); tE is taken as its
 * second of week so that no large total of seconds loses resolution]
 * (week, sow) on the representation of the header; ref is the reference
 * time of A1, event the end of day DN, before and after the window. */
static int32_t utc_from_scale(int32_t week, double sow,
                              double dt_ls, double dt_lsf,
                              double a0, double a1, int32_t ref_week, double ref_sow,
                              int32_t event_week, double event_sow,
                              double before, double after,
                              int32_t *days, double *sod)
{
    double poly = a0 + a1 * ssr_timescale_diff(week, sow, ref_week, ref_sow);
    double since = ssr_timescale_diff(week, sow, event_week, event_sow);
    double k, r, u;
    int32_t which;

    if (dt_ls == dt_lsf || since < -before) {
        which = 1;
    } else if (since <= after) {
        which = 2;
    } else {
        which = 3;
    }

    if (which == 2) {
        double w, m = SSR_DAY_SECONDS + dt_lsf - dt_ls;
        u = sow - (dt_ls + poly);
        split(u - 43200.0, SSR_DAY_SECONDS, &k, &r);
        w = r + 43200.0;
        if (w >= m) {
            *days = 7 * week + (int32_t)k + 1;
            *sod = w - m;
        } else {
            *days = 7 * week + (int32_t)k;
            *sod = w;
        }
        return which;
    }

    u = sow - (((which == 1) ? dt_ls : dt_lsf) + poly);
    split(u, SSR_DAY_SECONDS, &k, &r);
    *days = 7 * week + (int32_t)k;
    *sod = r;
    return which;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a, b, c] */
int32_t ssr_timescale_gps_utc(int32_t week, double sow, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t wn_t, double t_ot,
                    int32_t *days, double *sod, int32_t *which)
{
    if (dn < 1 || dn > 7) {
        return SSR_TIMESCALE_RANGE;
    }
    *which = utc_from_scale(week, sow, dt_ls, dt_lsf, a0, a1, wn_t, t_ot,
                            wn_lsf, (double)dn * SSR_DAY_SECONDS, 21600.0, 21600.0,
                            days, sod);
    return SSR_TIMESCALE_OK;
}

/* [spec: Galileo OS SIS ICD 2.2 5.1.7 Eq. 21 to 23; 5.1.2 for the GST week]
 * Case a takes delta t_LS (Eq. 21), as in IS-GPS-200N with IRN-IS-200N-003 and
 * -004 20.3.3.5.2.4 a. Issue 2.1 printed delta t_LSF there (its Eq. 20), which
 * this kernel already read as delta t_LS; issue 2.2 corrects it (Document
 * Change Record, "Correction of Eq. 21 updating the term delta t_LS"). */
int32_t ssr_timescale_gst_utc(int32_t wn, double tow, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t wn_0t, double t_0t,
                    int32_t *days, double *sod, int32_t *which)
{
    if (dn < 1 || dn > 7) {
        return SSR_TIMESCALE_RANGE;
    }
    *which = utc_from_scale(wn + GST_WEEK0, tow, dt_ls, dt_lsf, a0, a1,
                            wn_0t + GST_WEEK0, t_0t,
                            wn_lsf + GST_WEEK0, (double)dn * SSR_DAY_SECONDS,
                            21600.0, 21600.0, days, sod);
    return SSR_TIMESCALE_OK;
}

/* [spec: BDS ICD 5.2.4.18, 3.3] The leap second takes effect at the end of
 * day DN, DN + 1 days into week WN_LSF; the window DN + 2/3 to DN + 5/4
 * starts 8 h before and ends 6 h after it. tE is the SOW, so the reference
 * of A1UTC is the start of the current week. [design: DN is read as 0 to 6,
 * the range Table 5-18 gives ("6") and the only reading in which DN + 2/3
 * and DN + 5/4 lie around the end of day DN] */
int32_t ssr_timescale_bdt_utc(int32_t wn, double sow_bdt, double dt_ls, double dt_lsf,
                    int32_t wn_lsf, int32_t dn, double a0, double a1,
                    int32_t *days, double *sod, int32_t *which)
{
    int32_t w;
    double s;
    if (dn < 0 || dn > 6) {
        return SSR_TIMESCALE_RANGE;
    }
    ssr_timescale_normalize(wn + BDT_WEEK0, sow_bdt, &w, &s);
    *which = utc_from_scale(w, s, dt_ls, dt_lsf, a0, a1, w, 0.0,
                            wn_lsf + BDT_WEEK0, (double)(dn + 1) * SSR_DAY_SECONDS,
                            28800.0, 21600.0, days, sod);
    return SSR_TIMESCALE_OK;
}

/* Seconds to add to a reading of the scale's clock to obtain GPS time.
 * [spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8 (GST: none); RTCM DF428 (QZSST: its week
 * begins at the GPS week's midnight, and unlike DF427 for BDT no offset of
 * whole seconds is given); BDS ICD 3.3, RTCM DF427
 * (BDT: 14 s); GLONASS ICD 3.3.3 and IS-GPS-200N with IRN-IS-200N-003 and
 * -004 20.3.3.5.2.4 (GLONASS: - 3 h + GPS - UTC)] */
static int32_t scale_to_gps(int32_t scale, int32_t leap_seconds, double *offset)
{
    switch (scale) {
    case SSR_SYSTEM_GPS:
    case SSR_SYSTEM_QZSS:
    case SSR_SYSTEM_GALILEO:
        *offset = 0.0;
        return SSR_TIMESCALE_OK;
    case SSR_SYSTEM_BEIDOU:
        *offset = BDT_GPS_SECONDS;
        return SSR_TIMESCALE_OK;
    case SSR_SYSTEM_GLONASS:
        *offset = (double)leap_seconds - GLONASS_UTC_SECONDS;
        return SSR_TIMESCALE_OK;
    default:
        return SSR_TIMESCALE_RANGE;
    }
}

/* [spec: SPARTN Table 7.1 TF009: "seconds since January 1st 2010, 00:00:00"
 * in the system time of the message's constellation] [design: for GLONASS,
 * the epoch read on the GLONASS clock and the seconds counted as the GLONASS
 * calendar counts them, 86 400 to the day, the proposal of DESIGN Q2] */
int32_t ssr_timescale_spartn_tag32_to_gps(uint32_t tag, int32_t scale, int32_t leap_seconds,
                                int32_t *week, double *sow)
{
    double offset;
    int64_t total = (int64_t)tag + SPARTN_EPOCH_SOW;
    if (scale_to_gps(scale, leap_seconds, &offset) != SSR_TIMESCALE_OK) {
        return SSR_TIMESCALE_RANGE;
    }
    ssr_timescale_normalize(SPARTN_EPOCH_WEEK + (int32_t)(total / 604800),
                       (double)(total % 604800) + offset, week, sow);
    return SSR_TIMESCALE_OK;
}

/* [spec: SPARTN Table 7.1 TF009: 0 to 4 294 967 295 s from the epoch]
 * [design: a time between two whole seconds has the tag of the second
 * before it; GLONASS as in ssr_timescale_spartn_tag32_to_gps] */
int32_t ssr_timescale_spartn_tag32_from_gps(int32_t week, double sow, int32_t scale,
                                  int32_t leap_seconds, uint32_t *tag)
{
    double offset, t;
    if (scale_to_gps(scale, leap_seconds, &offset) != SSR_TIMESCALE_OK) {
        return SSR_TIMESCALE_RANGE;
    }
    /* [design: the fraction is dropped before the count from 2010, whose double
     * resolves only about 6e-8 s near 2^29 s, so a time just below a whole
     * second would round up to it and get a tag after the time; the offset is
     * whole seconds, so floor(sow) - offset is exact] */
    t = ssr_timescale_diff(week, floor(sow) - offset, SPARTN_EPOCH_WEEK, (double)SPARTN_EPOCH_SOW);
    if (!(t >= 0.0 && t <= 4294967295.0)) {
        return SSR_TIMESCALE_RANGE;
    }
    *tag = (uint32_t)t;
    return SSR_TIMESCALE_OK;
}

/* [spec: SPARTN Table 7.1 TF009] The 16-bit tag is the second within its half
 * day, so the 32-bit tag it expands to (SPARTN 8.15.1 d) is congruent to it
 * modulo 43 200. [design: of the candidates in the reference's half day and
 * the two beside it, the nearest the reference, so that the result is right
 * whenever the reference is within 6 h of the tag (DESIGN T2); at an exact
 * tie the earlier, because a message is received after its time tag;
 * candidates outside the 32-bit range are not taken] */
int32_t ssr_timescale_spartn_resolve_tag16(uint32_t tag16, uint32_t reference, uint32_t *tag)
{
    int64_t ref = (int64_t)reference;
    int64_t base = ref - ref % HALF_DAY;
    int64_t best = -1, best_distance = 0;
    int32_t i;
    if (tag16 >= HALF_DAY) {
        return SSR_TIMESCALE_RANGE;
    }
    for (i = -1; i <= 1; i++) {
        int64_t c = base + (int64_t)i * HALF_DAY + (int64_t)tag16;
        int64_t distance = (c > ref) ? c - ref : ref - c;
        if (c < 0 || c > (int64_t)UINT32_MAX) {
            continue;
        }
        if (best < 0 || distance < best_distance) {
            best = c;
            best_distance = distance;
        }
    }
    *tag = (uint32_t)best;   /* the reference's own half day is always a candidate */
    return SSR_TIMESCALE_OK;
}
