"""Host binding of the kernel ``kernels/src/timescale.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; status codes are returned
as they are, and turned into exceptions by :mod:`ssrcore.toolkit.timescale`.
No arithmetic is done here.
"""

from libc.stdint cimport int32_t, uint32_t


cdef extern from "ssrcore/timescale.h":
    int32_t SSR_TIMESCALE_OK
    int32_t SSR_TIMESCALE_RANGE
    int32_t SSR_SYSTEM_GPS
    int32_t SSR_SYSTEM_GLONASS
    int32_t SSR_SYSTEM_GALILEO
    int32_t SSR_SYSTEM_BEIDOU
    int32_t SSR_SYSTEM_QZSS
    double SSR_WEEK_SECONDS
    double SSR_DAY_SECONDS

    void ssr_timescale_normalize(int32_t week, double sow, int32_t *week_out, double *sow_out) nogil
    void ssr_timescale_add(int32_t week, double sow, double dt, int32_t *week_out, double *sow_out) nogil
    double ssr_timescale_diff(int32_t week_a, double sow_a, int32_t week_b, double sow_b) nogil
    int32_t ssr_timescale_week_resolve(int32_t truncated, int32_t bits, int32_t reference, int32_t *week) nogil

    int32_t ssr_timescale_days_from_date(int32_t year, int32_t month, int32_t day, int32_t *days) nogil
    int32_t ssr_timescale_date_from_days(int32_t days, int32_t *year, int32_t *month, int32_t *day) nogil
    int32_t ssr_timescale_day_of_year(int32_t days, int32_t *doy) nogil
    int32_t ssr_timescale_glonass_day(int32_t n4, int32_t nt, int32_t *days) nogil
    int32_t ssr_timescale_glonass_n4_nt(int32_t days, int32_t *n4, int32_t *nt) nogil

    void ssr_timescale_gst_to_gps(int32_t wn, double tow, int32_t *week, double *sow) nogil
    void ssr_timescale_gps_to_gst(int32_t week, double sow, int32_t *wn, double *tow) nogil
    double ssr_timescale_ggto(int32_t wn, double tow, double a0g, double a1g, int32_t wn0g, double t0g) nogil
    void ssr_timescale_bdt_to_gps(int32_t wn, double sow_bdt, int32_t *week, double *sow) nogil
    void ssr_timescale_gps_to_bdt(int32_t week, double sow, int32_t *wn, double *sow_bdt) nogil
    double ssr_timescale_bdt_offset(double a0, double a1, double sow_bdt) nogil
    void ssr_timescale_glonass_to_gps(int32_t days, double sod, int32_t leap_seconds,
                                      int32_t *week, double *sow) nogil
    void ssr_timescale_gps_to_glonass(int32_t week, double sow, int32_t leap_seconds,
                                      int32_t *days, double *sod) nogil
    void ssr_timescale_utc_to_gps(int32_t days, double sod, int32_t leap_seconds,
                                  int32_t *week, double *sow) nogil

    int32_t ssr_timescale_gps_utc(int32_t week, double sow, double dt_ls, double dt_lsf,
                                  int32_t wn_lsf, int32_t dn, double a0, double a1,
                                  int32_t wn_t, double t_ot,
                                  int32_t *days, double *sod, int32_t *which) nogil
    int32_t ssr_timescale_gst_utc(int32_t wn, double tow, double dt_ls, double dt_lsf,
                                  int32_t wn_lsf, int32_t dn, double a0, double a1,
                                  int32_t wn_0t, double t_0t,
                                  int32_t *days, double *sod, int32_t *which) nogil
    int32_t ssr_timescale_bdt_utc(int32_t wn, double sow_bdt, double dt_ls, double dt_lsf,
                                  int32_t wn_lsf, int32_t dn, double a0, double a1,
                                  int32_t *days, double *sod, int32_t *which) nogil

    int32_t ssr_timescale_spartn_tag32_to_gps(uint32_t tag, int32_t scale, int32_t leap_seconds,
                                              int32_t *week, double *sow) nogil
    int32_t ssr_timescale_spartn_tag32_from_gps(int32_t week, double sow, int32_t scale,
                                                int32_t leap_seconds, uint32_t *tag) nogil
    int32_t ssr_timescale_spartn_resolve_tag16(uint32_t tag16, uint32_t reference, uint32_t *tag) nogil


OK = SSR_TIMESCALE_OK
RANGE = SSR_TIMESCALE_RANGE
SYSTEM_GPS = SSR_SYSTEM_GPS
SYSTEM_GLONASS = SSR_SYSTEM_GLONASS
SYSTEM_GALILEO = SSR_SYSTEM_GALILEO
SYSTEM_BEIDOU = SSR_SYSTEM_BEIDOU
SYSTEM_QZSS = SSR_SYSTEM_QZSS
WEEK_SECONDS = SSR_WEEK_SECONDS
DAY_SECONDS = SSR_DAY_SECONDS


def normalize(int32_t week, double sow):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] (week, sow) with sow in [0, 604800)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_normalize(week, sow, &w, &s)
    return w, s


def add(int32_t week, double sow, double dt):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] (week, sow) + dt seconds."""
    cdef int32_t w
    cdef double s
    ssr_timescale_add(week, sow, dt, &w, &s)
    return w, s


def diff(int32_t week_a, double sow_a, int32_t week_b, double sow_b):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] a - b in seconds."""
    return ssr_timescale_diff(week_a, sow_a, week_b, sow_b)


def week_resolve(int32_t truncated, int32_t bits, int32_t reference):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1, 20.3.3.5.2.4; Galileo OS SIS
    ICD 2.2 5.1.2; BDS ICD 5.2.4.4] (status, full week) with these low bits nearest the
    reference."""
    cdef int32_t w = 0
    status = ssr_timescale_week_resolve(truncated, bits, reference, &w)
    return status, w


def days_from_date(int32_t year, int32_t month, int32_t day):
    """[spec: TM-23 A.1.4 (A.4), (A.5)] (status, days from 1980-01-06)."""
    cdef int32_t d = 0
    status = ssr_timescale_days_from_date(year, month, day, &d)
    return status, d


def date_from_days(int32_t days):
    """[spec: TM-23 A.1.4 (A.6), (A.7)] (status, (year, month, day)) of a day count."""
    cdef int32_t y = 0, m = 0, d = 0
    status = ssr_timescale_date_from_days(days, &y, &m, &d)
    return status, (y, m, d)


def day_of_year(int32_t days):
    """[design: ordinal day, 1 January = 1] (status, day of year) of a day count."""
    cdef int32_t doy = 0
    status = ssr_timescale_day_of_year(days, &doy)
    return status, doy


def glonass_day(int32_t n4, int32_t nt):
    """[spec: GLONASS ICD A.3.1.3] (status, day count) of a GLONASS (N4, NT)."""
    cdef int32_t d = 0
    status = ssr_timescale_glonass_day(n4, nt, &d)
    return status, d


def glonass_n4_nt(int32_t days):
    """[spec: GLONASS ICD A.3.1.3] (status, (N4, NT)) of a day count."""
    cdef int32_t n4 = 0, nt = 0
    status = ssr_timescale_glonass_n4_nt(days, &n4, &nt)
    return status, (n4, nt)


def gst_to_gps(int32_t wn, double tow):
    """[spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8] GPS (week, sow) of GST (WN, TOW)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_gst_to_gps(wn, tow, &w, &s)
    return w, s


def gps_to_gst(int32_t week, double sow):
    """[spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8] GST (WN, TOW) of GPS (week, sow)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_gps_to_gst(week, sow, &w, &s)
    return w, s


def ggto(int32_t wn, double tow, double a0g, double a1g, int32_t wn0g, double t0g):
    """[spec: Galileo OS SIS ICD 2.2 5.1.8 Eq. 24] t_Galileo - t_GPS in seconds."""
    return ssr_timescale_ggto(wn, tow, a0g, a1g, wn0g, t0g)


def bdt_to_gps(int32_t wn, double sow_bdt):
    """[spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] GPS (week, sow) of BDT (WN, SOW)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_bdt_to_gps(wn, sow_bdt, &w, &s)
    return w, s


def gps_to_bdt(int32_t week, double sow):
    """[spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] BDT (WN, SOW) of GPS (week, sow)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_gps_to_bdt(week, sow, &w, &s)
    return w, s


def bdt_offset(double a0, double a1, double sow_bdt):
    """[spec: BDS ICD 5.2.4.19 to 5.2.4.21] A0 + A1 tE in seconds."""
    return ssr_timescale_bdt_offset(a0, a1, sow_bdt)


def glonass_to_gps(int32_t days, double sod, int32_t leap_seconds):
    """[spec: GLONASS ICD 3.3.3] GPS (week, sow) of a GLONASS (day count, second of day)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_glonass_to_gps(days, sod, leap_seconds, &w, &s)
    return w, s


def gps_to_glonass(int32_t week, double sow, int32_t leap_seconds):
    """[spec: GLONASS ICD 3.3.3] GLONASS (day count, second of day) of GPS (week, sow)."""
    cdef int32_t d
    cdef double s
    ssr_timescale_gps_to_glonass(week, sow, leap_seconds, &d, &s)
    return d, s


def utc_to_gps(int32_t days, double sod, int32_t leap_seconds):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a] GPS (week, sow) of a UTC
    (day count, second of day)."""
    cdef int32_t w
    cdef double s
    ssr_timescale_utc_to_gps(days, sod, leap_seconds, &w, &s)
    return w, s


def gps_utc(int32_t week, double sow, double dt_ls, double dt_lsf, int32_t wn_lsf,
            int32_t dn, double a0, double a1, int32_t wn_t, double t_ot):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] (status, (day count, second of
    day, case 1-3))."""
    cdef int32_t d = 0, which = 0
    cdef double s = 0.0
    status = ssr_timescale_gps_utc(week, sow, dt_ls, dt_lsf, wn_lsf, dn, a0, a1, wn_t, t_ot,
                                   &d, &s, &which)
    return status, (d, s, which)


def gst_utc(int32_t wn, double tow, double dt_ls, double dt_lsf, int32_t wn_lsf,
            int32_t dn, double a0, double a1, int32_t wn_0t, double t_0t):
    """[spec: Galileo OS SIS ICD 2.2 5.1.7] (status, (day count, second of day, case 1-3))."""
    cdef int32_t d = 0, which = 0
    cdef double s = 0.0
    status = ssr_timescale_gst_utc(wn, tow, dt_ls, dt_lsf, wn_lsf, dn, a0, a1, wn_0t, t_0t,
                                   &d, &s, &which)
    return status, (d, s, which)


def bdt_utc(int32_t wn, double sow_bdt, double dt_ls, double dt_lsf, int32_t wn_lsf,
            int32_t dn, double a0, double a1):
    """[spec: BDS ICD 5.2.4.18] (status, (day count, second of day, case 1-3))."""
    cdef int32_t d = 0, which = 0
    cdef double s = 0.0
    status = ssr_timescale_bdt_utc(wn, sow_bdt, dt_ls, dt_lsf, wn_lsf, dn, a0, a1, &d, &s, &which)
    return status, (d, s, which)


def spartn_tag32_to_gps(uint32_t tag, int32_t scale, int32_t leap_seconds):
    """[spec: SPARTN Table 7.1 TF009] (status, GPS (week, sow)) of a 32-bit time tag."""
    cdef int32_t w = 0
    cdef double s = 0.0
    status = ssr_timescale_spartn_tag32_to_gps(tag, scale, leap_seconds, &w, &s)
    return status, (w, s)


def spartn_tag32_from_gps(int32_t week, double sow, int32_t scale, int32_t leap_seconds):
    """[spec: SPARTN Table 7.1 TF009] (status, 32-bit time tag) of a GPS (week, sow)."""
    cdef uint32_t t = 0
    status = ssr_timescale_spartn_tag32_from_gps(week, sow, scale, leap_seconds, &t)
    return status, t


def spartn_resolve_tag16(uint32_t tag16, uint32_t reference):
    """[spec: SPARTN Table 7.1 TF009] (status, 32-bit tag) of a 16-bit tag nearest the reference."""
    cdef uint32_t t = 0
    status = ssr_timescale_spartn_resolve_tag16(tag16, reference, &t)
    return status, t
