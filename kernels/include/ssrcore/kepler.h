/* Keplerian broadcast orbit and clock: GPS (and QZSS), Galileo, BeiDou.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * Units are SI and angles are radians (DESIGN 4). The broadcast angles are
 * in semicircles; ssr_kepler_semicircles converts them with the value of pi
 * the ICDs fix for this purpose.
 *
 * Times are an integer week and seconds of week. The time of evaluation and
 * the reference time (toe, toc) are in the same time scale and the same week
 * count: GPS time for GPS and QZSS, GST for Galileo, BDT for BeiDou. The
 * kernels do not judge whether an ephemeris is still valid at that time;
 * validity is the model's (DESIGN 2.5).
 */
#ifndef SSRCORE_KEPLER_H
#define SSRCORE_KEPLER_H

#include <stdint.h>

/* The system of an orbit, which selects the ICD's constants mu and omega_e,
 * is an SSR_SYSTEM_ code of this header: GPS, Galileo, BeiDou, and QZSS with
 * the IS-GPS-200N with IRN-IS-200N-003 and -004 constants until IS-QZSS-PNT is
 * staged (DESIGN Q14). */
#include "ssrcore/timescale.h"

/* Status of an evaluation (DESIGN 3.1). */
#define SSR_KEPLER_OK             0
/* [design: the equations of Table 20-IV hold for 0 <= e < 1 and A > 0, and
 * for the systems whose constants the kernel has] */
#define SSR_KEPLER_BAD_ELEMENTS   1
/* [design: see KEPLER_TOLERANCE in kepler.c] */
#define SSR_KEPLER_NOT_CONVERGED  2

/* Broadcast orbit elements in SI units, angles in radians. */
typedef struct {
    int32_t system;      /* SSR_SYSTEM_GPS, _GALILEO, _BEIDOU or _QZSS */
    int32_t geo;         /* 1 for a BeiDou GEO satellite (ssr_kepler_beidou_geo), else 0 */
    int32_t toe_week;    /* week of toe, in the system's own week count */
    double toe;          /* s of week */
    double sqrt_a;       /* m^1/2 */
    double e;
    double m0;           /* rad */
    double delta_n;      /* rad/s */
    double omega;        /* argument of perigee, rad */
    double omega0;       /* longitude of ascending node at the weekly epoch, rad */
    double omega_dot;    /* rad/s */
    double i0;           /* rad */
    double idot;         /* rad/s */
    double cuc, cus;     /* rad */
    double crc, crs;     /* m */
    double cic, cis;     /* rad */
} ssr_kepler_orbit;

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2; Galileo OS SIS
 * ICD 2.2 5.1.1 Table 68; BDS-SIS-ICD-B1I-3.0 5.2.4.12 Table 5-11] Radians from
 * semicircles, with pi = 3.1415926535898 as the three ICDs fix it. */
double ssr_kepler_semicircles(double semicircles);

/* [spec: BDS-SIS-ICD-B1I-3.0 4.3 Table 4-1] 1 when the BeiDou SV ID is a GEO
 * satellite (1 to 5, 59 to 63), else 0. */
int32_t ssr_kepler_beidou_geo(int32_t prn);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV]
 * Eccentric anomaly E of mean anomaly m by the table's iteration. Returns
 * SSR_KEPLER_OK or SSR_KEPLER_NOT_CONVERGED. */
int32_t ssr_kepler_eccentric_anomaly(double m, double e, double *ecc_anomaly);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV;
 * Galileo OS SIS ICD 2.2 5.1.1 Table 68; BDS-SIS-ICD-B1I-3.0 5.2.4.12 Table
 * 5-11] Antenna phase centre position (m) and velocity (m/s) in the system's
 * ECEF frame at (week, sow), and the eccentric anomaly (rad) used, for the
 * relativistic clock term. */
int32_t ssr_kepler_state(const ssr_kepler_orbit *orbit, int32_t week, double sow,
                         double position[3], double velocity[3], double *ecc_anomaly);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS
 * ICD 2.2 5.1.4 Eq. 14; BDS-SIS-ICD-B1I-3.0 5.2.4.9] Satellite clock offset (s)
 * of the broadcast polynomial af0 + af1 (t - toc) + af2 (t - toc)^2, without
 * the relativistic term (DESIGN 2.4 keeps dt_rel separate) and without group
 * delays (E4). */
double ssr_kepler_clock(double af0, double af1, double af2,
                        int32_t toc_week, double toc, int32_t week, double sow);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] Clock drift
 * (s/s) of the same polynomial: af1 + 2 af2 (t - toc). */
double ssr_kepler_clock_drift(double af1, double af2,
                              int32_t toc_week, double toc, int32_t week, double sow);

#endif
