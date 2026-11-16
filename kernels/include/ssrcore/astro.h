/* Positions of the Sun and the Moon in the terrestrial frame.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * Sources (DESIGN Q11), short names used in the citations below:
 *   ELP 2000-82B  M. Chapront-Touze, J. Chapront, "Lunar solution ELP
 *                 2000-82B", notice and series files ELP1 to ELP36, IMCCE /
 *                 Bureau des longitudes (1995-1996); the theory of A&A 124,
 *                 50 (1983) with the arguments of A&A 190, 342 (1988).
 *   VSOP87        P. Bretagnon, G. Francou, "VSOP87 files" (BDL-9502),
 *                 version D, file VSOP87D.ear; A&A 202, 309 (1988).
 *   TN36 ch. 5    IERS Conventions (2010), IERS Technical Note 36, chapter 5.
 *
 * Time arguments: t is in Julian centuries of TT from J2000.0 (JD 2451545.0
 * TT); TT is used for TDB [spec: TN36 ch. 5 section 5.7.1, the difference
 * is below 1.7 ms]. Angles are radians, distances metres.
 *
 * Accuracy, from the truncation of the series (see astro.c), against the
 * complete series over 1990-2060: Moon within 4.3 arcsec and 3.3 km, Sun
 * within 2 arcsec and 820 km. Polar motion is not applied (below 1 arcsec).
 */
#ifndef SSRCORE_ASTRO_H
#define SSRCORE_ASTRO_H

#include <stdint.h>

/* [spec: TN36 ch. 1 Table 1.1] The astronomical unit, m. */
#define SSR_IERS_AU 149597870700.0

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4, GPS time and its
 * epoch; TM-23 3.1.1, TAI - GPS = 19 s; TN36 ch. 5 section 5.5.7, TT - TAI =
 * 32.184 s] TT - GPS, s. */
#define SSR_TT_MINUS_GPS 51.184

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4, GPS time and its
 * epoch; TM-23 3.1.1, TAI - GPS = 19 s; TN36 ch. 5 section 5.5.7, TT - TAI =
 * 32.184 s; TN36 ch. 5 section 5.7.1, TT for TDB] TT in Julian centuries from
 * J2000.0 for a GPS week and second of week (weeks counted from 1980-01-06,
 * without roll-over): TT = GPS + 51.184 s. */
double ssr_astro_tt_centuries(int32_t week, double seconds);

/* [spec: TN36 ch. 5 eq. (5.14), (5.15)] Earth Rotation Angle (rad, in
 * [0, 2 pi)) at the GPS week and second of week; ut1_minus_gps is UT1 - GPS
 * in seconds (UT1 - UTC minus the GPS - UTC leap seconds). */
double ssr_astro_earth_rotation_angle(int32_t week, double seconds, double ut1_minus_gps);

/* [spec: TN36 ch. 5 eq. (5.32)] Greenwich mean sidereal time (rad, in
 * [0, 2 pi)), IAU 2006. */
double ssr_astro_gmst(int32_t week, double seconds, double ut1_minus_gps);

/* [spec: TN36 ch. 5 eq. (5.40)] Mean obliquity of the ecliptic of date,
 * epsilon_A (rad). */
double ssr_astro_mean_obliquity(double t);

/* [spec: TN36 ch. 5 Tables 5.3a and 5.3b, the extract printed there; eq.
 * (5.43)] Nutation in longitude and obliquity (rad). */
void ssr_astro_nutation(double t, double *dpsi, double *deps);

/* [spec: ELP 2000-82B sections 2 to 4, 7 and 8] Geocentric Moon: longitude
 * and latitude referred to the mean ecliptic and equinox of date (rad) and
 * distance (m), in out[0], out[1], out[2]. */
void ssr_astro_moon_ecliptic(double t, double out[3]);

/* [spec: VSOP87, version D, "COMPUTATION"] Heliocentric Earth: longitude and
 * latitude referred to the mean ecliptic and equinox of date (rad) and
 * radius (au), in out[0], out[1], out[2]. Longitude is not reduced. */
void ssr_astro_earth_heliocentric(double t, double out[3]);

/* [spec: VSOP87 version D; TN36 ch. 1 Table 1.1, the astronomical unit]
 * Geocentric Sun: longitude and latitude referred to the mean ecliptic and
 * equinox of date (rad) and distance (m): the heliocentric Earth reversed. */
void ssr_astro_sun_ecliptic(double t, double out[3]);

/* [spec: TN36 ch. 5 section 5.4.5, the classical nutation matrix] From
 * longitude, latitude (mean ecliptic and equinox of date, rad) and distance
 * to rectangular coordinates on the true equator and equinox of date: the
 * longitude is moved by dpsi and the result rotated by the true obliquity
 * eps_a + deps. */
void ssr_astro_ecliptic_to_true_equator(const double lbr[3], double dpsi, double eps_a,
                                  double deps, double out[3]);

/* [spec: TN36 ch. 5 section 5.4.3, eq. (5.31) without its complementary
 * terms] From the true equator and equinox of date to the terrestrial frame:
 * a rotation about z by Greenwich apparent sidereal time,
 * GST = GMST + dpsi cos(eps_a). Polar motion is not applied [design: below
 * 1 arcsec, which moves a tide by under 1 micrometre]. */
void ssr_astro_true_equator_to_ecef(const double v[3], double gst, double out[3]);

/* [spec: the functions above, composed] Sun and Moon in ECEF (m) at the GPS
 * week and second of week, with ut1_minus_gps as in
 * ssr_astro_earth_rotation_angle. */
void ssr_astro_sun_moon_ecef(int32_t week, double seconds, double ut1_minus_gps,
                          double sun[3], double moon[3]);

#endif
