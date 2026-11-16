/* Solid Earth tide displacement of a site.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * "TN36" is IERS Conventions (2010), IERS Technical Note 36; its chapter 7
 * section 7.1.1 is the model, chapter 1 Table 1.1 the constants, chapter 5
 * eq. (5.43) and chapter 6 section 6.2 the tide arguments of step 2.
 *
 * Positions are ECEF metres; the Sun and the Moon are geocentric in the same
 * frame and at the same time as the site (taken as arguments, DESIGN 3.1;
 * the astro kernel computes them). Displacements are ECEF metres, the
 * conventional tide-free model of 7.1.1, which includes the permanent tide
 * [spec: TN36 7.1.1.2]. The site latitude and longitude in the formulae are
 * geocentric [spec: TN36 7.1.1, "geocentric latitude of station"].
 */
#ifndef SSRCORE_TIDES_H
#define SSRCORE_TIDES_H

#include <stdint.h>

/* [spec: TN36 ch. 1 Table 1.1] GM of the Sun and of the Earth (m^3/s^2), the
 * Moon-Earth mass ratio, and the equatorial radius of the Earth (m) of the
 * IERS Conventions, which the formulae of TN36 7.1.1 are written with. The
 * value of GM of the Earth equals SSR_WGS84_GM; it is named for the IERS
 * system because TN36 states it (DESIGN 4). */
#define SSR_IERS_GM_SUN 1.32712442099e20
#define SSR_IERS_GM_EARTH 3.986004418e14
#define SSR_IERS_MOON_EARTH_MASS_RATIO 0.0123000371
#define SSR_IERS_EARTH_RADIUS 6378136.6
/* [spec: TN36 ch. 1 Table 1.1] GM_sun / GM_earth, the mass_ratio of the Sun
 * in the functions below; the Moon's is SSR_IERS_MOON_EARTH_MASS_RATIO. */
#define SSR_IERS_SUN_EARTH_MASS_RATIO (SSR_IERS_GM_SUN / SSR_IERS_GM_EARTH)

/* [spec: TN36 7.1.1 eq. (5) with eq. (2)] Degree 2 in-phase displacement due
 * to one body: Love and Shida numbers h(phi) = h0 + h2 (3 sin^2 phi - 1)/2,
 * l(phi) = l0 + l2 (3 sin^2 phi - 1)/2 (TN36: h0 = 0.6078, h2 = -0.0006,
 * l0 = 0.0847, l2 = 0.0002). */
void ssr_tides_degree2(const double site[3], const double body[3], double mass_ratio,
                      double h0, double h2, double l0, double l2, double out[3]);

/* [spec: TN36 7.1.1 eq. (6)] Degree 3 displacement due to one body, with h3
 * = 0.292 and l3 = 0.015. */
void ssr_tides_degree3(const double site[3], const double body[3], double mass_ratio,
                      double out[3]);

/* [spec: TN36 7.1.1 eqs. (7), (8), (9)] Transverse displacement from the
 * l(1) term, diurnal band with l1_diurnal (TN36: 0.0012) and semidiurnal
 * band with l1_semidiurnal (TN36: 0.0024). */
void ssr_tides_l1(const double site[3], const double body[3], double mass_ratio,
                 double l1_diurnal, double l1_semidiurnal, double out[3]);

/* [spec: TN36 7.1.1 eqs. (10), (11)] Out-of-phase displacement from the
 * imaginary parts of the Love and Shida numbers; diurnal hI = -0.0025,
 * lI = -0.0007 and semidiurnal hI = -0.0022, lI = -0.0007 in TN36. */
void ssr_tides_out_of_phase(const double site[3], const double body[3], double mass_ratio,
                           double hi_diurnal, double li_diurnal,
                           double hi_semidiurnal, double li_semidiurnal, double out[3]);

/* [spec: TN36 7.1.1 eqs. (12), (13), Tables 7.3a and 7.3b; ch. 6 section 6.2
 * for theta_f; ch. 5 eq. (5.43)] Step 2: corrections for the frequency
 * dependence of the Love and Shida numbers. t is TT in Julian centuries
 * from J2000.0, gmst Greenwich mean sidereal time (rad). */
void ssr_tides_step2(const double site[3], double t, double gmst, double out[3]);

/* [spec: TN36 7.1.1, steps 1 and 2 with the nominal values of its summary
 * table] The displacement of the site: eq. (5) for the Moon and the Sun,
 * eq. (6) for the Moon only, eqs. (8) to (11) for both, and step 2. */
void ssr_tides_displacement(const double site[3], const double sun[3], const double moon[3],
                             double t, double gmst, double out[3]);

/* [spec: TN36 7.1.1.2 eqs. (14a), (14b)] The permanent deformation, to add
 * to a conventional tide-free position to obtain a mean-tide position. */
void ssr_tides_permanent(const double site[3], double out[3]);

#endif
