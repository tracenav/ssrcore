/* State-space corrections of the satellite orbit and clock (SPARTN 8.7, 8.8).
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers). Vectors are ECEF,
 * metres and metres per second, as caller-owned double[3].
 */
#ifndef SSRCORE_CORRECTION_H
#define SSRCORE_CORRECTION_H

#include <stdint.h>

/* [spec: SPARTN v2.0.3 Table 8.0] Speed of light, m/s. */
#define SSR_SPEED_OF_LIGHT 299792458.0

/* Status. [design: DESIGN 3.1, errors are int32_t status codes] */
#define SSR_CORRECTION_OK 0
#define SSR_CORRECTION_DEGENERATE 1 /* |r_dot| = 0 or |r x r_dot| = 0: the frame is undefined; outputs untouched */

/* [spec: SPARTN v2.0.3 section 8.8] Unit vectors of the orbit frame from the
 * broadcast position r and velocity r_dot:
 *   e_along = r_dot / |r_dot|
 *   e_cross = (r x r_dot) / |r x r_dot|
 *   e_radial = e_along x e_cross */
int32_t ssr_correction_orbit_frame(const double r[3], const double r_dot[3],
                                   double e_radial[3], double e_along[3], double e_cross[3]);

/* [spec: SPARTN v2.0.3 section 8.8] Satellite position correction in ECEF,
 * dX = [e_radial e_along e_cross] . dO, for the orbit correction
 * dO = (radial, along, cross) in metres. */
int32_t ssr_correction_orbit(const double r[3], const double r_dot[3],
                             const double d_o[3], double d_x[3]);

/* [spec: SPARTN v2.0.3 section 8.8] Corrected satellite nominal antenna
 * phase centre position X = X_b - dX, with X_b = r. */
int32_t ssr_correction_position(const double r[3], const double r_dot[3],
                                const double d_o[3], double x[3]);

/* [spec: SPARTN v2.0.3 section 8.7] Corrected satellite clock, s:
 * t_s = t_b - dC / c, for the broadcast clock t_b (s) and the clock
 * correction dC (m). The relativistic term of 8.7 is not included: it is
 * the relativity kernel's, and the record keeps it apart (DESIGN 2.4). */
double ssr_correction_clock(double t_b, double d_c);

/* [spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6); IGS SSR v1.00 4.1] Orbit
 * correction at a time: dO(t) = dO + dO_dot * dt, for dO = (radial, along,
 * cross) in metres at the reference time t0, its rates dO_dot in m/s and
 * dt = t - t0 in seconds. */
void ssr_correction_orbit_at(const double d_o[3], const double rate[3], double dt, double out[3]);

/* [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8); IGS SSR v1.00 4.2] Clock
 * correction at a time, m: dC(t) = C0 + C1 dt + C2 dt^2, for C0 (m), C1
 * (m/s), C2 (m/s^2) at the reference time t0 and dt = t - t0 in seconds. */
double ssr_correction_clock_at(double c0, double c1, double c2, double dt);

/* [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8), its derivative in t] Rate of
 * the clock correction, m/s: C1 + 2 C2 dt. */
double ssr_correction_clock_rate_at(double c1, double c2, double dt);

#endif
