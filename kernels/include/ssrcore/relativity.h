/* Relativistic terms of a satellite signal: the eccentricity clock term, the
 * Shapiro path delay and the Sagnac (Earth rotation) term.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers). Positions are ECEF
 * metres, velocities ECEF metres per second, times seconds.
 */
#ifndef SSRCORE_RELATIVITY_H
#define SSRCORE_RELATIVITY_H

#include <stdint.h>

/* c is SSR_SPEED_OF_LIGHT (SPARTN Table 8.0); GM and omega_E are
 * SSR_WGS84_GM and SSR_WGS84_OMEGA_E, the values Ashby 2003 evaluates
 * sections 2 to 10 with (DESIGN 4). */
#include "ssrcore/correction.h"
#include "ssrcore/geodesy.h"

/* Status codes. */
#define SSR_RELATIVITY_OK 0
/* The straight path passes through the Earth's centre (r1 + r2 - l <= 0):
 * Ashby 2003 Eq. (96) is undefined there. */
#define SSR_RELATIVITY_THROUGH_CENTRE 1
/* Satellite and receiver coincide: there is no path. */
#define SSR_RELATIVITY_ZERO_RANGE 2

/* [spec: SPARTN v2.0.3 section 8.7; IS-GPS-200N with IRN-IS-200N-003 and -004
 * 20.3.3.3.3.1, the Control Segment's equivalent expression] Relativistic clock
 * correction dt_r = -2 r.v / c^2, seconds, of a satellite at position r with
 * velocity v (ECEF or ECI alike, IS-GPS-200N with IRN-IS-200N-003 and -004
 * 20.3.3.3.3.1). */
double ssr_relativity_clock(const double position[3], const double velocity[3]);

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS
 * ICD 2.2 section 5.1.4 Eqs. 15 and 16; BDS-SIS-ICD-B1I-3.0 section 5.2.4.9]
 * The same correction in the broadcast-element form dt_r = F e sqrt(A)
 * sin(E_k), seconds, with F = -2 sqrt(mu) / c^2 for the system's mu (m^3/s^2),
 * eccentricity e, sqrt_a = sqrt(A) (m^1/2) and eccentric anomaly E_k (rad). */
double ssr_relativity_clock_kepler(double mu, double e, double sqrt_a, double eccentric_anomaly);

/* [spec: Ashby 2003 section 10, Eq. (96), second term] Shapiro path delay,
 * metres: (2 GM / c^2) ln((r1 + r2 + l) / (r1 + r2 - l)) with r1 = |satellite|,
 * r2 = |receiver|, l = |satellite - receiver|. Returns a status; the delay is
 * written to *delay only when the status is SSR_RELATIVITY_OK. */
int32_t ssr_relativity_shapiro(const double satellite[3], const double receiver[3], double *delay);

/* [spec: Ashby 2003 section 2: the light path is straight in the inertial
 * frame, and the receiver turns with the Earth during the travel time tau]
 * Sagnac term, metres: E = c tau - |r_s - r_r|, with r_s the satellite at
 * transmission in the ECEF frame of its own epoch, r_r the receiver at
 * reception in the ECEF frame of its own, and c tau = |R_z(omega_E tau) r_s
 * - r_r| the distance in the frame of reception, R_z(a) the rotation of the
 * frame by a about z (x' = x cos a + y sin a, y' = -x sin a + y cos a).
 * Eq. (7), omega_E (x_s y_r - y_s x_r) / c, is the first-order term of E
 * in omega_E tau; the exact form is evaluated [design: the remainder
 * reaches 0.10 mm for GPS, 0.13 mm for Galileo and 0.28 mm for GEO and
 * IGSO orbits above 5 degrees elevation (tests/kernels/test_relativity.py),
 * at or above 0.1 mm]. Returns
 * SSR_RELATIVITY_OK and writes *delta, or SSR_RELATIVITY_ZERO_RANGE. */
int32_t ssr_relativity_sagnac(const double satellite[3], const double receiver[3], double *delta);

#endif
