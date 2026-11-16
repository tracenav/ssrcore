/* GLONASS broadcast orbit and clock: integration of the broadcast state.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * A state is position (m) and velocity (m/s) in PZ-90, in the order
 * x, y, z, vx, vy, vz. The luni-solar acceleration (m/s^2) is the one
 * broadcast with the state, in PZ-90.
 *
 * Times are an integer week and seconds of week; tb and the time of
 * evaluation are in the same time scale and week count. Bringing tb, which
 * GLONASS broadcasts in Moscow time (GLONASS ICD 3.3.3, 4.4), to that scale
 * is the time scale's work (DESIGN T4), not this kernel's. The kernels do
 * not judge whether an ephemeris is still valid at that time; validity is
 * the model's (DESIGN 2.5).
 */
#ifndef SSRCORE_GLONASS_H
#define SSRCORE_GLONASS_H

#include <stdint.h>

/* Status of an evaluation. */
#define SSR_GLONASS_OK       0
#define SSR_GLONASS_TOO_FAR  1  /* |t - tb| above SSR_GLONASS_MAX_SPAN */

/* [design: the integration runs one step per 30 s of t - tb, so a bound
 * on the span bounds the work of one call, which a microcontroller caller
 * needs. One day is 2880 steps and far beyond the 15 min (+ P1/2) the
 * model allows (DESIGN 2.5); a longer span is an error of the caller.] */
#define SSR_GLONASS_MAX_SPAN 86400.0

/* [spec: RTCM 10403.3 3.5.13: RK4 "with the step of numerical integration
 * that does not exceed 30 seconds"] Largest integration step, s. */
#define SSR_GLONASS_MAX_STEP 30.0

/* [spec: GLONASS ICD 5.1 A.3.1.2, with the two corrections of RTCM
 * 10403.3 3.5.13] Time derivative of a state: (velocity, acceleration). */
void ssr_glonass_derivative(const double state[6], const double lunisolar[3], double rate[6]);

/* [spec: GLONASS ICD 5.1 A.3.1.1 (4th order Runge-Kutta), A.3.1.2] One
 * classical Runge-Kutta step of h seconds from state to out. */
void ssr_glonass_rk4_step(const double state[6], const double lunisolar[3], double h,
                          double out[6]);

/* [spec: GLONASS ICD 5.1 A.3.1.2; RTCM 10403.3 3.5.13] State at (week, sow)
 * from the broadcast state at tb, integrated with RK4 in equal steps of at
 * most SSR_GLONASS_MAX_STEP. Returns SSR_GLONASS_OK or SSR_GLONASS_TOO_FAR. */
int32_t ssr_glonass_state(const double state_tb[6], const double lunisolar[3],
                          int32_t tb_week, double tb, int32_t week, double sow,
                          double out[6]);

/* [spec: GLONASS ICD 5.1 3.3.3, 4.4] Satellite clock offset (s) from
 * GLONASS time, -tau_n + gamma_n (t - tb). */
double ssr_glonass_clock(double tau_n, double gamma_n,
                         int32_t tb_week, double tb, int32_t week, double sow);

#endif
