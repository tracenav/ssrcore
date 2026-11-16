/* GLONASS broadcast orbit and clock. See ssrcore/glonass.h. */

#include <math.h>
#include <stdint.h>

#include "ssrcore/glonass.h"

/* [spec: GLONASS ICD 5.1 A.3.1.2] PZ-90 constants of the simplified
 * algorithm. */
static const double MU = 398600.44e9;       /* m^3/s^2 */
static const double AE = 6378136.0;         /* m */
static const double J2 = 1082625.7e-9;      /* J0^2 of the ICD */
static const double OMEGA = 7.292115e-5;    /* rad/s */

static const double WEEK_SECONDS = 604800.0;

/* [design: GLONASS ICD 5.1 4.4 gives tb as a time of day in Moscow time;
 * the caller brings tb and t to one scale with weeks (DESIGN T4), so the
 * difference t - tb is exact across days and weeks, taken in seconds before
 * any multiplication (DESIGN 4)] */
static double since(int32_t week, double sow, int32_t ref_week, double ref_sow)
{
    return (double)(week - ref_week) * WEEK_SECONDS + (sow - ref_sow);
}

/* [spec: GLONASS ICD 5.1 A.3.1.2, corrected by RTCM 10403.3 3.5.13: the
 * fourth term of dVy/dt is -2 omega Vx, and the parenthesis of dVz/dt is
 * (3 - 5 z^2 / r^2)]
 *   dVx/dt = -mu/r^3 x - 3/2 J2 mu ae^2/r^5 x (1 - 5 z^2/r^2) + omega^2 x + 2 omega Vy + x''
 *   dVy/dt = -mu/r^3 y - 3/2 J2 mu ae^2/r^5 y (1 - 5 z^2/r^2) + omega^2 y - 2 omega Vx + y''
 *   dVz/dt = -mu/r^3 z - 3/2 J2 mu ae^2/r^5 z (3 - 5 z^2/r^2) + z''                        */
void ssr_glonass_derivative(const double state[6], const double lunisolar[3], double rate[6])
{
    const double x = state[0], y = state[1], z = state[2];
    const double vx = state[3], vy = state[4];
    const double r2 = x * x + y * y + z * z;
    const double r = sqrt(r2);
    const double mu_r3 = MU / (r2 * r);
    const double j2_r5 = 1.5 * J2 * MU * AE * AE / (r2 * r2 * r);
    const double zz = 5.0 * z * z / r2;
    const double omega2 = OMEGA * OMEGA;

    rate[0] = state[3];
    rate[1] = state[4];
    rate[2] = state[5];
    rate[3] = -mu_r3 * x - j2_r5 * x * (1.0 - zz) + omega2 * x + 2.0 * OMEGA * vy + lunisolar[0];
    rate[4] = -mu_r3 * y - j2_r5 * y * (1.0 - zz) + omega2 * y - 2.0 * OMEGA * vx + lunisolar[1];
    rate[5] = -mu_r3 * z - j2_r5 * z * (3.0 - zz) + lunisolar[2];
}

/* [spec: GLONASS ICD 5.1 A.3.1.1, 4th order Runge-Kutta; RTCM 10403.3
 * 3.5.13] k1 = f(s), k2 = f(s + h/2 k1), k3 = f(s + h/2 k2), k4 = f(s + h k3),
 * s + h/6 (k1 + 2 k2 + 2 k3 + k4). The luni-solar acceleration is constant
 * over the interval (A.3.1.2). */
void ssr_glonass_rk4_step(const double state[6], const double lunisolar[3], double h,
                          double out[6])
{
    double k1[6], k2[6], k3[6], k4[6], s[6];
    int32_t j;
    ssr_glonass_derivative(state, lunisolar, k1);
    for (j = 0; j < 6; j++)
        s[j] = state[j] + 0.5 * h * k1[j];
    ssr_glonass_derivative(s, lunisolar, k2);
    for (j = 0; j < 6; j++)
        s[j] = state[j] + 0.5 * h * k2[j];
    ssr_glonass_derivative(s, lunisolar, k3);
    for (j = 0; j < 6; j++)
        s[j] = state[j] + h * k3[j];
    ssr_glonass_derivative(s, lunisolar, k4);
    for (j = 0; j < 6; j++)
        out[j] = state[j] + h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
}

/* [spec: GLONASS ICD 5.1 A.3.1.2: the initial conditions are the broadcast
 * state at tb; RTCM 10403.3 3.5.13: RK4, steps of at most 30 s]
 * [design: the span t - tb is cut into the fewest equal steps of at most
 * 30 s, so every step is as long as allowed and no short remainder step
 * is taken; the result does not depend on the direction of a split] */
int32_t ssr_glonass_state(const double state_tb[6], const double lunisolar[3],
                          int32_t tb_week, double tb, int32_t week, double sow,
                          double out[6])
{
    const double span = since(week, sow, tb_week, tb);
    int32_t j;
    if (!(fabs(span) <= SSR_GLONASS_MAX_SPAN))
        return SSR_GLONASS_TOO_FAR;
    for (j = 0; j < 6; j++)
        out[j] = state_tb[j];
    const int32_t steps = (int32_t)ceil(fabs(span) / SSR_GLONASS_MAX_STEP);
    if (steps == 0)
        return SSR_GLONASS_OK;
    const double h = span / (double)steps;
    for (j = 0; j < steps; j++) {
        double next[6];
        int32_t k;
        ssr_glonass_rk4_step(out, lunisolar, h, next);
        for (k = 0; k < 6; k++)
            out[k] = next[k];
    }
    return SSR_GLONASS_OK;
}

/* [spec: GLONASS ICD 5.1 3.3.3: t_GLONASS = t + tau_n(tb) - gamma_n(tb) (t - tb)
 * for the onboard time t, with tau_c left to the time scale (DESIGN 2.7,
 * timescale); 4.4 for tau_n and gamma_n] The satellite clock offset
 * t - t_GLONASS is -tau_n + gamma_n (t - tb). */
double ssr_glonass_clock(double tau_n, double gamma_n,
                         int32_t tb_week, double tb, int32_t week, double sow)
{
    return -tau_n + gamma_n * since(week, sow, tb_week, tb);
}
