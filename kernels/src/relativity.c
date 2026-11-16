/* Relativistic terms of a satellite signal. See ssrcore/relativity.h. */

#include <math.h>

#include "ssrcore/relativity.h"

/* [spec: SPARTN v2.0.3 section 8.7] dt_r = -2 (r . r_dot) / c^2. */
double ssr_relativity_clock(const double position[3], const double velocity[3])
{
    double rv = position[0] * velocity[0] + position[1] * velocity[1] + position[2] * velocity[2];
    return -2.0 * rv / (SSR_SPEED_OF_LIGHT * SSR_SPEED_OF_LIGHT);
}

/* [spec: Galileo OS SIS ICD 2.2 5.1.4 Eqs. 15, 16; IS-GPS-200N with
 * IRN-IS-200N-003 and -004 20.3.3.3.3.1] dt_r = F e sqrt(A) sin(E_k), F = -2
 * sqrt(mu) / c^2. */
double ssr_relativity_clock_kepler(double mu, double e, double sqrt_a, double eccentric_anomaly)
{
    double f = -2.0 * sqrt(mu) / (SSR_SPEED_OF_LIGHT * SSR_SPEED_OF_LIGHT);
    return f * e * sqrt_a * sin(eccentric_anomaly);
}

/* [spec: Ashby 2003 section 10, Eq. (96)] c times the Shapiro term
 * (2 GM / c^3) ln((r1 + r2 + l) / (r1 + r2 - l)). */
int32_t ssr_relativity_shapiro(const double satellite[3], const double receiver[3], double *delay)
{
    double dx = satellite[0] - receiver[0];
    double dy = satellite[1] - receiver[1];
    double dz = satellite[2] - receiver[2];
    double r1 = sqrt(satellite[0] * satellite[0] + satellite[1] * satellite[1] + satellite[2] * satellite[2]);
    double r2 = sqrt(receiver[0] * receiver[0] + receiver[1] * receiver[1] + receiver[2] * receiver[2]);
    double l = sqrt(dx * dx + dy * dy + dz * dz);
    double denominator = r1 + r2 - l;
    if (!(denominator > 0.0)) {
        return SSR_RELATIVITY_THROUGH_CENTRE;
    }
    *delay = 2.0 * SSR_WGS84_GM / (SSR_SPEED_OF_LIGHT * SSR_SPEED_OF_LIGHT) * log((r1 + r2 + l) / denominator);
    return SSR_RELATIVITY_OK;
}

/* [design: the light time is a fixed point tau = |R_z(omega_E tau) r_s -
 * r_r| / c whose map contracts by about |dE/dtau| / c, under 1e-5 for any
 * Earth orbit; from tau = |r_s - r_r| / c the first step leaves under 1e-12
 * s and the second reaches the rounding of a double, so three steps are a
 * fixed count with margin] */
#define SAGNAC_STEPS 3

/* [spec: Ashby 2003 section 2] See ssrcore/relativity.h.
 * [design: E is formed without subtracting two ranges of 2e7 m:
 * |r' - r_r|^2 - |r_s - r_r|^2 = 2 r_r . (r_s - r') since |r'| = |r_s|, so
 * E = 2 r_r . (r_s - r') / (|r' - r_r| + |r_s - r_r|), and r_s - r' is
 * taken with 1 - cos a = 2 sin^2(a / 2)] */
int32_t ssr_relativity_sagnac(const double satellite[3], const double receiver[3], double *delta)
{
    double dx = satellite[0] - receiver[0];
    double dy = satellite[1] - receiver[1];
    double dz = satellite[2] - receiver[2];
    double range = sqrt(dx * dx + dy * dy + dz * dz);
    double rotated = range;
    double sx = 0.0, sy = 0.0;
    int32_t step;

    if (!(range > 0.0)) {
        return SSR_RELATIVITY_ZERO_RANGE;
    }
    for (step = 0; step < SAGNAC_STEPS; step++) {
        double a = SSR_WGS84_OMEGA_E * rotated / SSR_SPEED_OF_LIGHT;
        double h = sin(0.5 * a);
        double versine = 2.0 * h * h;
        double s = sin(a);
        double ex, ey;
        /* r_s - r' for r' = R_z(a) r_s */
        sx = versine * satellite[0] - s * satellite[1];
        sy = versine * satellite[1] + s * satellite[0];
        ex = dx - sx;
        ey = dy - sy;
        rotated = sqrt(ex * ex + ey * ey + dz * dz);
    }
    *delta = 2.0 * (receiver[0] * sx + receiver[1] * sy) / (rotated + range);
    return SSR_RELATIVITY_OK;
}
