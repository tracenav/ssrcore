/* State-space corrections of the satellite orbit and clock. See ssrcore/correction.h. */

#include <math.h>

#include "ssrcore/correction.h"

/* [design: arithmetic] a x b. */
static void cross3(const double a[3], const double b[3], double out[3])
{
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

/* [design: arithmetic] |a|, as sqrt(x x + y y + z z). */
static double norm3(const double a[3])
{
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* [spec: SPARTN v2.0.3 section 8.8] See correction.h. Each unit vector is its
 * vector divided by its norm, as 8.8 writes it. [design: a zero norm is
 * reported, not divided by] */
int32_t ssr_correction_orbit_frame(const double r[3], const double r_dot[3],
                                   double e_radial[3], double e_along[3], double e_cross[3])
{
    double rc[3];
    double nv, nc;
    double a[3], c[3];
    int32_t i;

    cross3(r, r_dot, rc);
    nv = norm3(r_dot);
    nc = norm3(rc);
    if (!(nv > 0.0) || !(nc > 0.0)) {
        return SSR_CORRECTION_DEGENERATE;
    }
    for (i = 0; i < 3; ++i) {
        a[i] = r_dot[i] / nv;
        c[i] = rc[i] / nc;
    }
    cross3(a, c, e_radial);
    for (i = 0; i < 3; ++i) {
        e_along[i] = a[i];
        e_cross[i] = c[i];
    }
    return SSR_CORRECTION_OK;
}

/* [spec: SPARTN v2.0.3 section 8.8] See correction.h. */
int32_t ssr_correction_orbit(const double r[3], const double r_dot[3],
                             const double d_o[3], double d_x[3])
{
    double er[3], ea[3], ec[3];
    int32_t i;
    int32_t status = ssr_correction_orbit_frame(r, r_dot, er, ea, ec);
    if (status != SSR_CORRECTION_OK) {
        return status;
    }
    for (i = 0; i < 3; ++i) {
        d_x[i] = er[i] * d_o[0] + ea[i] * d_o[1] + ec[i] * d_o[2];
    }
    return SSR_CORRECTION_OK;
}

/* [spec: SPARTN v2.0.3 section 8.8] See correction.h. */
int32_t ssr_correction_position(const double r[3], const double r_dot[3],
                                const double d_o[3], double x[3])
{
    double d_x[3];
    int32_t i;
    int32_t status = ssr_correction_orbit(r, r_dot, d_o, d_x);
    if (status != SSR_CORRECTION_OK) {
        return status;
    }
    for (i = 0; i < 3; ++i) {
        x[i] = r[i] - d_x[i];
    }
    return SSR_CORRECTION_OK;
}

/* [spec: SPARTN v2.0.3 section 8.7, Table 8.0 for c] See correction.h. */
double ssr_correction_clock(double t_b, double d_c)
{
    return t_b - d_c / SSR_SPEED_OF_LIGHT;
}

/* [spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6); IGS SSR v1.00 4.1] See correction.h. */
void ssr_correction_orbit_at(const double d_o[3], const double rate[3], double dt, double out[3])
{
    int32_t i;

    for (i = 0; i < 3; ++i) {
        out[i] = d_o[i] + rate[i] * dt;
    }
}

/* [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8); IGS SSR v1.00 4.2] See correction.h.
 * Evaluated in the printed order, C0 + C1 dt + C2 dt dt. */
double ssr_correction_clock_at(double c0, double c1, double c2, double dt)
{
    return c0 + c1 * dt + c2 * dt * dt;
}

/* [spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8), its derivative] See correction.h. */
double ssr_correction_clock_rate_at(double c1, double c2, double dt)
{
    return c1 + 2.0 * c2 * dt;
}
