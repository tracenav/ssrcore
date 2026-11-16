/* Carrier phase wind-up and the satellite antenna frame. See ssrcore/windup.h. */

#include <math.h>

#include "ssrcore/geodesy.h"
#include "ssrcore/windup.h"

/* [design: vector arithmetic] u . v. */
static double dot3(const double u[3], const double v[3])
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

/* [design: vector arithmetic] out = u x v; out may not alias u or v. */
static void cross3(const double u[3], const double v[3], double out[3])
{
    out[0] = u[1] * v[2] - u[2] * v[1];
    out[1] = u[2] * v[0] - u[0] * v[2];
    out[2] = u[0] * v[1] - u[1] * v[0];
}

/* [design: vector arithmetic] out = s * u / |u|; 0 if |u| is zero, else 1. */
static int32_t unit3(const double u[3], double s, double out[3])
{
    double n = sqrt(dot3(u, u));
    if (n == 0.0) {
        return 0;
    }
    out[0] = s * u[0] / n;
    out[1] = s * u[1] / n;
    out[2] = s * u[2] / n;
    return 1;
}

/* [spec: SPARTN v2.0.3 section 8.14] See ssrcore/windup.h. */
int32_t ssr_windup_frame_spartn(const double r[3], const double v_ecef[3],
                                double yaw, double x[3], double y[3],
                                double z[3])
{
    double v[3], ezv[3], ex[3], ey[3], c, s;
    int32_t n;

    if (!unit3(r, -1.0, z)) {
        return SSR_WINDUP_DEGENERATE;
    }
    /* v_eci in ECEF axes: v_ecef + (0, 0, omega_E) x r. */
    v[0] = v_ecef[0] - SSR_WGS84_OMEGA_E * r[1];
    v[1] = v_ecef[1] + SSR_WGS84_OMEGA_E * r[0];
    v[2] = v_ecef[2];
    cross3(z, v, ezv);
    if (!unit3(ezv, 1.0, ey)) {
        return SSR_WINDUP_DEGENERATE;
    }
    cross3(ey, z, ex);
    c = cos(yaw);
    s = sin(yaw);
    for (n = 0; n < 3; n++) {
        x[n] = c * ex[n] + s * ey[n];
        y[n] = c * ey[n] - s * ex[n];
    }
    return SSR_WINDUP_OK;
}

/* [spec: TM-23/1 section 5.6.3, eqs. (5.76) to (5.78)] See ssrcore/windup.h. */
int32_t ssr_windup_frame_sun(const double r[3], const double r_sun[3],
                             double i[3], double j[3], double k[3])
{
    double d[3], e[3], ke[3];

    if (!unit3(r, -1.0, k)) {
        return SSR_WINDUP_DEGENERATE;
    }
    d[0] = r_sun[0] - r[0];
    d[1] = r_sun[1] - r[1];
    d[2] = r_sun[2] - r[2];
    if (!unit3(d, 1.0, e)) {
        return SSR_WINDUP_DEGENERATE;
    }
    cross3(k, e, ke);
    if (!unit3(ke, 1.0, j)) {
        return SSR_WINDUP_DEGENERATE;
    }
    cross3(j, k, i);
    return SSR_WINDUP_OK;
}

/* [spec: TM-23/1 section 5.5, eqs. (5.72), (5.74)] See ssrcore/windup.h. */
int32_t ssr_windup_cycles(const double r_sat[3], const double r_rcv[3],
                          const double a[3], const double b[3],
                          const double a_s[3], const double b_s[3],
                          double *cycles)
{
    double los[3], rho[3], rb[3], rbs[3], d[3], ds[3], dsd[3];
    double pa, pas;
    int32_t n;

    los[0] = r_rcv[0] - r_sat[0];
    los[1] = r_rcv[1] - r_sat[1];
    los[2] = r_rcv[2] - r_sat[2];
    if (!unit3(los, 1.0, rho)) {
        return SSR_WINDUP_DEGENERATE;
    }
    pa = dot3(rho, a);
    pas = dot3(rho, a_s);
    cross3(rho, b, rb);
    cross3(rho, b_s, rbs);
    for (n = 0; n < 3; n++) {
        d[n] = a[n] - rho[n] * pa + rb[n];
        ds[n] = a_s[n] - rho[n] * pas - rbs[n];
    }
    if (dot3(d, d) == 0.0 || dot3(ds, ds) == 0.0) {
        return SSR_WINDUP_DEGENERATE;
    }
    cross3(ds, d, dsd);
    *cycles = atan2(dot3(rho, dsd), dot3(ds, d)) / (2.0 * SSR_PI);
    /* atan2(-0, negative) is -pi: fold it to +1/2, the interval's end. */
    if (*cycles <= -0.5) {
        *cycles += 1.0;
    }
    return SSR_WINDUP_OK;
}

/* [spec: TM-23/1 section 5.5, eqs. (5.71), (5.73)] See ssrcore/windup.h. */
double ssr_windup_unwrap(double previous, double current)
{
    return current + round(previous - current);
}
