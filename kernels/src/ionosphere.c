/* Ionosphere: SPARTN slant delay and IGS SSR VTEC. See ssrcore/ionosphere.h. */

#include <math.h>
#include <stdint.h>

#include "ssrcore/geodesy.h"
#include "ssrcore/ionosphere.h"

/* [spec: SPARTN v2.0.3 8.11; TM-23 5.4.1 eqs. 5.37 to 5.40] Delay of one
 * TECU at frequency f is 40.3e16 / f^2 metres (f in Hz). SPARTN 8.11 prints
 * the factor as 40.3 x 10^-16; with TEC in TECU (10^16 electrons/m^2,
 * SF057 to SF067) and f in Hz the factor is 40.3 x 10^16, as TM-23 eq. 5.40
 * and its Table 4.2 (0.1624 m/TECU at GPS L1) give. */
#define IONOSPHERE_FACTOR 40.3e16

/* [spec: SPARTN v2.0.3 8.11] delta I_PR = 40.3e16 / f^2 (I_P + I_G),
 * delta I_CP = -delta I_PR. */
int32_t ssr_ionosphere_delay(double tec_polynomial, double tec_grid,
                             double frequency, double *code, double *phase)
{
    double d;
    if (!(frequency > 0.0))
        return SSR_IONOSPHERE_BAD_FREQUENCY;
    d = IONOSPHERE_FACTOR / (frequency * frequency)
        * (tec_polynomial + tec_grid);
    *code = d;
    *phase = -d;
    return SSR_IONOSPHERE_OK;
}

/* ---- IGS SSR VTEC (IGS SSR v1.00 4.5.1) --------------------------------- */

/* [spec: IGS SSR v1.00 4.5.1, "a spherical Earth model, with a mean radius of
 * 6370 km"] metres. */
#define SPHERE_RADIUS 6370.0e3
/* [spec: IGS SSR v1.00 Table 25, IDF037 and IDF038 are 4 bits of N - 1 and
 * M - 1] the largest degree and order a message can carry. */
#define MAX_DEGREE 16
/* [spec: IGS SSR v1.00 4.5.1, "phase shifted by 2 h to the approximate TEC
 * maximum at 14:00 local time (resp. 50400 s)", pi / 43200 rad per second;
 * design: SSR_PI, the double nearest pi (DESIGN 4)] */
#define SUN_SHIFT 50400.0
#define SUN_RATE (SSR_PI / 43200.0)

/* [design: rounding can carry the argument of an arcsine just past 1] */
static double clamped_asin(double x)
{
    return asin(x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x));
}

/* [spec: IGS SSR v1.00 4.5.1, latitude, longitude, height, azimuth and
 * elevation with respect to the 6370 km sphere; design: on a sphere the
 * local vertical is the radius, so the geocentric latitude and longitude
 * give the east, north, up basis as TM-23 B.1.2 gives it for the geodetic
 * ones, and the elevation and azimuth of the unit line of sight follow as
 * in TM-23 B.3] */
int32_t ssr_ionosphere_spherical(const double receiver[3], const double unit[3],
                                 double out[5])
{
    double p = hypot(receiver[0], receiver[1]);
    double r = hypot(p, receiver[2]);
    double n = sqrt(unit[0] * unit[0] + unit[1] * unit[1] + unit[2] * unit[2]);
    double lat, lon, sl, cl, so, co, e, nn, u, az;
    if (!(r > 0.0) || !(n > 0.0))
        return SSR_IONOSPHERE_BAD_POSITION;
    lat = atan2(receiver[2], p);
    lon = atan2(receiver[1], receiver[0]);
    sl = sin(lat);
    cl = cos(lat);
    so = sin(lon);
    co = cos(lon);
    e = (-so * unit[0] + co * unit[1]) / n;
    nn = (-sl * co * unit[0] - sl * so * unit[1] + cl * unit[2]) / n;
    u = (cl * co * unit[0] + cl * so * unit[1] + sl * unit[2]) / n;
    az = atan2(e, nn);
    if (az < 0.0)
        az += 2.0 * SSR_PI;
    out[0] = lat;
    out[1] = lon;
    out[2] = r - SPHERE_RADIUS;
    out[3] = clamped_asin(u);
    out[4] = az;
    return SSR_IONOSPHERE_OK;
}

/* [spec: IGS SSR v1.00 4.5.1: psi = pi/2 - E - asin((Re + hR)/(Re + hI) cos E);
 * phi_PP = asin(sin phi_R cos psi + cos phi_R sin psi cos A);
 * lambda_PP = lambda_R + pi - asin(sin psi sin A / cos phi_PP) when
 * phi_R >= 0 and tan psi cos A > tan(pi/2 - phi_R), or phi_R < 0 and
 * -tan psi cos A > tan(pi/2 + phi_R); lambda_R + asin(...) otherwise;
 * STEC_i = VTEC / sin(E + psi). design: a layer the line of sight does not
 * cross outward (the arcsine undefined, or E + psi outside (0, pi)) has no
 * pierce point and is refused] */
int32_t ssr_ionosphere_pierce(double latitude, double longitude, double height,
                              double elevation, double azimuth,
                              double layer_height, double out[4])
{
    double ratio = (SPHERE_RADIUS + height) / (SPHERE_RADIUS + layer_height);
    double ce = cos(elevation);
    double psi, lat_pp, s;
    int32_t far;
    if (!(SPHERE_RADIUS + layer_height > 0.0) || !(ratio * ce <= 1.0))
        return SSR_IONOSPHERE_BELOW_LAYER;
    psi = SSR_PI / 2.0 - elevation - asin(ratio * ce);
    s = sin(elevation + psi);
    if (!(s > 0.0))
        return SSR_IONOSPHERE_BELOW_LAYER;
    lat_pp = clamped_asin(sin(latitude) * cos(psi) + cos(latitude) * sin(psi) * cos(azimuth));
    if (latitude >= 0.0)
        far = tan(psi) * cos(azimuth) > tan(SSR_PI / 2.0 - latitude);
    else
        far = -tan(psi) * cos(azimuth) > tan(SSR_PI / 2.0 + latitude);
    out[0] = lat_pp;
    out[2] = psi;
    out[3] = 1.0 / s;
    s = clamped_asin(sin(psi) * sin(azimuth) / cos(lat_pp));
    out[1] = far ? longitude + SSR_PI - s : longitude + s;
    return SSR_IONOSPHERE_OK;
}

/* [spec: IGS SSR v1.00 4.5.1, VTEC = sum over n = 0..N, m = 0..min(n, M) of
 * (C_nm cos m lambda_S + S_nm sin m lambda_S) P_nm(sin phi_PP), with
 * lambda_S = (lambda_PP + (t - 50400) pi / 43200) modulo 2 pi; 8.4.1 for the
 * order of the coefficients. design: the document says "fully normalized"
 * and defines no normalization; P_nm is taken as the geodetic 4 pi
 * normalization, sqrt((2 - d_m0)(2n + 1)(n - m)! / (n + m)!) times the
 * associated Legendre function without the Condon-Shortley phase. The
 * recursions follow from that definition: P_00 = 1, P_11 = sqrt(3) cos phi,
 * P_mm = sqrt((2m + 1)/(2m)) cos phi P_(m-1)(m-1), P_(m+1)m =
 * sqrt(2m + 3) sin phi P_mm, and P_nm = a sin phi P_(n-1)m - b P_(n-2)m with
 * a = sqrt((2n - 1)(2n + 1) / ((n - m)(n + m))),
 * b = sqrt((2n + 1)(n + m - 1)(n - m - 1) / ((n - m)(n + m)(2n - 3))).
 * One column of m at a time, in the order of the coefficients, so no table
 * is kept] */
int32_t ssr_ionosphere_vtec(int32_t degree, int32_t order, const double *cosine,
                            const double *sine, double latitude, double longitude,
                            double t, double *vtec)
{
    double x = sin(latitude), u = cos(latitude);
    double lambda_s, sum = 0.0, p_mm = 1.0;
    int32_t m, n, ic = 0, is = 0;
    if (degree < 0 || degree > MAX_DEGREE || order < 0 || order > MAX_DEGREE)
        return SSR_IONOSPHERE_BAD_DEGREE;
    if (order > degree)
        order = degree;
    lambda_s = fmod(longitude + (t - SUN_SHIFT) * SUN_RATE, 2.0 * SSR_PI);
    if (lambda_s < 0.0)
        lambda_s += 2.0 * SSR_PI;
    for (m = 0; m <= order; m++) {
        double cm = cos(m * lambda_s), sm = sin(m * lambda_s);
        double p1, p2 = 0.0, p;
        if (m == 1)
            p_mm = sqrt(3.0) * u;
        else if (m > 1)
            p_mm *= sqrt((2.0 * m + 1.0) / (2.0 * m)) * u;
        p1 = p_mm;
        for (n = m; n <= degree; n++) {
            if (n == m) {
                p = p_mm;
            } else if (n == m + 1) {
                p = sqrt(2.0 * m + 3.0) * x * p_mm;
            } else {
                double a = sqrt((2.0 * n - 1.0) * (2.0 * n + 1.0) / ((double)(n - m) * (n + m)));
                double b = sqrt((2.0 * n + 1.0) * (n + m - 1.0) * (n - m - 1.0)
                                / ((double)(n - m) * (n + m) * (2.0 * n - 3.0)));
                p = a * x * p1 - b * p2;
            }
            if (n > m) {
                p2 = p1;
                p1 = p;
            }
            sum += cosine[ic++] * cm * p;
            if (m > 0)
                sum += sine[is++] * sm * p;
        }
    }
    *vtec = sum;
    return SSR_IONOSPHERE_OK;
}
