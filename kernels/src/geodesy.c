/* Geodesy: ECEF and ellipsoidal coordinates, local ENU frame, elevation
 * and azimuth. See ssrcore/geodesy.h. */

#include <math.h>

#include "ssrcore/geodesy.h"

/* [design: 1e-14 rad is 0.06 micrometres on the Earth's surface, far below
 * the 0.1 mm target of DESIGN 4, and some 30 units in the last place of a
 * latitude near pi/2, so rounding cannot keep the iteration from reaching
 * it] */
#define LATITUDE_TOLERANCE 1e-14

/* [design: the iteration contracts by a factor of at most about e^2 per
 * step (0.0067 for WGS-84), so a few steps reach the tolerance for any
 * point from the surface to beyond the Moon (at most 4 for WGS-84 heights
 * from -6 km to 400 000 km at every latitude); 16 leaves a wide margin and
 * bounds the time a microcontroller spends] */
#define MAX_ITERATIONS 16

/* [spec: TM-23 B.3] e^2 = 2f - f^2. */
double ssr_geodesy_eccentricity_squared(double f)
{
    return 2.0 * f - f * f;
}

/* [spec: TM-23 B.1.1] x = (N + h) cos(lat) cos(lon),
 * y = (N + h) cos(lat) sin(lon), z = ((1 - e^2) N + h) sin(lat), with
 * N = a / sqrt(1 - e^2 sin^2(lat)) (B.2). */
void ssr_geodesy_geodetic_to_ecef(double lat, double lon, double h, double a,
                                  double f, double ecef[3])
{
    const double e2 = ssr_geodesy_eccentricity_squared(f);
    const double sin_lat = sin(lat);
    const double cos_lat = cos(lat);
    const double n = a / sqrt(1.0 - e2 * sin_lat * sin_lat);

    ecef[0] = (n + h) * cos_lat * cos(lon);
    ecef[1] = (n + h) * cos_lat * sin(lon);
    ecef[2] = ((1.0 - e2) * n + h) * sin_lat;
}

/* [design: height as the projection of the position on the ellipsoid
 * normal, h = p cos(lat) + z sin(lat) - a sqrt(1 - e^2 sin^2(lat)), which
 * follows from B.1 and equals the h of B.6, p / cos(lat) - N, wherever
 * cos(lat) is not zero. B.6 divides by cos(lat), which vanishes at the
 * poles; this form holds at every latitude, and its error is of second
 * order in the latitude's, because its derivative with respect to latitude
 * vanishes on the normal] */
static double height(double p, double z, double sin_lat, double cos_lat,
                     double a, double e2)
{
    return p * cos_lat + z * sin_lat - a * sqrt(1.0 - e2 * sin_lat * sin_lat);
}

/* [spec: TM-23 B.1.2] lon = arctan(y / x) (B.4); lat(0) =
 * arctan((z / p) / (1 - e^2)) (B.5); then N(i), h(i) from lat(i-1) and
 * lat(i) = arctan((z / p) / (1 - e^2 N(i) / (N(i) + h(i)))) (B.6), until two
 * successive latitudes differ by less than the precision required.
 * [design: the arctangents of B.4 to B.6 are taken with atan2 of the
 * numerator and denominator, which gives the quadrant of the longitude and
 * lat = +-pi/2 on the axis (p = 0, where the longitude is 0); h is taken
 * in the form of height(), and computed once more from the final latitude
 * so that (lat, h) belong together] */
int32_t ssr_geodesy_ecef_to_geodetic(const double ecef[3], double a, double f,
                                     double geodetic[3])
{
    const double e2 = ssr_geodesy_eccentricity_squared(f);
    const double x = ecef[0];
    const double y = ecef[1];
    const double z = ecef[2];
    const double p = sqrt(x * x + y * y);
    double lat = atan2(z, p * (1.0 - e2));
    int32_t i;

    for (i = 0; i < MAX_ITERATIONS; ++i) {
        const double sin_lat = sin(lat);
        const double cos_lat = cos(lat);
        const double n = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
        const double h = height(p, z, sin_lat, cos_lat, a, e2);
        const double next = atan2(z, p * (1.0 - e2 * n / (n + h)));
        const double change = fabs(next - lat);

        lat = next;
        if (change < LATITUDE_TOLERANCE) {
            const double h_final = height(p, z, sin(lat), cos(lat), a, e2);
            if (!isfinite(h_final)) {
                return SSR_GEODESY_NO_SOLUTION;
            }
            geodetic[0] = lat;
            geodetic[1] = atan2(y, x);
            geodetic[2] = h_final;
            return SSR_GEODESY_OK;
        }
    }
    return SSR_GEODESY_NO_SOLUTION;
}

/* [spec: TM-23 B.2.2] (de, dn, du) = R1[pi/2 - lat] R3[pi/2 + lon] d, whose
 * rows (B.11) are e = (-sin lon, cos lon, 0),
 * n = (-cos lon sin lat, -sin lon sin lat, cos lat),
 * u = (cos lon cos lat, sin lon cos lat, sin lat) (B.9). */
void ssr_geodesy_ecef_to_enu(double lat, double lon, const double d[3],
                             double enu[3])
{
    const double sin_lat = sin(lat);
    const double cos_lat = cos(lat);
    const double sin_lon = sin(lon);
    const double cos_lon = cos(lon);

    enu[0] = -sin_lon * d[0] + cos_lon * d[1];
    enu[1] = -cos_lon * sin_lat * d[0] - sin_lon * sin_lat * d[1]
             + cos_lat * d[2];
    enu[2] = cos_lon * cos_lat * d[0] + sin_lon * cos_lat * d[1]
             + sin_lat * d[2];
}

/* [spec: TM-23 B.2.1] d = R3[-(pi/2 + lon)] R1[-(pi/2 - lat)] (de, dn, du),
 * whose columns (B.8) are the unit vectors e, n, u of B.9. */
void ssr_geodesy_enu_to_ecef(double lat, double lon, const double enu[3],
                             double d[3])
{
    const double sin_lat = sin(lat);
    const double cos_lat = cos(lat);
    const double sin_lon = sin(lon);
    const double cos_lon = cos(lon);

    d[0] = -sin_lon * enu[0] - cos_lon * sin_lat * enu[1]
           + cos_lon * cos_lat * enu[2];
    d[1] = cos_lon * enu[0] - sin_lon * sin_lat * enu[1]
           + sin_lon * cos_lat * enu[2];
    d[2] = cos_lat * enu[1] + sin_lat * enu[2];
}

/* [spec: TM-23 B.3] rho = (r_sat - r_rcv) / |r_sat - r_rcv| (B.13), and
 * rho.e = cos E sin A, rho.n = cos E cos A, rho.u = sin E (B.14).
 * [design: E and A are taken with atan2 from all three relations of B.14,
 * E = atan2(rho.u, sqrt((rho.e)^2 + (rho.n)^2)) and A = atan2(rho.e, rho.n),
 * rather than arcsin and arctan of B.15 and B.16: the same angles, with
 * the azimuth's quadrant, and no arcsin of a rounded value above 1 near
 * the zenith. A is wrapped to [0, 2 pi), the range of a bearing clockwise
 * from north (Fig. B.3)] */
int32_t ssr_geodesy_line_of_sight(const double satellite[3],
                                  const double receiver[3], double lat,
                                  double lon, double unit[3],
                                  double *elevation, double *azimuth)
{
    double d[3];
    double enu[3];
    double range;
    double az;

    d[0] = satellite[0] - receiver[0];
    d[1] = satellite[1] - receiver[1];
    d[2] = satellite[2] - receiver[2];
    range = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (!(range > 0.0)) {
        return SSR_GEODESY_ZERO_RANGE;
    }
    d[0] /= range;
    d[1] /= range;
    d[2] /= range;
    ssr_geodesy_ecef_to_enu(lat, lon, d, enu);

    az = atan2(enu[0], enu[1]);
    if (az < 0.0) {
        /* [design: 2 SSR_PI, so that azimuth wraps where atan2, which
         * works with the double nearest pi, does] */
        az += 2.0 * SSR_PI;
        if (az >= 2.0 * SSR_PI) {
            az = 0.0;  /* a negative azimuth smaller than half an ulp of 2 pi */
        }
    }
    unit[0] = d[0];
    unit[1] = d[1];
    unit[2] = d[2];
    *elevation = atan2(enu[2], sqrt(enu[0] * enu[0] + enu[1] * enu[1]));
    *azimuth = az;
    return SSR_GEODESY_OK;
}

/* [spec: TM-23 B.3, eq. B.13; design: the rate for a receiver at rest] See
 * ssrcore/geodesy.h. */
int32_t ssr_geodesy_range(const double satellite[3], const double velocity[3],
                          const double receiver[3], double *range,
                          double *rate)
{
    double d0 = satellite[0] - receiver[0];
    double d1 = satellite[1] - receiver[1];
    double d2 = satellite[2] - receiver[2];
    double r = sqrt(d0 * d0 + d1 * d1 + d2 * d2);

    if (!(r > 0.0)) {
        return SSR_GEODESY_ZERO_RANGE;
    }
    *range = r;
    *rate = (d0 * velocity[0] + d1 * velocity[1] + d2 * velocity[2]) / r;
    return SSR_GEODESY_OK;
}

/* [spec: Ashby 2003 section 2] See ssrcore/geodesy.h. */
void ssr_geodesy_rotate_z(const double v[3], double angle, double out[3])
{
    double c = cos(angle);
    double s = sin(angle);
    double x = v[0];
    double y = v[1];

    out[0] = x * c + y * s;
    out[1] = -x * s + y * c;
    out[2] = v[2];
}
