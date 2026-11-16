/* Geodesy: Cartesian (ECEF) and ellipsoidal coordinates, the local
 * East-North-Up frame, elevation and azimuth.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * Angles are in radians, lengths in metres. Ellipsoidal coordinates are
 * (latitude, longitude, height) in that order. The ellipsoid is passed as
 * its semi-major axis a and flattening f = 1 - b/a, the two parameters
 * TM-23 3.1.4 names, so that the same functions serve WGS-84 and the other
 * ellipsoids of TM-23 3.1.3.
 */
#ifndef SSRCORE_GEODESY_H
#define SSRCORE_GEODESY_H

#include <stdint.h>

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV "WGS 84 Earth
 * Equatorial Radius"; TM-23 Table 3.2] Semi-major axis a of the WGS-84
 * ellipsoid, m. */
#define SSR_WGS84_A 6378137.0

/* [spec: TM-23 Table 3.2] Flattening f of the WGS-84 ellipsoid. */
#define SSR_WGS84_F (1.0 / 298.257223563)

/* [spec: TM-23 Table 3.2; Ashby 2003 section 3, Eq. (13), WGS-84 (G873)]
 * Geocentric gravitational constant GM of WGS-84, m^3/s^2. The orbit
 * kernels use the mu of each system's ICD for its own algorithm (DESIGN 4). */
#define SSR_WGS84_GM 3.986004418e14

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV, "WGS 84 value
 * of the earth's rotation rate"; the value of Ashby 2003 section 2] Rotation
 * rate of the ECEF frame, rad/s. TM-23 Table 3.2 and the other ICDs print
 * 7.2921150e-5 to 7.292115e-5; the relative difference, 2e-8, is under a
 * micrometre in every term that uses this constant (DESIGN 4). */
#define SSR_WGS84_OMEGA_E 7.2921151467e-5

/* [design: the double nearest pi, the scale of the angles the C library's
 * trigonometric functions take and return; SPARTN Table 8.0 and IS-GPS-200N
 * with IRN-IS-200N-003 and -004 print 3.1415926535898, which only the orbit
 * kernels use, for the ICD curve fit (DESIGN 4)] */
#define SSR_PI 3.141592653589793

/* Status codes. */
#define SSR_GEODESY_OK 0
/* [design: the iteration of TM-23 B.1.2 has no end condition other than a
 * precision; a point where it does not reach that precision in a bounded
 * number of steps, or reaches it with a height that is not finite (the
 * centre of the ellipsoid, non-finite input), is reported, not returned
 * as a value] */
#define SSR_GEODESY_NO_SOLUTION 1
/* [design: B.13 divides by the distance; coincident points have no line
 * of sight] */
#define SSR_GEODESY_ZERO_RANGE 2

/* [spec: TM-23 B.3] Squared first eccentricity e^2 = 2f - f^2. */
double ssr_geodesy_eccentricity_squared(double f);

/* [spec: TM-23 B.1.1, eqs. B.1 to B.3] ECEF position of the ellipsoidal
 * coordinates (lat, lon, h) on the ellipsoid (a, f). */
void ssr_geodesy_geodetic_to_ecef(double lat, double lon, double h, double a,
                                  double f, double ecef[3]);

/* [spec: TM-23 B.1.2, eqs. B.4 to B.6] Ellipsoidal coordinates
 * (lat, lon, h) of an ECEF position on the ellipsoid (a, f). Returns
 * SSR_GEODESY_OK, or SSR_GEODESY_NO_SOLUTION with geodetic unchanged. */
int32_t ssr_geodesy_ecef_to_geodetic(const double ecef[3], double a, double f,
                                     double geodetic[3]);

/* [spec: TM-23 B.2.2, eq. B.10 with the matrix B.11] East, north and up
 * components of an ECEF vector d, in the local frame at (lat, lon). */
void ssr_geodesy_ecef_to_enu(double lat, double lon, const double d[3],
                             double enu[3]);

/* [spec: TM-23 B.2.1, eq. B.7 with the matrix B.8] ECEF components of a
 * vector given by its east, north and up components at (lat, lon). */
void ssr_geodesy_enu_to_ecef(double lat, double lon, const double enu[3],
                             double d[3]);

/* [spec: TM-23 B.3, eqs. B.13 and B.14] Unit vector from the receiver to
 * the satellite (ECEF), elevation in [-pi/2, pi/2] and azimuth in
 * [0, 2 pi), clockwise from north, in the local frame at the receiver's
 * ellipsoidal latitude and longitude (lat, lon). The satellite comes first,
 * as in every function of a line of sight (DESIGN 3.1). Returns
 * SSR_GEODESY_OK, or SSR_GEODESY_ZERO_RANGE with the outputs unchanged. */
int32_t ssr_geodesy_line_of_sight(const double satellite[3],
                                  const double receiver[3], double lat,
                                  double lon, double unit[3],
                                  double *elevation, double *azimuth);

/* [spec: TM-23 B.3, eq. B.13: rho = |r_s - r_r|, the norm that makes the
 * unit vector; DESIGN 2.4, row rho] Range (m) from the receiver to the
 * satellite, both ECEF, and its rate (m/s) u . v_s, the derivative of rho
 * for the satellite velocity v_s and a receiver at rest in the frame
 * [design: the query has the receiver's position only; a moving receiver
 * adds -u . v_r, which is not an input]. Returns SSR_GEODESY_OK, or
 * SSR_GEODESY_ZERO_RANGE with the outputs unchanged. */
int32_t ssr_geodesy_range(const double satellite[3], const double velocity[3],
                          const double receiver[3], double *range,
                          double *rate);

/* [spec: Ashby 2003 section 2: the frame turns by omega_E tau during the
 * travel time, R_z(a) with x' = x cos a + y sin a, y' = -x sin a + y cos a,
 * z' = z] The coordinates of v in a frame turned by angle a (rad) about z. */
void ssr_geodesy_rotate_z(const double v[3], double angle, double out[3]);

#endif
