/* Ionosphere: the slant delay of a signal in metres from slant TEC (SPARTN
 * v2.0.3 8.11), and the slant TEC of an IGS SSR VTEC layer (IGS SSR v1.00
 * 4.5.1). The slant TEC of a SPARTN area polynomial and grid residual are
 * ssr_area_polynomial and ssr_area_interpolate.
 *
 * TEC values are in TEC units (TECU, 10^16 electrons per square metre), the
 * unit of SF057 to SF067.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 */
#ifndef SSRCORE_IONOSPHERE_H
#define SSRCORE_IONOSPHERE_H

#include <stdint.h>

/* Status codes. [design: DESIGN 3.1, errors as int32_t status codes] */
#define SSR_IONOSPHERE_OK 0
#define SSR_IONOSPHERE_BAD_FREQUENCY 1 /* frequency not positive */
#define SSR_IONOSPHERE_BAD_POSITION 2  /* receiver at the centre, or no line of sight */
#define SSR_IONOSPHERE_BELOW_LAYER 3   /* no pierce point: layer below the receiver, or E + psi not in (0, pi) */
#define SSR_IONOSPHERE_BAD_DEGREE 4    /* degree or order outside 0 to 16 */

/* [spec: SPARTN v2.0.3 8.11] Slant ionosphere delay, metres, of a signal of
 * frequency f (Hz) for slant TEC I_P + I_G (TECU):
 *   pseudorange   *code  = +40.3e16 / f^2 (I_P + I_G)
 *   carrier phase *phase = -40.3e16 / f^2 (I_P + I_G)
 * Pass tec_grid = 0 when no grid is present (8.11). */
int32_t ssr_ionosphere_delay(double tec_polynomial, double tec_grid,
                             double frequency, double *code, double *phase);

/* IGS SSR VTEC (IGS SSR v1.00 4.5.1, message IM201): a thin-layer vertical TEC
 * of spherical harmonics on a spherical Earth of radius 6370 km, mapped to
 * slant TEC at the layer's pierce point. Angles in radians, heights in
 * metres, TEC in TECU. */

/* [spec: IGS SSR v1.00 4.5.1, "Latitude, longitude, height as well as
 * azimuth and elevation of the satellites are defined with respect to this
 * model"] The receiver (ECEF, m) and the unit line of sight to the satellite
 * (ECEF, in the frame of reception) in the spherical model:
 * out = {geocentric latitude, longitude, height above the 6370 km sphere,
 * elevation, azimuth clockwise from north in [0, 2 pi)}. */
int32_t ssr_ionosphere_spherical(const double receiver[3], const double unit[3],
                                 double out[5]);

/* [spec: IGS SSR v1.00 4.5.1] The pierce point of a layer of height
 * layer_height above the sphere, seen from a receiver at geocentric latitude,
 * longitude and height in the spherical model along elevation and azimuth:
 * out = {latitude, longitude, central angle psi, slant factor 1/sin(E + psi)}. */
int32_t ssr_ionosphere_pierce(double latitude, double longitude, double height,
                              double elevation, double azimuth,
                              double layer_height, double out[4]);

/* [spec: IGS SSR v1.00 4.5.1, 8.4.1] VTEC (TECU) of one layer at a pierce
 * point (latitude, longitude) for t, the seconds of day: degree N and order M
 * (each 0 to 16, M taken as min(M, N)), cosine C(n,m) for m = 0..M, n = m..N
 * and sine S(n,m) for m = 1..M, n = m..N, in that order (IGS SSR 8.4.1).
 * Negative values are returned as they are; the caller replaces them by 0. */
int32_t ssr_ionosphere_vtec(int32_t degree, int32_t order, const double *cosine,
                            const double *sine, double latitude, double longitude,
                            double t, double *vtec);

#endif
