/* Carrier phase wind-up and the satellite antenna frame.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers). Vectors are
 * caller-owned double[3] in ECEF metres (positions), metres per second
 * (velocities) or unit vectors; angles in radians; wind-up in cycles.
 */
#ifndef SSRCORE_WINDUP_H
#define SSRCORE_WINDUP_H

#include <stdint.h>

/* Status codes. */
#define SSR_WINDUP_OK 0
/* [design: no frame or dipole can be built] A vector the formula
 * normalises has zero length: position at the origin, velocity parallel to
 * the position, Sun along the satellite's radius, or an effective dipole
 * that vanishes because the line of sight lies in the antenna plane. */
#define SSR_WINDUP_DEGENERATE 1

/* [spec: SPARTN v2.0.3 section 8.14, figure] Satellite-fixed antenna frame
 * used for yaw:
 *   e_z = -r/|r|,  e_y = (e_z x v)/|e_z x v|,  e_x = e_y x e_z,
 * with r the position in ECEF and v the velocity in ECI. The ECI velocity
 * is formed from the ECEF velocity as v_ecef + omega_e x r, expressed in
 * the ECEF axes of the same instant, so that r and v are in one frame
 * [design: the figure names two frames for r and v; at the instant of
 * evaluation the inertial velocity in Earth-fixed axes is v + omega_e x r,
 * omega_e = (0, 0, SSR_WGS84_OMEGA_E) for every system, since SPARTN
 * 8.2 gives every orbit in one Earth-fixed frame; the ICDs' rates differ
 * by 2e-8 relative, which turns the frame by under 1e-8 rad].
 * The frame is then rotated by the yaw angle about e_z [spec: SPARTN
 * Table 6.2 SF021, 8.14; a zero angle when yaw is not provided]; the sense
 * is right-handed about e_z, x' = cos(yaw) e_x + sin(yaw) e_y
 * [design: SPARTN does not state the sense; a positive angle by the
 * right-hand rule about the axis is the usual sense of a rotation angle;
 * DESIGN Q23].
 * Writes x', y', z (= e_z). Returns SSR_WINDUP_OK or SSR_WINDUP_DEGENERATE. */
int32_t ssr_windup_frame_spartn(const double r[3], const double v_ecef[3],
                                double yaw, double x[3], double y[3],
                                double z[3]);

/* [spec: ESA GNSS Data Processing Vol. I (TM-23/1) section 5.6.3, eqs.
 * (5.76) to (5.78)] Nominal satellite-fixed frame from the Sun:
 *   k = -r/|r|,  j = unit(k x e),  e = (r_sun - r)/|r_sun - r|,  i = j x k.
 * Writes i, j, k. Returns SSR_WINDUP_OK or SSR_WINDUP_DEGENERATE. */
int32_t ssr_windup_frame_sun(const double r[3], const double r_sun[3],
                             double i[3], double j[3], double k[3]);

/* [spec: TM-23/1 section 5.5, eqs. (5.72) and (5.74), after Wu et al. 1993]
 * Fractional wind-up, in cycles in (-1/2, 1/2], of the line of sight from
 * the satellite at r_sat to the receiver at r_rcv, for receiver dipoles a, b
 * (east and north for the local frame) and satellite dipoles a_s, b_s (x
 * and y of the satellite frame):
 *   d  = a   - rho (rho . a)   + rho x b
 *   d' = a_s - rho (rho . a_s) - rho x b_s
 *   dphi = sign(rho . (d' x d)) arccos(d' . d / (|d'| |d|)).
 * Evaluated as atan2(rho . (d' x d), d' . d) [design: d and d' are both
 * perpendicular to rho, so rho . (d' x d) = sign * |d'| |d| sin(dphi) and
 * the two forms are equal; atan2 keeps full precision near 0 and 1/2 cycle
 * where arccos does not; at exactly half a cycle, where sign() is 0 and the
 * literal form gives 0, it gives +1/2, the end of the interval].
 * Returns SSR_WINDUP_OK or SSR_WINDUP_DEGENERATE. */
int32_t ssr_windup_cycles(const double r_sat[3], const double r_rcv[3],
                          const double a[3], const double b[3],
                          const double a_s[3], const double b_s[3],
                          double *cycles);

/* [spec: TM-23/1 section 5.5, eqs. (5.71) and (5.73)] Wind-up made
 * continuous with the previous value: current + N, N = nint(previous -
 * current), all in cycles. Ties of nint are rounded away from zero, as C99
 * round() does [design: the specification does not resolve ties; a tie
 * means a step of exactly half a cycle, which either choice explains]. */
double ssr_windup_unwrap(double previous, double current);

#endif
