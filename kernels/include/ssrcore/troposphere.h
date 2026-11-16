/* Troposphere: Niell (1996) mapping functions as TM-23 5.4.2.2.1 gives them,
 * the SPARTN 8.12 slant delay (its polynomial and grid residuals are
 * ssr_area_polynomial and ssr_area_interpolate), and
 * the UNB3m height reduction of the hydrostatic zenith delay (DESIGN Q7).
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers). Angles in radians,
 * heights in metres, day of year as days from January 0.0 (1.0 at 00:00 of
 * 1 January), all as DESIGN section 4 states.
 */
#ifndef SSRCORE_TROPOSPHERE_H
#define SSRCORE_TROPOSPHERE_H

#include <stdint.h>

/* Status codes. [design: DESIGN 3.1, errors as int32_t status codes] */
#define SSR_TROPOSPHERE_OK 0
/* elevation not in (0, pi/2] */
#define SSR_TROPOSPHERE_ELEVATION 1
/* height at or above T0/beta, where the UNB3m lapse-rate profile ends */
#define SSR_TROPOSPHERE_HEIGHT 2

/* [spec: TM-23 5.4.2.2.1 eq. (5.69)] Marini (1972) continued fraction
 * normalised to unity at zenith, m(E, a, b, c), for sin E = sin_e. */
double ssr_troposphere_marini(double sin_e, double a, double b, double c);

/* [spec: TM-23 5.4.2.2.1 eq. (5.70), Tables 5.2 and 5.3] Niell hydrostatic
 * (a_d, b_d, c_d) and wet (a_w, b_w, c_w) coefficients at a latitude (rad)
 * and day of year. */
void ssr_troposphere_niell_coefficients(double latitude, double doy, double hydrostatic[3], double wet[3]);

/* [spec: TM-23 5.4.2.2.1 eqs. (5.67) to (5.70); SPARTN v2.0.3 8.12 for the
 * ellipsoidal height] Niell hydrostatic mapping M_dry(E, h), with the height
 * correction, and wet mapping M_wet(E). Returns SSR_TROPOSPHERE_ELEVATION, and
 * leaves the outputs untouched, when the elevation is not in (0, pi/2]. */
int32_t ssr_troposphere_niell_mapping(double latitude, double height, double doy, double elevation,
                                      double *hydrostatic, double *wet);

/* [spec: SPARTN v2.0.3 8.12] Slant delay
 * dT = T_h * m_h(E, h) + (T_P + T_G) * m_w(E). */
double ssr_troposphere_slant(double hydrostatic_zenith, double polynomial_residual,
                             double grid_residual, double hydrostatic_mapping,
                             double wet_mapping);

/* [spec: UNB3m eqs. (1) to (3), Table 1; TM-23 5.4.2.1 eq. (5.62), Table 5.1
 * for |latitude| and the southern phase] UNB3m surface temperature T0 (K) and
 * temperature lapse rate beta (K/m) at a latitude (rad) and day of year. */
void ssr_troposphere_unb3m_temperature(double latitude, double doy, double *temperature, double *lapse_rate);

/* [spec: UNB3m eqs. (4) and (6); TM-23 5.4.2.1 eq. (5.64) for R_d and g;
 * design: DESIGN Q7] Ratio of the hydrostatic zenith delay at height h to
 * that at height 0, for surface temperature T0 and lapse rate beta:
 * (g_m(0) / g_m(h)) * (1 - beta h / T0)^(g / (R_d beta)). Returns
 * SSR_TROPOSPHERE_HEIGHT when 1 - beta h / T0 <= 0. */
int32_t ssr_troposphere_unb3m_height_factor(double latitude, double temperature,
                                            double lapse_rate, double height, double *factor);

#endif
