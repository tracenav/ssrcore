/* UNB3m a-priori neutral atmosphere model: the meteorological parameters of
 * its look-up table and the hydrostatic and non-hydrostatic (wet) zenith
 * delays, after Leandro, Santos and Langley (2006), "UNB Neutral Atmosphere
 * Models: Development and Performance" (UNB3m).
 *
 * The temperature T0 and lapse rate beta of the same table are the
 * troposphere kernel's (ssr_troposphere_unb3m_temperature), which the caller
 * passes in [design: DESIGN 3.1, a kernel that needs another's result takes
 * it as an argument; one definition of the shared table rows]. This kernel
 * holds the rows only it needs: pressure and lambda (UNB3m Table 1) and
 * relative humidity (UNB3m Table 2).
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers). Latitude in
 * radians, heights in metres, pressures in mbar (hPa), temperatures in K,
 * day of year from January 0.0 (1.0 at 00:00 of 1 January), as DESIGN 2.7
 * states for the troposphere kernel.
 */
#ifndef SSRCORE_UNB3M_H
#define SSRCORE_UNB3M_H

#include <stdint.h>

/* Status codes. [design: DESIGN 3.1, errors as int32_t status codes] */
#define SSR_UNB3M_OK 0
/* height at or above T0/beta, where the lapse-rate profile of UNB3m eq. (4) ends */
#define SSR_UNB3M_HEIGHT 1
/* temperature or lapse rate not positive, or lambda' g_m - beta R not positive (UNB3m eqs. (4), (5), (7)) */
#define SSR_UNB3M_PARAMETERS 2

/* [spec: UNB3m eqs. (1) to (3), Tables 1 and 2; TM-23 5.4.2.1 eq. (5.62) for
 * |latitude| and the southern D_min of 211] Sea-level pressure P0 (mbar),
 * relative humidity RH (%) and water vapour lapse rate lambda (unitless) at a
 * latitude (rad) and day of year. */
void ssr_unb3m_parameters(double latitude, double doy, double *pressure, double *humidity, double *lambda);

/* [spec: UNB3m eqs. (9) to (11)] Water vapour pressure e0 (mbar) from
 * relative humidity RH (%), temperature T0 (K) and pressure P0 (mbar):
 * e0 = RH/100 * e_s(T0) * f_w(P0, T0). */
double ssr_unb3m_vapour_pressure(double humidity, double temperature, double pressure);

/* [spec: UNB3m eqs. (4) to (7); TM-23 5.4.2.1 eq. (5.64) for g] Hydrostatic
 * d_zh and non-hydrostatic d_znh zenith delays (m) at height H (m) for the
 * sea-level parameters P0 (mbar), T0 (K), e0 (mbar), beta (K/m) and lambda,
 * with g_m of eq. (6) at the latitude (rad) and H. Returns SSR_UNB3M_HEIGHT
 * when 1 - beta H / T0 <= 0, SSR_UNB3M_PARAMETERS when T0 <= 0, beta <= 0
 * or (lambda + 1) g_m - beta R <= 0; the outputs are then untouched. */
int32_t ssr_unb3m_zenith(double latitude, double height, double pressure, double temperature,
                         double vapour_pressure, double lapse_rate, double lambda,
                         double *hydrostatic, double *wet);

#endif
