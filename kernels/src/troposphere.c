/* Troposphere: Niell mapping, SPARTN slant delay, UNB3m hydrostatic height
 * reduction. See ssrcore/troposphere.h. */

#include <math.h>

#include "ssrcore/geodesy.h"
#include "ssrcore/troposphere.h"

/* [spec: TM-23 5.4.2.2.1 eq. (5.70)] Length of the seasonal period, days. */
static const double YEAR = 365.25;

/* [spec: TM-23 5.4.2.2.1 eq. (5.70)] Day of year of the seasonal extreme,
 * T0 = 28, for the Niell coefficients. */
static const double NIELL_T0 = 28.0;

/* [spec: TM-23 5.4.2.1 eq. (5.62)] D_min of the UNB3 table: 28 for northern
 * latitudes, 211 for southern latitudes. UNB3m eq. (3) gives the northern
 * value; its Table 1 is the UNB3 table that TM-23 Table 5.1 reproduces. */
static const double UNB3_DMIN_NORTH = 28.0;
static const double UNB3_DMIN_SOUTH = 211.0;

/* [spec: TM-23 5.4.2.1 eq. (5.64)] Gas constant of dry air R_d, J/(kg K),
 * also UNB3m eq. (5) text; surface gravity g, m/s^2, which UNB3m eq. (4)
 * uses and does not state. */
static const double R_DRY = 287.054;
static const double G_SURFACE = 9.80665;

/* [spec: TM-23 Table 5.2, Niell (1996)] Hydrostatic averages a, b, c and
 * amplitudes a, b, c at latitudes 15, 30, 45, 60 and 75 degrees. */
static const double NIELL_HYDROSTATIC_AVG[3][5] = {
    {1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3},
    {2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3},
    {62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3},
};
static const double NIELL_HYDROSTATIC_AMP[3][5] = {
    {0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5},
    {0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5},
    {0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5},
};

/* [spec: TM-23 Table 5.2] Height correction a_ht, b_ht, c_ht. */
static const double NIELL_HEIGHT[3] = {2.53e-5, 5.49e-3, 1.14e-3};

/* [spec: TM-23 Table 5.3, Niell (1996)] Wet a, b, c at 15 to 75 degrees. */
static const double NIELL_WET[3][5] = {
    {5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4},
    {1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3},
    {4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2},
};

/* [spec: UNB3m Table 1; TM-23 Table 5.1] Temperature T0 (K) and lapse rate
 * beta (K/m): averages, then amplitudes, at 15 to 75 degrees. */
static const double UNB3_T_AVG[5] = {299.65, 294.15, 283.15, 272.15, 263.65};
static const double UNB3_T_AMP[5] = {0.00, 7.00, 11.00, 15.00, 14.50};
static const double UNB3_BETA_AVG[5] = {6.30e-3, 6.05e-3, 5.58e-3, 5.39e-3, 4.53e-3};
static const double UNB3_BETA_AMP[5] = {0.00, 0.25e-3, 0.32e-3, 0.81e-3, 0.62e-3};

/* [spec: TM-23 Table 5.1 caption; UNB3m eqs. (1), (2)] Linear interpolation
 * of a row tabled at 15, 30, ..., 75 degrees, at |latitude| in degrees;
 * the end values hold at or below 15 and at or above 75 degrees. TM-23
 * 5.4.2.2.1 says only "linearly interpolated between the nearest" for the
 * Niell tables; the same rule is applied to them [design: the tables span
 * the same latitudes and TM-23 states no other rule outside them]. */
static double interpolate(const double row[5], double abs_lat_deg)
{
    double x;
    int32_t i;
    if (abs_lat_deg <= 15.0) {
        return row[0];
    }
    if (abs_lat_deg >= 75.0) {
        return row[4];
    }
    x = (abs_lat_deg - 15.0) / 15.0;
    i = (int32_t)x;
    if (i > 3) {
        i = 3;
    }
    return row[i] + (row[i + 1] - row[i]) * (x - (double)i);
}

/* [spec: TM-23 5.4.2.2.1 eq. (5.69)] m(E, a, b, c). */
double ssr_troposphere_marini(double sin_e, double a, double b, double c)
{
    return (1.0 + a / (1.0 + b / (1.0 + c))) / (sin_e + a / (sin_e + b / (sin_e + c)));
}

/* [spec: TM-23 5.4.2.2.1 eq. (5.70), Tables 5.2, 5.3] xi = xi_avg -
 * xi_amp cos(2 pi (t - T0) / 365.25). In the southern hemisphere the phase
 * is advanced by half a year [design: TM-23 5.4.2.2.1 states no southern
 * case; the seasons are opposite there, and TM-23 5.4.2.1 eq. (5.62) shifts
 * its seasonal table the same way]. */
void ssr_troposphere_niell_coefficients(double latitude, double doy, double hydrostatic[3], double wet[3])
{
    double lat_deg = latitude * 180.0 / SSR_PI;
    double abs_lat = fabs(lat_deg);
    double phase = (doy - NIELL_T0) / YEAR;
    double season;
    int32_t k;
    if (lat_deg < 0.0) {
        phase += 0.5;
    }
    season = cos(2.0 * SSR_PI * phase);
    for (k = 0; k < 3; k++) {
        hydrostatic[k] = interpolate(NIELL_HYDROSTATIC_AVG[k], abs_lat)
                         - interpolate(NIELL_HYDROSTATIC_AMP[k], abs_lat) * season;
        wet[k] = interpolate(NIELL_WET[k], abs_lat);
    }
}

/* [spec: TM-23 5.4.2.2.1 eqs. (5.67), (5.68)] M_dry = m(E, a_d, b_d, c_d) +
 * (1/sin E - m(E, a_ht, b_ht, c_ht)) H with H in km; M_wet = m(E, a_w, b_w,
 * c_w). H is the ellipsoidal height [spec: SPARTN v2.0.3 8.12, m_h(e, h)
 * "for ... ellipsoidal height (h)"], where TM-23 says height above sea
 * level. */
int32_t ssr_troposphere_niell_mapping(double latitude, double height, double doy, double elevation,
                                      double *hydrostatic, double *wet)
{
    double hc[3], wc[3];
    double sin_e, height_km;
    if (!(elevation > 0.0 && elevation <= 0.5 * SSR_PI)) {
        return SSR_TROPOSPHERE_ELEVATION;
    }
    ssr_troposphere_niell_coefficients(latitude, doy, hc, wc);
    sin_e = sin(elevation);
    height_km = height / 1000.0;
    *hydrostatic = ssr_troposphere_marini(sin_e, hc[0], hc[1], hc[2])
                   + (1.0 / sin_e - ssr_troposphere_marini(sin_e, NIELL_HEIGHT[0], NIELL_HEIGHT[1], NIELL_HEIGHT[2]))
                         * height_km;
    *wet = ssr_troposphere_marini(sin_e, wc[0], wc[1], wc[2]);
    return SSR_TROPOSPHERE_OK;
}

/* [spec: SPARTN v2.0.3 8.12] dT = T_h m_h(E, h) + (T_P + T_G) m_w(E). */
double ssr_troposphere_slant(double hydrostatic_zenith, double polynomial_residual,
                             double grid_residual, double hydrostatic_mapping,
                             double wet_mapping)
{
    return hydrostatic_zenith * hydrostatic_mapping
           + (polynomial_residual + grid_residual) * wet_mapping;
}

/* [spec: UNB3m eq. (3); TM-23 5.4.2.1 eq. (5.62) for D_min and |latitude|]
 * X = X_avg - X_amp cos(2 pi (doy - D_min) / 365.25). */
void ssr_troposphere_unb3m_temperature(double latitude, double doy, double *temperature, double *lapse_rate)
{
    double lat_deg = latitude * 180.0 / SSR_PI;
    double abs_lat = fabs(lat_deg);
    double dmin = lat_deg < 0.0 ? UNB3_DMIN_SOUTH : UNB3_DMIN_NORTH;
    double season = cos(2.0 * SSR_PI * (doy - dmin) / YEAR);
    *temperature = interpolate(UNB3_T_AVG, abs_lat) - interpolate(UNB3_T_AMP, abs_lat) * season;
    *lapse_rate = interpolate(UNB3_BETA_AVG, abs_lat) - interpolate(UNB3_BETA_AMP, abs_lat) * season;
}

/* [spec: UNB3m eq. (6)] g_m = 9.784 (1 - 2.66e-3 cos(2 phi) - 2.8e-7 H). */
static double centroid_gravity(double latitude, double height)
{
    return 9.784 * (1.0 - 2.66e-3 * cos(2.0 * latitude) - 2.8e-7 * height);
}

/* [spec: UNB3m eq. (4), d_zh = 1e-6 k1 R / g_m * P0 * (1 - beta H / T0)^(g /
 * (R beta)), divided by its value at H = 0, so that P0 and k1 cancel; eq. (6)
 * for g_m] [design: DESIGN Q7. SPARTN 8.12 refers the area's hydrostatic
 * delay to 0 m ellipsoidal height and gives no reduction to the user's
 * height; this ratio applies the UNB3m height dependence to it. Heights are
 * counted from the ellipsoid, the reference SPARTN names, where UNB3m uses
 * orthometric height with its table at sea level] */
int32_t ssr_troposphere_unb3m_height_factor(double latitude, double temperature,
                                            double lapse_rate, double height, double *factor)
{
    double base = 1.0 - lapse_rate * height / temperature;
    if (!(base > 0.0)) {
        return SSR_TROPOSPHERE_HEIGHT;
    }
    *factor = centroid_gravity(latitude, 0.0) / centroid_gravity(latitude, height)
              * pow(base, G_SURFACE / (R_DRY * lapse_rate));
    return SSR_TROPOSPHERE_OK;
}
