/* UNB3m zenith delays and the look-up rows they need. See ssrcore/unb3m.h. */

#include <math.h>

#include "ssrcore/geodesy.h"
#include "ssrcore/unb3m.h"

/* [spec: UNB3m eq. (3)] Length of the seasonal period, days. */
static const double YEAR = 365.25;

/* [spec: UNB3m eq. (3), "the origin of the yearly variation is day of year
 * (doy) 28"; TM-23 5.4.2.1 eq. (5.62), D_min = 211 for southern latitudes]
 * [design: UNB3m was validated in North America only and states no southern
 * case; TM-23 shifts the same table by half a year there, as the troposphere
 * kernel's ssr_troposphere_unb3m_temperature does, so that P0, RH and lambda
 * keep the phase of the T0 and beta they are combined with] */
static const double DMIN_NORTH = 28.0;
static const double DMIN_SOUTH = 211.0;

/* [spec: UNB3m eq. (4) and (5) text] Gas constant for dry air R,
 * J/(kg K). */
static const double R_DRY = 287.054;

/* [spec: TM-23 5.4.2.1 text after eq. (5.64), g = 9.80665 m/s^2] Surface
 * acceleration of gravity g of UNB3m eqs. (4), (5), which UNB3m names and
 * does not state. */
static const double G_SURFACE = 9.80665;

/* [spec: UNB3m eq. (5) text] Refractivity constants k1 (K/mbar), k'2
 * (K/mbar) and k3 (K^2/mbar). */
static const double K1 = 77.60;
static const double K2_PRIME = 16.6;
static const double K3 = 377600.0;

/* [spec: UNB3m Table 1] Pressure P (mbar) and water vapour lapse rate lambda:
 * averages, then amplitudes, at 15, 30, 45, 60 and 75 degrees. */
static const double PRESSURE_AVG[5] = {1013.25, 1017.25, 1015.75, 1011.75, 1013.00};
static const double PRESSURE_AMP[5] = {0.00, -3.75, -2.25, -1.75, -0.50};
static const double LAMBDA_AVG[5] = {2.77, 3.15, 2.57, 1.81, 1.55};
static const double LAMBDA_AMP[5] = {0.00, 0.33, 0.46, 0.74, 0.30};

/* [spec: UNB3m Table 2] Relative humidity (%), average and amplitude, at 15
 * to 75 degrees; it replaces the water vapour pressure columns of Table 1. */
static const double HUMIDITY_AVG[5] = {75.0, 80.0, 76.0, 77.5, 82.5};
static const double HUMIDITY_AMP[5] = {0.0, 0.0, -1.0, -2.5, 2.5};

/* [spec: UNB3m eqs. (1), (2)] Linear interpolation of a row tabled at 15,
 * 30, ..., 75 degrees, at |latitude| in degrees; the end values at or below
 * 15 and at or above 75 degrees. */
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

/* [spec: UNB3m eqs. (1) to (3), Tables 1 and 2; TM-23 5.4.2.1 eq. (5.62)]
 * X = Avg(phi) - Amp(phi) cos(2 pi (doy - D_min) / 365.25). */
void ssr_unb3m_parameters(double latitude, double doy, double *pressure, double *humidity, double *lambda)
{
    double lat_deg = latitude * 180.0 / SSR_PI;
    double abs_lat = fabs(lat_deg);
    double dmin = lat_deg < 0.0 ? DMIN_SOUTH : DMIN_NORTH;
    double season = cos(2.0 * SSR_PI * (doy - dmin) / YEAR);
    *pressure = interpolate(PRESSURE_AVG, abs_lat) - interpolate(PRESSURE_AMP, abs_lat) * season;
    *humidity = interpolate(HUMIDITY_AVG, abs_lat) - interpolate(HUMIDITY_AMP, abs_lat) * season;
    *lambda = interpolate(LAMBDA_AVG, abs_lat) - interpolate(LAMBDA_AMP, abs_lat) * season;
}

/* [spec: UNB3m eqs. (9) to (11)] e0 = RH/100 e_s f_w, with
 * e_s = 0.01 exp(1.2378847e-5 T0^2 - 1.9121316e-2 T0 + 33.93711047
 * - 6.3431645e3 / T0) and f_w = 1.00062 + 3.14e-6 P0 + 5.6e-7 (T0 - 273.15)^2. */
double ssr_unb3m_vapour_pressure(double humidity, double temperature, double pressure)
{
    double t_c = temperature - 273.15;
    double saturation = 0.01 * exp(1.2378847e-5 * temperature * temperature - 1.9121316e-2 * temperature
                                   + 33.93711047 - 6.3431645e3 / temperature);
    double enhancement = 1.00062 + 3.14e-6 * pressure + 5.6e-7 * t_c * t_c;
    return humidity / 100.0 * saturation * enhancement;
}

/* [spec: UNB3m eq. (6)] g_m = 9.784 (1 - 2.66e-3 cos(2 phi) - 2.8e-7 H). */
static double centroid_gravity(double latitude, double height)
{
    return 9.784 * (1.0 - 2.66e-3 * cos(2.0 * latitude) - 2.8e-7 * height);
}

/* [spec: UNB3m eq. (4), d_zh = 1e-6 k1 R / g_m P0 (1 - beta H / T0)^(g / (R
 * beta)); eq. (5), d_znh = 1e-6 (T_m k'2 + k3) R / (g_m lambda' - beta R)
 * e0 / T0 (1 - beta H / T0)^(lambda' g / (R beta) - 1); eq. (6) g_m; eq. (7)
 * T_m = T0 (1 - beta R / (g_m lambda')); lambda' = lambda + 1] H is the
 * height the caller passes; UNB3m defines it as orthometric. */
int32_t ssr_unb3m_zenith(double latitude, double height, double pressure, double temperature,
                         double vapour_pressure, double lapse_rate, double lambda,
                         double *hydrostatic, double *wet)
{
    double g_m, base, exponent, lambda_1, denominator, t_m;
    if (!(temperature > 0.0) || !(lapse_rate > 0.0)) {
        return SSR_UNB3M_PARAMETERS;
    }
    base = 1.0 - lapse_rate * height / temperature;
    if (!(base > 0.0)) {
        return SSR_UNB3M_HEIGHT;
    }
    g_m = centroid_gravity(latitude, height);
    lambda_1 = lambda + 1.0;
    denominator = g_m * lambda_1 - lapse_rate * R_DRY;
    if (!(denominator > 0.0)) {
        return SSR_UNB3M_PARAMETERS;
    }
    exponent = G_SURFACE / (R_DRY * lapse_rate);
    t_m = temperature * (1.0 - lapse_rate * R_DRY / (g_m * lambda_1));
    *hydrostatic = 1e-6 * K1 * R_DRY / g_m * pressure * pow(base, exponent);
    *wet = 1e-6 * (t_m * K2_PRIME + K3) * R_DRY / denominator * vapour_pressure / temperature
           * pow(base, lambda_1 * exponent - 1.0);
    return SSR_UNB3M_OK;
}
