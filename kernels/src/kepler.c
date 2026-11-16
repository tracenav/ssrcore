/* Keplerian broadcast orbit and clock. See ssrcore/kepler.h. */

#include <math.h>
#include <stdint.h>

#include "ssrcore/kepler.h"

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2; Galileo OS SIS
 * ICD 2.2 Table 68; BDS-SIS-ICD-B1I-3.0 Table 5-11] pi as the ICDs fix it for
 * orbit computation. */
static const double PI_ICD = 3.1415926535898;

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV] WGS
 * 84 values for the GPS user; also used for QZSS [design: IS-QZSS-PNT is not
 * staged (DESIGN Q14), and QZSS broadcasts the GPS LNAV format] */
static const double MU_GPS = 3.986005e14;
static const double OMEGA_E_GPS = 7.2921151467e-5;
/* [spec: Galileo OS SIS ICD 2.2 5.1.1 Table 68] */
static const double MU_GALILEO = 3.986004418e14;
static const double OMEGA_E_GALILEO = 7.2921151467e-5;
/* [spec: BDS-SIS-ICD-B1I-3.0 5.2.4.12 Table 5-11] BDCS values. */
static const double MU_BEIDOU = 3.986004418e14;
static const double OMEGA_E_BEIDOU = 7.2921150e-5;

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV]
 * seconds in a week, for the actual total time difference across weeks. */
static const double WEEK_SECONDS = 604800.0;

/* [design: Table 20-IV asks for a minimum of three iterations and gives no
 * stopping rule. Newton's iteration on Kepler's equation converges
 * quadratically for GNSS eccentricities (below 0.03 for every broadcast
 * orbit, and the field ranges allow at most 0.5 for BeiDou, Table 5-10);
 * 1e-14 rad is 0.3 um along a GNSS orbit, and 30 iterations is far beyond
 * what e < 0.5 needs, so a failure means bad elements, not slow convergence.] */
#define KEPLER_MIN_ITERATIONS 3
#define KEPLER_MAX_ITERATIONS 30
static const double KEPLER_TOLERANCE = 1e-14;

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV note
 * *; Galileo OS SIS ICD 2.2 Table 68 note *; BDS-SIS-ICD-B1I-3.0 Table 5-11
 * note *] The actual total time difference t - ref. With full weeks on both
 * sides the difference is exact, so the half-week crossover rule the notes give
 * for a time of week alone is not needed. The difference is taken in seconds
 * before any multiplication (DESIGN 4). */
static double since(int32_t week, double sow, int32_t ref_week, double ref_sow)
{
    return (double)(week - ref_week) * WEEK_SECONDS + (sow - ref_sow);
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2] Radians from semicircles. */
double ssr_kepler_semicircles(double semicircles)
{
    return semicircles * PI_ICD;
}

/* [spec: BDS-SIS-ICD-B1I-3.0 4.3 Table 4-1] SV IDs 1-5 and 59-63 are GEO. */
int32_t ssr_kepler_beidou_geo(int32_t prn)
{
    return (prn >= 1 && prn <= 5) || (prn >= 59 && prn <= 63);
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV] E0 =
 * M; Ej = Ej-1 + (M - Ej-1 + e sin Ej-1) / (1 - e cos Ej-1), at least three
 * times; stopping rule [design: see KEPLER_TOLERANCE]. */
int32_t ssr_kepler_eccentric_anomaly(double m, double e, double *ecc_anomaly)
{
    double ek = m;
    int32_t j;
    for (j = 1; j <= KEPLER_MAX_ITERATIONS; j++) {
        double step = (m - ek + e * sin(ek)) / (1.0 - e * cos(ek));
        ek += step;
        if (j >= KEPLER_MIN_ITERATIONS && fabs(step) <= KEPLER_TOLERANCE) {
            *ecc_anomaly = ek;
            return SSR_KEPLER_OK;
        }
    }
    *ecc_anomaly = ek;
    return SSR_KEPLER_NOT_CONVERGED;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV
 * sheets 1 to 3 (position and SV velocity); Galileo OS SIS ICD 2.2 5.1.1 Table
 * 68; BDS-SIS-ICD-B1I-3.0 5.2.4.12 Table 5-11, MEO/IGSO and GEO]
 *
 * The true anomaly is atan2(sqrt(1 - e^2) sin E, cos E - e), the form of
 * the Galileo and BeiDou tables; the 2 atan(sqrt((1+e)/(1-e)) tan(E/2)) of
 * IS-GPS-200N with IRN-IS-200N-003 and -004 is the same angle, and atan2 has
 * no singularity at E = pi.
 *
 * Velocity: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV sheet 3. The
 * Galileo and BeiDou ICDs give no velocity [design: their position equations
 * are those of Table 20-IV with other constants, so the time derivative of the
 * position is the same formulae; for a BeiDou GEO the derivative of the table's
 * position is the same formulae in the inertial plane, with the node rate
 * omega_dot alone, then rotated, plus the derivative of the rotation
 * Rz(omega_e tk), which is omega_e (y, -x, 0) of the rotated position]. */
int32_t ssr_kepler_state(const ssr_kepler_orbit *orbit, int32_t week, double sow,
                         double position[3], double velocity[3], double *ecc_anomaly)
{
    double mu, omega_e;
    switch (orbit->system) {
    case SSR_SYSTEM_GPS:
    case SSR_SYSTEM_QZSS:    mu = MU_GPS;     omega_e = OMEGA_E_GPS;     break;
    case SSR_SYSTEM_GALILEO: mu = MU_GALILEO; omega_e = OMEGA_E_GALILEO; break;
    case SSR_SYSTEM_BEIDOU:  mu = MU_BEIDOU;  omega_e = OMEGA_E_BEIDOU;  break;
    default: return SSR_KEPLER_BAD_ELEMENTS;
    }
    const double e = orbit->e;
    if (!(e >= 0.0 && e < 1.0) || !(orbit->sqrt_a > 0.0))
        return SSR_KEPLER_BAD_ELEMENTS;
    const int32_t geo = orbit->system == SSR_SYSTEM_BEIDOU && orbit->geo;

    const double a = orbit->sqrt_a * orbit->sqrt_a;
    const double n0 = sqrt(mu / (a * a * a));
    const double tk = since(week, sow, orbit->toe_week, orbit->toe);
    const double n = n0 + orbit->delta_n;
    const double mk = orbit->m0 + n * tk;

    double ek;
    int32_t status = ssr_kepler_eccentric_anomaly(mk, e, &ek);
    if (status != SSR_KEPLER_OK)
        return status;
    *ecc_anomaly = ek;

    const double sin_e = sin(ek), cos_e = cos(ek);
    const double one_minus_ecos = 1.0 - e * cos_e;
    const double root = sqrt(1.0 - e * e);
    const double nu = atan2(root * sin_e, cos_e - e);
    const double phi = nu + orbit->omega;
    const double sin2 = sin(2.0 * phi), cos2 = cos(2.0 * phi);

    const double du = orbit->cus * sin2 + orbit->cuc * cos2;
    const double dr = orbit->crs * sin2 + orbit->crc * cos2;
    const double di = orbit->cis * sin2 + orbit->cic * cos2;
    const double u = phi + du;
    const double r = a * one_minus_ecos + dr;
    const double i = orbit->i0 + di + orbit->idot * tk;
    const double sin_u = sin(u), cos_u = cos(u);
    const double xp = r * cos_u;
    const double yp = r * sin_u;

    /* Node: Earth-fixed for GPS, Galileo and BeiDou MEO/IGSO; inertial for
     * a BeiDou GEO (Table 5-11). */
    const double node_rate = geo ? orbit->omega_dot : orbit->omega_dot - omega_e;
    const double node = orbit->omega0 + node_rate * tk - omega_e * orbit->toe;
    const double sin_o = sin(node), cos_o = cos(node);
    const double sin_i = sin(i), cos_i = cos(i);

    double x = xp * cos_o - yp * cos_i * sin_o;
    double y = xp * sin_o + yp * cos_i * cos_o;
    double z = yp * sin_i;

    /* Table 20-IV sheet 3. */
    const double e_rate = n / one_minus_ecos;
    const double nu_rate = e_rate * root / one_minus_ecos;
    const double i_rate = orbit->idot + 2.0 * nu_rate * (orbit->cis * cos2 - orbit->cic * sin2);
    const double u_rate = nu_rate + 2.0 * nu_rate * (orbit->cus * cos2 - orbit->cuc * sin2);
    const double r_rate = e * a * e_rate * sin_e
                          + 2.0 * nu_rate * (orbit->crs * cos2 - orbit->crc * sin2);
    const double xp_rate = r_rate * cos_u - r * u_rate * sin_u;
    const double yp_rate = r_rate * sin_u + r * u_rate * cos_u;

    double vx = -xp * node_rate * sin_o + xp_rate * cos_o - yp_rate * sin_o * cos_i
                - yp * (node_rate * cos_o * cos_i - i_rate * sin_o * sin_i);
    double vy = xp * node_rate * cos_o + xp_rate * sin_o + yp_rate * cos_o * cos_i
                - yp * (node_rate * sin_o * cos_i + i_rate * cos_o * sin_i);
    double vz = yp_rate * sin_i + yp * i_rate * cos_i;

    if (geo) {
        /* [spec: BDS-SIS-ICD-B1I-3.0 Table 5-11] X = Rz(omega_e tk) Rx(-5 deg) X_GK,
         * Rx(p) = [1 0 0; 0 cos p sin p; 0 -sin p cos p],
         * Rz(p) = [cos p sin p 0; -sin p cos p 0; 0 0 1]. */
        const double tilt = -5.0 * PI_ICD / 180.0;
        const double sin_t = sin(tilt), cos_t = cos(tilt);
        const double spin = omega_e * tk;
        const double sin_s = sin(spin), cos_s = cos(spin);

        const double x1 = x;
        const double y1 = cos_t * y + sin_t * z;
        const double z1 = -sin_t * y + cos_t * z;
        const double vx1 = vx;
        const double vy1 = cos_t * vy + sin_t * vz;
        const double vz1 = -sin_t * vy + cos_t * vz;

        x = cos_s * x1 + sin_s * y1;
        y = -sin_s * x1 + cos_s * y1;
        z = z1;
        vx = cos_s * vx1 + sin_s * vy1 + omega_e * y;
        vy = -sin_s * vx1 + cos_s * vy1 - omega_e * x;
        vz = vz1;
    }

    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    return SSR_KEPLER_OK;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1 equation (2)
 * without delta t_r; Galileo OS SIS ICD 2.2 5.1.4 Eq. 14; BDS-SIS-ICD-B1I-3.0
 * 5.2.4.9] Evaluated as af0 + dt (af1 + dt af2). */
double ssr_kepler_clock(double af0, double af1, double af2,
                        int32_t toc_week, double toc, int32_t week, double sow)
{
    const double dt = since(week, sow, toc_week, toc);
    return af0 + dt * (af1 + dt * af2);
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1, the first
 * derivative, without the relativistic part] af1 + 2 af2 (t - toc). */
double ssr_kepler_clock_drift(double af1, double af2,
                              int32_t toc_week, double toc, int32_t week, double sow)
{
    const double dt = since(week, sow, toc_week, toc);
    return af1 + 2.0 * af2 * dt;
}
