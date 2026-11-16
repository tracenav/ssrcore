/* Solid Earth tide displacement of a site. See ssrcore/tides.h. */

#include <math.h>
#include <stdint.h>

#include "ssrcore/geodesy.h"
#include "ssrcore/tides.h"

/* [design: SSR_PI, the double nearest pi (DESIGN 4); TN36 is not a SPARTN or
 * ICD formula] */
#define PI SSR_PI
#define DEG (PI / 180.0)
#define ARCSEC (PI / 648000.0)
#define RE SSR_IERS_EARTH_RADIUS

/* [spec: TN36 7.1.1, the definitions after eqs. (1) and (5)] The site's
 * geocentric frame: r, the unit radial vector; n, perpendicular to r and
 * northward; e, east; with sin and cos of the geocentric latitude phi and
 * the longitude lambda. */
typedef struct {
    double r[3], n[3], e[3];
    double sin_phi, cos_phi, lambda;
} site_frame;

/* [spec: TN36 7.1.1, as above] the frame of a site given in ECEF. */
static site_frame frame_of(const double site[3])
{
    site_frame f;
    double p = sqrt(site[0] * site[0] + site[1] * site[1]);
    double rho = sqrt(p * p + site[2] * site[2]);
    double sin_l;
    double cos_l;
    f.sin_phi = site[2] / rho;
    f.cos_phi = p / rho;
    f.lambda = atan2(site[1], site[0]);
    sin_l = sin(f.lambda);
    cos_l = cos(f.lambda);
    f.r[0] = site[0] / rho;
    f.r[1] = site[1] / rho;
    f.r[2] = site[2] / rho;
    f.n[0] = -f.sin_phi * cos_l;
    f.n[1] = -f.sin_phi * sin_l;
    f.n[2] = f.cos_phi;
    f.e[0] = -sin_l;
    f.e[1] = cos_l;
    f.e[2] = 0.0;
    return f;
}

/* [design: the radial, north and east components of TN36's displacements
 * added to out in ECEF] */
static void add_local(const site_frame *f, double radial, double north, double east, double out[3])
{
    int32_t i;
    for (i = 0; i < 3; i++) {
        out[i] += radial * f->r[i] + north * f->n[i] + east * f->e[i];
    }
}

/* [spec: TN36 7.1.1 eq. (5) with eq. (2)]
 * dr = GMj Re^4 / (GM Rj^3) { h2 r (3 (R.r)^2 - 1)/2 + 3 l2 (R.r) [R - (R.r) r] }. */
void ssr_tides_degree2(const double site[3], const double body[3], double mass_ratio,
                      double h0, double h2, double l0, double l2, double out[3])
{
    site_frame f = frame_of(site);
    double distance = sqrt(body[0] * body[0] + body[1] * body[1] + body[2] * body[2]);
    double k = mass_ratio * (RE * RE * RE * RE) / (distance * distance * distance);
    double p2 = (3.0 * f.sin_phi * f.sin_phi - 1.0) / 2.0;
    double h = h0 + h2 * p2;
    double l = l0 + l2 * p2;
    double u[3];
    double s;
    int32_t i;
    for (i = 0; i < 3; i++) {
        u[i] = body[i] / distance;
    }
    s = u[0] * f.r[0] + u[1] * f.r[1] + u[2] * f.r[2];
    for (i = 0; i < 3; i++) {
        out[i] = k * (h * f.r[i] * (3.0 * s * s - 1.0) / 2.0 + 3.0 * l * s * (u[i] - s * f.r[i]));
    }
}

/* [spec: TN36 7.1.1 eq. (6)] dr = GMj Re^5 / (GM Rj^4) { h3 r (5/2 (R.r)^3
 * - 3/2 (R.r)) + l3 (15/2 (R.r)^2 - 3/2) [R - (R.r) r] }, h3 = 0.292,
 * l3 = 0.015. */
void ssr_tides_degree3(const double site[3], const double body[3], double mass_ratio,
                      double out[3])
{
    const double h3 = 0.292;
    const double l3 = 0.015;
    site_frame f = frame_of(site);
    double distance = sqrt(body[0] * body[0] + body[1] * body[1] + body[2] * body[2]);
    double k = mass_ratio * (RE * RE * RE * RE * RE) / (distance * distance * distance * distance);
    double u[3];
    double s;
    int32_t i;
    for (i = 0; i < 3; i++) {
        u[i] = body[i] / distance;
    }
    s = u[0] * f.r[0] + u[1] * f.r[1] + u[2] * f.r[2];
    for (i = 0; i < 3; i++) {
        out[i] = k * (h3 * f.r[i] * (2.5 * s * s * s - 1.5 * s)
                      + l3 * (7.5 * s * s - 1.5) * (u[i] - s * f.r[i]));
    }
}

/* [spec: TN36 7.1.1 "In the expressions given below", Phi_j and lambda_j]
 * the body's geocentric latitude and east longitude, and GMj Re^4 / (GM Rj^3). */
static void body_angles(const double body[3], double mass_ratio,
                        double *sin_lat, double *cos_lat, double *lon, double *k)
{
    double p = sqrt(body[0] * body[0] + body[1] * body[1]);
    double distance = sqrt(p * p + body[2] * body[2]);
    *sin_lat = body[2] / distance;
    *cos_lat = p / distance;
    *lon = atan2(body[1], body[0]);
    *k = mass_ratio * (RE * RE * RE * RE) / (distance * distance * distance);
}

/* [spec: TN36 7.1.1 eqs. (7), (8), (9)]
 * diurnal:     dt = -l1 sin phi K P21(sin Phi) [sin phi cos(lambda - lambda_j) n
 *                                              - cos 2phi sin(lambda - lambda_j) e],
 * semidiurnal: dt = -1/2 l1 sin phi cos phi K P22(sin Phi) [cos 2(lambda - lambda_j) n
 *                                              + sin phi sin 2(lambda - lambda_j) e],
 * with P21(sin Phi) = 3 sin Phi cos Phi and P22(sin Phi) = 3 cos^2 Phi. */
void ssr_tides_l1(const double site[3], const double body[3], double mass_ratio,
                 double l1_diurnal, double l1_semidiurnal, double out[3])
{
    site_frame f = frame_of(site);
    double sin_lat, cos_lat, lon, k;
    double p21, p22, dl, cos_2phi;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    body_angles(body, mass_ratio, &sin_lat, &cos_lat, &lon, &k);
    p21 = 3.0 * sin_lat * cos_lat;
    p22 = 3.0 * cos_lat * cos_lat;
    dl = f.lambda - lon;
    cos_2phi = f.cos_phi * f.cos_phi - f.sin_phi * f.sin_phi;
    add_local(&f, 0.0,
              -l1_diurnal * f.sin_phi * k * p21 * f.sin_phi * cos(dl),
              l1_diurnal * f.sin_phi * k * p21 * cos_2phi * sin(dl), out);
    add_local(&f, 0.0,
              -0.5 * l1_semidiurnal * f.sin_phi * f.cos_phi * k * p22 * cos(2.0 * dl),
              -0.5 * l1_semidiurnal * f.sin_phi * f.cos_phi * k * p22 * f.sin_phi * sin(2.0 * dl), out);
}

/* [spec: TN36 7.1.1 eqs. (10a), (10b), (11a), (11b)]
 * diurnal:     dr = -3/4 hI K sin 2Phi sin 2phi sin(lambda - lambda_j),
 *              dt = -3/2 lI K sin 2Phi [cos 2phi sin(lambda - lambda_j) n
 *                                       + sin phi cos(lambda - lambda_j) e];
 * semidiurnal: dr = -3/4 hI K cos^2 Phi cos^2 phi sin 2(lambda - lambda_j),
 *              dt = 3/4 lI K cos^2 Phi [sin 2phi sin 2(lambda - lambda_j) n
 *                                       - 2 cos phi cos 2(lambda - lambda_j) e]. */
void ssr_tides_out_of_phase(const double site[3], const double body[3], double mass_ratio,
                           double hi_diurnal, double li_diurnal,
                           double hi_semidiurnal, double li_semidiurnal, double out[3])
{
    site_frame f = frame_of(site);
    double sin_lat, cos_lat, lon, k;
    double dl, sin_2lat, cos2_lat, sin_2phi, cos_2phi;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    body_angles(body, mass_ratio, &sin_lat, &cos_lat, &lon, &k);
    dl = f.lambda - lon;
    sin_2lat = 2.0 * sin_lat * cos_lat;
    cos2_lat = cos_lat * cos_lat;
    sin_2phi = 2.0 * f.sin_phi * f.cos_phi;
    cos_2phi = f.cos_phi * f.cos_phi - f.sin_phi * f.sin_phi;
    add_local(&f,
              -0.75 * hi_diurnal * k * sin_2lat * sin_2phi * sin(dl),
              -1.5 * li_diurnal * k * sin_2lat * cos_2phi * sin(dl),
              -1.5 * li_diurnal * k * sin_2lat * f.sin_phi * cos(dl), out);
    add_local(&f,
              -0.75 * hi_semidiurnal * k * cos2_lat * f.cos_phi * f.cos_phi * sin(2.0 * dl),
              0.75 * li_semidiurnal * k * cos2_lat * sin_2phi * sin(2.0 * dl),
              -0.75 * li_semidiurnal * k * cos2_lat * 2.0 * f.cos_phi * cos(2.0 * dl), out);
}

typedef struct {
    int8_t n[5];          /* multipliers of l, l', F, D, Omega */
    double dr_ip, dr_op, dt_ip, dt_op;   /* mm */
} step2_term;

/* [spec: TN36 7.1.1 Table 7.3a] Diurnal band: Delaunay multipliers and
 * Delta R_f (ip, op), Delta T_f (ip, op) in mm. */
static const step2_term DIURNAL[11] = {
    {{1, 0, 2, 0, 2}, -0.08, 0.00, -0.01, 0.01},     /* Q1 */
    {{0, 0, 2, 0, 1}, -0.10, 0.00, 0.00, 0.00},      /* 145,545 */
    {{0, 0, 2, 0, 2}, -0.51, 0.00, -0.02, 0.03},     /* O1 */
    {{1, 0, 0, 0, 0}, 0.06, 0.00, 0.00, 0.00},       /* NO1 */
    {{0, 1, 2, -2, 2}, -0.06, 0.00, 0.00, 0.00},     /* pi1 */
    {{0, 0, 2, -2, 2}, -1.23, -0.07, 0.06, 0.01},    /* P1 */
    {{0, 0, 0, 0, -1}, -0.22, 0.01, 0.01, 0.00},     /* 165,545 */
    {{0, 0, 0, 0, 0}, 12.00, -0.78, -0.67, -0.03},   /* K1 */
    {{0, 0, 0, 0, 1}, 1.73, -0.12, -0.10, 0.00},     /* 165,565 */
    {{0, -1, 0, 0, 0}, -0.50, -0.01, 0.03, 0.00},    /* psi1 */
    {{0, 0, -2, 2, -2}, -0.11, 0.01, 0.01, 0.00},    /* phi1 */
};

/* [spec: TN36 7.1.1 Table 7.3b] Long-period band, as Table 7.3a. */
static const step2_term LONG_PERIOD[5] = {
    {{0, 0, 0, 0, 1}, 0.47, 0.16, 0.23, 0.07},       /* 55,565 */
    {{0, 0, -2, 2, -2}, -0.20, -0.11, -0.12, -0.05}, /* Ssa */
    {{-1, 0, 0, 0, 0}, -0.11, -0.09, -0.08, -0.04},  /* Mm */
    {{0, 0, -2, 0, -2}, -0.13, -0.15, -0.11, -0.07}, /* Mf */
    {{0, 0, -2, 0, -1}, -0.05, -0.06, -0.05, -0.03}, /* 75,565 */
};

/* [spec: TN36 ch. 5 eq. (5.43)] The Delaunay arguments l, l', F, D, Omega
 * (rad), t in TT centuries. */
static void delaunay(double t, double f[5])
{
    f[0] = 134.96340251 * DEG + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * -0.00024470))) * ARCSEC;
    f[1] = 357.52910918 * DEG + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * -0.00001149))) * ARCSEC;
    f[2] = 93.27209062 * DEG + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * 0.00000417))) * ARCSEC;
    f[3] = 297.85019547 * DEG + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * -0.00003169))) * ARCSEC;
    f[4] = 125.04455501 * DEG + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * -0.00005939))) * ARCSEC;
}

/* [spec: TN36 ch. 6 section 6.2, theta_f = m (theta_g + pi) - sum N_j F_j] */
static double tide_argument(const step2_term *term, int32_t m, double gmst, const double f[5])
{
    double theta = m * (gmst + PI);
    int32_t j;
    for (j = 0; j < 5; j++) {
        theta -= term->n[j] * f[j];
    }
    return theta;
}

/* [spec: TN36 7.1.1 eqs. (12a), (12b) for Table 7.3a; (13a), (13b) for
 * Table 7.3b]
 * diurnal:      dr = [dR_ip sin(theta + lambda) + dR_op cos(theta + lambda)] sin 2phi,
 *               dt = [dT_ip cos(theta + lambda) - dT_op sin(theta + lambda)] sin phi e
 *                  + [dT_ip sin(theta + lambda) + dT_op cos(theta + lambda)] cos 2phi n;
 * long period:  dr = (3/2 sin^2 phi - 1/2)(dR_ip cos theta + dR_op sin theta),
 *               dt = (dT_ip cos theta + dT_op sin theta) sin 2phi n. */
void ssr_tides_step2(const double site[3], double t, double gmst, double out[3])
{
    site_frame s = frame_of(site);
    double f[5];
    double sin_2phi = 2.0 * s.sin_phi * s.cos_phi;
    double cos_2phi = s.cos_phi * s.cos_phi - s.sin_phi * s.sin_phi;
    double radial = 0.0;
    double north = 0.0;
    double east = 0.0;
    int32_t i;
    delaunay(t, f);
    for (i = 0; i < 11; i++) {
        const step2_term *d = &DIURNAL[i];
        double a = tide_argument(d, 1, gmst, f) + s.lambda;
        radial += (d->dr_ip * sin(a) + d->dr_op * cos(a)) * sin_2phi;
        east += (d->dt_ip * cos(a) - d->dt_op * sin(a)) * s.sin_phi;
        north += (d->dt_ip * sin(a) + d->dt_op * cos(a)) * cos_2phi;
    }
    for (i = 0; i < 5; i++) {
        const step2_term *z = &LONG_PERIOD[i];
        double a = tide_argument(z, 0, gmst, f);
        radial += (1.5 * s.sin_phi * s.sin_phi - 0.5) * (z->dr_ip * cos(a) + z->dr_op * sin(a));
        north += (z->dt_ip * cos(a) + z->dt_op * sin(a)) * sin_2phi;
    }
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    add_local(&s, radial * 1e-3, north * 1e-3, east * 1e-3, out);
}

/* [spec: TN36 7.1.1, steps 1 and 2, nominal values of the summary table:
 * h(0) = 0.6078, h(2) = -0.0006, l(0) = 0.0847, l(2) = 0.0002; l(1) 0.0012
 * diurnal, 0.0024 semidiurnal; hI, lI -0.0025, -0.0007 diurnal and -0.0022,
 * -0.0007 semidiurnal; degree 3 for the Moon only] */
void ssr_tides_displacement(const double site[3], const double sun[3], const double moon[3],
                            double t, double gmst, double out[3])
{
    const double *bodies[2];
    double ratios[2];
    double part[3];
    int32_t b;
    int32_t i;
    bodies[0] = moon;
    bodies[1] = sun;
    ratios[0] = SSR_IERS_MOON_EARTH_MASS_RATIO;
    ratios[1] = SSR_IERS_SUN_EARTH_MASS_RATIO;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    for (b = 0; b < 2; b++) {
        ssr_tides_degree2(site, bodies[b], ratios[b], 0.6078, -0.0006, 0.0847, 0.0002, part);
        for (i = 0; i < 3; i++) {
            out[i] += part[i];
        }
        ssr_tides_l1(site, bodies[b], ratios[b], 0.0012, 0.0024, part);
        for (i = 0; i < 3; i++) {
            out[i] += part[i];
        }
        ssr_tides_out_of_phase(site, bodies[b], ratios[b], -0.0025, -0.0007, -0.0022, -0.0007, part);
        for (i = 0; i < 3; i++) {
            out[i] += part[i];
        }
    }
    ssr_tides_degree3(site, moon, SSR_IERS_MOON_EARTH_MASS_RATIO, part);
    for (i = 0; i < 3; i++) {
        out[i] += part[i];
    }
    ssr_tides_step2(site, t, gmst, part);
    for (i = 0; i < 3; i++) {
        out[i] += part[i];
    }
}

/* [spec: TN36 7.1.1.2 eqs. (14a), (14b)] radial [-0.1206 + 0.0001 P2] P2 and
 * northward [-0.0252 - 0.0001 P2] sin 2phi, metres, P2 = (3 sin^2 phi - 1)/2. */
void ssr_tides_permanent(const double site[3], double out[3])
{
    site_frame f = frame_of(site);
    double p2 = (3.0 * f.sin_phi * f.sin_phi - 1.0) / 2.0;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    add_local(&f, (-0.1206 + 0.0001 * p2) * p2,
              (-0.0252 - 0.0001 * p2) * 2.0 * f.sin_phi * f.cos_phi, 0.0, out);
}
