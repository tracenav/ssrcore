/* One query's kernel calls, end to end (query.h). */
#include <math.h>
#include <stdint.h>

#include "query.h"
#include "ssrcore/area.h"
#include "ssrcore/astro.h"
#include "ssrcore/correction.h"
#include "ssrcore/geodesy.h"
#include "ssrcore/ionosphere.h"
#include "ssrcore/kepler.h"
#include "ssrcore/relativity.h"
#include "ssrcore/scale.h"
#include "ssrcore/tides.h"
#include "ssrcore/timescale.h"
#include "ssrcore/troposphere.h"
#include "ssrcore/windup.h"

/* [design: model/selection.py, FIRST_GUESS, TOLERANCE, MAX_ITERATIONS] */
#define SSRQ_FIRST_GUESS 0.075
#define SSRQ_TOLERANCE 1e-12
#define SSRQ_MAX_ITERATIONS 5
/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] one week in whole seconds */
#define SSRQ_WEEK 604800
/* SPARTN Table 6.2: URE SF024, troposphere quality SF042, ionosphere SF055 */
#define SSRQ_SF_URE 24u
#define SSRQ_SF_TROPO_QUALITY 42u
#define SSRQ_SF_IONO_QUALITY 55u

/* [design: model.types._split and GnssTime.normalised] t + dt: the whole
 * seconds of dt added exactly, its fraction to the fraction, the week
 * carried by floor division. */
static ssrq_time time_add(ssrq_time t, double dt)
{
    double whole_dt = floor(dt);
    double frac_dt = dt - whole_dt;
    int64_t seconds = (int64_t)t.second + (int64_t)whole_dt;
    double f = 0.0 + (t.fraction + frac_dt);
    double extra = floor(f);
    int64_t total = seconds + (int64_t)extra;
    int64_t weeks = total / SSRQ_WEEK;
    ssrq_time out;
    if (total % SSRQ_WEEK < 0) {
        weeks -= 1;
    }
    out.week = t.week + (int32_t)weeks;
    out.second = (int32_t)(total - weeks * SSRQ_WEEK);
    out.fraction = (f - extra) + 0.0;
    return out;
}

/* [design: model.types.GnssTime.__sub__] a - b, s: weeks and whole seconds
 * differenced as integers, the fractions apart. */
static double time_diff(ssrq_time a, ssrq_time b)
{
    int64_t whole = (int64_t)(a.week - b.week) * SSRQ_WEEK + (int64_t)(a.second - b.second);
    return (double)whole + (a.fraction - b.fraction);
}

/* [design: model.types.GnssTime.seconds] the second of week as one double */
static double time_seconds(ssrq_time t)
{
    return (double)t.second + t.fraction;
}

/* [spec: TM-23 B.1.2, B.2.1; design: engine._receiver_of, CPython's
 * math.degrees, x * (180 / pi)] */
int32_t ssrq_receiver(const ssrq_inputs *in, ssrq_receiver_state *r)
{
    static const double east[3] = {1.0, 0.0, 0.0};
    static const double north[3] = {0.0, 1.0, 0.0};
    static const double up[3] = {0.0, 0.0, 1.0};
    if (ssr_geodesy_ecef_to_geodetic(in->receiver, SSR_WGS84_A, SSR_WGS84_F, r->geodetic) != 0) {
        return SSRQ_GEODETIC;
    }
    ssr_geodesy_enu_to_ecef(r->geodetic[0], r->geodetic[1], east, r->east);
    ssr_geodesy_enu_to_ecef(r->geodetic[0], r->geodetic[1], north, r->north);
    ssr_geodesy_enu_to_ecef(r->geodetic[0], r->geodetic[1], up, r->up);
    r->latitude_deg = r->geodetic[0] * (180.0 / SSR_PI);
    r->longitude_deg = r->geodetic[1] * (180.0 / SSR_PI);
    return SSRQ_OK;
}

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2 through
 * ssr_kepler_semicircles; design: selection.prepare, once per choice] */
int32_t ssrq_choice(const ssrq_inputs *in, ssrq_choice_state *c)
{
    c->m0 = ssr_kepler_semicircles(in->m0);
    c->delta_n = ssr_kepler_semicircles(in->delta_n);
    c->omega = ssr_kepler_semicircles(in->omega);
    c->omega0 = ssr_kepler_semicircles(in->omega0);
    c->omega_dot = ssr_kepler_semicircles(in->omega_dot);
    c->i0 = ssr_kepler_semicircles(in->i0);
    c->idot = ssr_kepler_semicircles(in->idot);
    return SSRQ_OK;
}

/* [spec: TM-23 eq. (5.70); TN36 7.1.1 with the Sun and Moon of astro; SPARTN
 * 8.10 to 8.12 and rule A4; UNB3m for DESIGN Q7; design: engine._epoch_of,
 * _tide, the SPARTN provider's per-epoch values] */
int32_t ssrq_epoch(const ssrq_inputs *in, const ssrq_receiver_state *r, ssrq_epoch_state *e)
{
    int32_t day = in->t_rx.second / 86400;
    int32_t valid[4] = {1, 1, 1, 1};
    double residuals[4];
    double offsets[2];
    double t_rx = time_seconds(in->t_rx);
    int32_t k;
    if (ssr_timescale_day_of_year(in->t_rx.week * 7 + day, &e->day_of_year) != 0) {
        return SSRQ_MODELLED;
    }
    e->doy = (double)e->day_of_year + (t_rx - (double)day * 86400.0) / 86400.0;

    /* the Sun, the Moon and the tide: astro.sky, tides.solid_earth_tide */
    ssr_astro_sun_moon_ecef(in->t_rx.week, t_rx, -in->leap_seconds, e->sun, e->moon);
    e->tt_centuries = ssr_astro_tt_centuries(in->t_rx.week, t_rx);
    e->gmst = ssr_astro_gmst(in->t_rx.week, t_rx, -in->leap_seconds);
    ssr_tides_displacement(in->receiver, e->sun, e->moon, e->tt_centuries, e->gmst, e->tide_displacement);

    /* the area choice (rule A4): inside, then the nearest centroid */
    if (ssr_area_contains(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                          in->lon_nodes, in->lat_spacing, in->lon_spacing, &e->inside) != 0 || !e->inside) {
        return SSRQ_ATMOSPHERE;
    }
    e->cos_lat = cos(r->latitude_deg * (SSR_PI / 180.0));
    if (ssr_area_offsets(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                         in->lon_nodes, in->lat_spacing, in->lon_spacing, offsets) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    e->area_distance = offsets[0] * offsets[0] + pow(offsets[1] * e->cos_lat, 2.0);

    /* the zenith troposphere of the area: SPARTN 8.12, the grid of 8.13 */
    if (ssr_area_offsets(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                         in->lon_nodes, in->lat_spacing, in->lon_spacing, offsets) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    e->tropo_dlat = offsets[0];
    e->tropo_dlon = offsets[1];
    if (ssr_area_polynomial(in->tropo_equation, in->tropo_c, e->tropo_dlat, e->tropo_dlon,
                            &e->tropo_polynomial) != 0
        || ssr_area_grid(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                         in->lon_nodes, in->lat_spacing, in->lon_spacing, e->tropo_nodes,
                         e->tropo_weights) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    for (k = 0; k < 4; k++) {
        residuals[k] = in->tropo_residuals[e->tropo_nodes[k]];
    }
    if (ssr_area_interpolate(e->tropo_weights, residuals, valid, &e->tropo_grid) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    ssr_troposphere_unb3m_temperature(r->geodetic[0], e->doy, &e->temperature, &e->lapse_rate);
    if (ssr_troposphere_unb3m_height_factor(r->geodetic[0], e->temperature, e->lapse_rate, r->geodetic[2],
                                            &e->height_factor) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    e->hydrostatic_zenith = in->hydrostatic * e->height_factor;
    if (ssr_scale_field(SSRQ_SF_TROPO_QUALITY, in->tropo_quality, &e->tropo_sigma) != 0) {
        return SSRQ_SCALE;
    }
    e->wet_zenith = e->tropo_polynomial + e->tropo_grid;
    /* the satellite's ionosphere quality, kept with its block */
    if (ssr_scale_field(SSRQ_SF_IONO_QUALITY, in->iono_quality, &e->iono_sigma) != 0) {
        return SSRQ_SCALE;
    }
    return SSRQ_OK;
}

/* [spec: DESIGN 3.4 steps 4 to 6: the light time of Ashby 2003 section 2
 * with SPARTN 8.8; SPARTN 8.7; TM-23 B.3; Ashby 2003 Eq. (96); SPARTN 8.7
 * dt_rel; TN36 7.1.1 projected; SPARTN 8.14 and TM-23 5.5; SPARTN 8.11,
 * 8.12 with Niell; design: engine.query, selection.transmit] */
int32_t ssrq_query(const ssrq_inputs *in, const ssrq_receiver_state *r, const ssrq_choice_state *c,
                   const ssrq_epoch_state *e, ssrq_query_state *q)
{
    ssr_kepler_orbit orbit;
    int32_t valid[4] = {1, 1, 1, 1};
    double residuals[4];
    double offsets[2];
    double tau = SSRQ_FIRST_GUESS;
    double light, sow = 0.0;
    int32_t step, k, converged = 0;

    orbit.system = in->system;
    orbit.geo = 0;
    orbit.toe_week = in->toe_week;
    orbit.toe = in->toe;
    orbit.sqrt_a = in->sqrt_a;
    orbit.e = in->e;
    orbit.m0 = c->m0;
    orbit.delta_n = c->delta_n;
    orbit.omega = c->omega;
    orbit.omega0 = c->omega0;
    orbit.omega_dot = c->omega_dot;
    orbit.i0 = c->i0;
    orbit.idot = c->idot;
    orbit.cuc = in->cuc;
    orbit.cus = in->cus;
    orbit.crc = in->crc;
    orbit.crs = in->crs;
    orbit.cic = in->cic;
    orbit.cis = in->cis;

    /* transmit time: the light-time iteration */
    for (step = 1; step <= SSRQ_MAX_ITERATIONS; step++) {
        q->t_tx = time_add(in->t_rx, -tau);
        sow = time_seconds(q->t_tx);
        if (ssr_kepler_state(&orbit, q->t_tx.week, sow, q->broadcast_position, q->velocity,
                             &q->ecc_anomaly) != 0) {
            return SSRQ_ORBIT;
        }
        q->broadcast_clock = ssr_kepler_clock(in->af0, in->af1, in->af2, in->toc_week, in->toc,
                                              q->t_tx.week, sow);
        if (ssr_correction_position(q->broadcast_position, q->velocity, in->d_o, q->position) != 0
            || ssr_geodesy_range(q->position, q->velocity, in->receiver, &q->range, &q->range_rate) != 0
            || ssr_relativity_sagnac(q->position, in->receiver, &q->sagnac) != 0) {
            return SSRQ_LIGHT_TIME;
        }
        light = (q->range + q->sagnac) / SSR_SPEED_OF_LIGHT;
        if (fabs(light - tau) < SSRQ_TOLERANCE) {
            converged = 1;
            break;
        }
        tau = light;
    }
    if (!converged) {
        return SSRQ_LIGHT_TIME;
    }
    q->iterations = step;
    if (ssr_correction_orbit(q->broadcast_position, q->velocity, in->d_o, q->orbit_ecef) != 0) {
        return SSRQ_LIGHT_TIME;
    }

    /* the satellite's clock: SPARTN 8.7 */
    q->clock_drift = ssr_kepler_clock_drift(in->af1, in->af2, in->toc_week, in->toc, q->t_tx.week, sow);
    q->clock = ssr_correction_clock(q->broadcast_clock, in->d_c);
    if (ssr_scale_field(SSRQ_SF_URE, in->ure_index, &q->ure) != 0) {
        return SSRQ_SCALE;
    }

    /* geometry in the frame of reception */
    q->tau = time_diff(in->t_rx, q->t_tx);
    ssr_geodesy_rotate_z(q->position, SSR_WGS84_OMEGA_E * q->tau, q->arrival);
    if (ssr_geodesy_line_of_sight(q->arrival, in->receiver, r->geodetic[0], r->geodetic[1], q->unit,
                                  &q->elevation, &q->azimuth) != 0) {
        return SSRQ_GEOMETRY;
    }

    /* modelled terms */
    if (ssr_relativity_shapiro(q->position, in->receiver, &q->shapiro) != 0) {
        return SSRQ_MODELLED;
    }
    q->relativity_clock = ssr_relativity_clock(q->broadcast_position, q->velocity);
    q->solid_tide = -(q->unit[0] * e->tide_displacement[0] + q->unit[1] * e->tide_displacement[1]
                      + q->unit[2] * e->tide_displacement[2]);
    q->yaw = in->yaw_deg * (SSR_PI / 180.0);
    if (ssr_windup_frame_spartn(q->position, q->velocity, q->yaw, q->frame_x, q->frame_y, q->frame_z) != 0
        || ssr_windup_cycles(q->position, in->receiver, r->east, r->north, q->frame_x, q->frame_y,
                             &q->windup) != 0) {
        return SSRQ_MODELLED;
    }

    /* troposphere: SPARTN 8.12 with Niell's mapping */
    if (ssr_troposphere_niell_mapping(r->geodetic[0], r->geodetic[2], e->doy, q->elevation,
                                      &q->hydrostatic_mapping, &q->wet_mapping) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    q->tropo_slant = ssr_troposphere_slant(e->hydrostatic_zenith, e->tropo_polynomial, e->tropo_grid,
                                           q->hydrostatic_mapping, q->wet_mapping);

    /* ionosphere: SPARTN 8.11 */
    if (ssr_area_offsets(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                         in->lon_nodes, in->lat_spacing, in->lon_spacing, offsets) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    q->iono_dlat = offsets[0];
    q->iono_dlon = offsets[1];
    if (ssr_area_polynomial(in->iono_equation, in->iono_c, q->iono_dlat, q->iono_dlon,
                            &q->iono_polynomial) != 0
        || ssr_area_grid(r->latitude_deg, r->longitude_deg, in->ref_lat, in->ref_lon, in->lat_nodes,
                         in->lon_nodes, in->lat_spacing, in->lon_spacing, q->iono_nodes,
                         q->iono_weights) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    for (k = 0; k < 4; k++) {
        residuals[k] = in->iono_residuals[q->iono_nodes[k]];
    }
    if (ssr_area_interpolate(q->iono_weights, residuals, valid, &q->iono_grid) != 0) {
        return SSRQ_ATMOSPHERE;
    }
    q->slant_tec = q->iono_polynomial + q->iono_grid;
    for (k = 0; k < in->signals; k++) {
        if (ssr_ionosphere_delay(q->iono_polynomial, q->iono_grid, in->frequency[k], &q->code_delay[k],
                                 &q->phase_delay[k]) != 0) {
            return SSRQ_ATMOSPHERE;
        }
    }
    return SSRQ_OK;
}
