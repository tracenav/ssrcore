/* One query's kernel calls, end to end, in plain C99.
 *
 * The calls the engine makes for one Keplerian satellite (GPS, Galileo,
 * BeiDou MEO, QZSS) with SPARTN orbit, clock and atmosphere corrections, at
 * one receiver and one epoch, in the engine's order (src/ssrcore/engine.py
 * Engine.query, model/selection.py transmit, model/providers.py and
 * model/atmosphere.py for the SPARTN provider), and the arithmetic the engine
 * does between them in Python, written here as the engine writes it so that
 * every value is the record's to the last bit (tests/test_portability.py).
 *
 * The work is split as the engine caches it (DESIGN 2.8):
 *   ssrq_receiver  once per receiver position: geodetic, local frame;
 *   ssrq_choice    once per selection: the orbit elements to radians;
 *   ssrq_epoch     once per epoch: day of year, Sun and Moon, the solid
 *                  Earth tide, the area choice, the zenith troposphere and
 *                  the satellite's ionosphere quality;
 *   ssrq_query     every query: light time, corrected state, geometry,
 *                  relativity, wind-up, slant troposphere and ionosphere.
 *
 * No allocation, no I/O, no global state: every value is in the caller's
 * structures. [design: the artefact a firmware project takes; the host build
 * checks it against the engine and the cross build measures its size] */
#ifndef SSRCORE_QUERY_H
#define SSRCORE_QUERY_H

#include <stdint.h>

/* A GPS time as the engine keeps it (model.types.GnssTime): week, whole
 * second of week, fraction of the second in [0, 1). [design: DESIGN 4,
 * time; differences are taken in whole seconds and fractions apart] */
typedef struct {
    int32_t week;
    int32_t second;
    double fraction;
} ssrq_time;

/* The inputs of one query, as the engine holds them after its adapters. */
typedef struct {
    /* receiver, ECEF m, and the reception time (GPS) */
    double receiver[3];
    ssrq_time t_rx;
    double leap_seconds;               /* GPS - UTC, s */
    /* broadcast ephemeris (KeplerEphemeris: angles in semicircles) */
    int32_t system;                    /* SSR_SYSTEM_GPS ... */
    int32_t toe_week;
    double toe, sqrt_a, e, m0, delta_n, omega, omega0, omega_dot, i0, idot;
    double cuc, cus, crc, crs, cic, cis;
    int32_t toc_week;
    double toc, af0, af1, af2;
    /* SPARTN orbit and clock correction (OcbInput) */
    double d_o[3];                     /* radial, along, cross, m */
    double yaw_deg;                    /* SF021, degrees */
    double d_c;                        /* clock correction, m */
    uint32_t ure_index;                /* SF024 */
    /* GAD area (degrees) */
    double ref_lat, ref_lon, lat_spacing, lon_spacing;
    int32_t lat_nodes, lon_nodes;
    /* HPAC troposphere: SF043 hydrostatic, polynomial, residuals, quality */
    double hydrostatic;
    int32_t tropo_equation;
    double tropo_c[4];
    double tropo_residuals[64];
    uint32_t tropo_quality;
    /* HPAC ionosphere of the satellite */
    int32_t iono_equation;
    double iono_c[4];
    double iono_residuals[64];
    uint32_t iono_quality;
    /* signals: carrier frequencies, Hz */
    int32_t signals;
    double frequency[4];
} ssrq_inputs;

/* Once per receiver. */
typedef struct {
    double geodetic[3];                /* lat, lon (rad), height (m) */
    double east[3], north[3], up[3];
    double latitude_deg, longitude_deg;
} ssrq_receiver_state;

/* Once per selection: the kernel's orbit (radians). */
typedef struct {
    double m0, delta_n, omega, omega0, omega_dot, i0, idot;
} ssrq_choice_state;

/* Once per epoch. */
typedef struct {
    int32_t day_of_year;
    double doy;                        /* from January 0.0 */
    double sun[3], moon[3], tt_centuries, gmst;
    double tide_displacement[3];
    int32_t inside;
    double cos_lat, area_distance;     /* the area choice of rule A4 */
    double tropo_dlat, tropo_dlon, tropo_polynomial;
    int32_t tropo_nodes[4];
    double tropo_weights[4], tropo_grid;
    double temperature, lapse_rate, height_factor;
    double hydrostatic_zenith, tropo_sigma, wet_zenith;
    double iono_sigma;
} ssrq_epoch_state;

/* Every query. */
typedef struct {
    int32_t iterations;
    ssrq_time t_tx;
    double broadcast_position[3], velocity[3], ecc_anomaly, broadcast_clock;
    double position[3], range, range_rate, sagnac;
    double orbit_ecef[3];
    double clock_drift, clock;
    double ure;
    double tau, arrival[3];
    double unit[3], elevation, azimuth;
    double shapiro, relativity_clock, solid_tide;
    double yaw, frame_x[3], frame_y[3], frame_z[3], windup;
    double hydrostatic_mapping, wet_mapping, tropo_slant;
    double iono_dlat, iono_dlon, iono_polynomial;
    int32_t iono_nodes[4];
    double iono_weights[4], iono_grid, slant_tec;
    double code_delay[4], phase_delay[4];
} ssrq_query_state;

/* Status: 0, or the step that failed (its kernel's own status is not kept). */
#define SSRQ_OK 0
#define SSRQ_GEODETIC 1
#define SSRQ_ORBIT 2
#define SSRQ_LIGHT_TIME 3
#define SSRQ_GEOMETRY 4
#define SSRQ_ATMOSPHERE 5
#define SSRQ_MODELLED 6
#define SSRQ_SCALE 7

/* [spec: the engine's order of DESIGN 3.4 and 2.8] */
int32_t ssrq_receiver(const ssrq_inputs *in, ssrq_receiver_state *r);
int32_t ssrq_choice(const ssrq_inputs *in, ssrq_choice_state *c);
int32_t ssrq_epoch(const ssrq_inputs *in, const ssrq_receiver_state *r, ssrq_epoch_state *e);
int32_t ssrq_query(const ssrq_inputs *in, const ssrq_receiver_state *r, const ssrq_choice_state *c,
                   const ssrq_epoch_state *e, ssrq_query_state *q);

/* The fixed inputs of tests/test_engine.py, test_every_term_of_a_gps_record_
 * is_the_toolkits_value: GPS 5 of the public broadcast record of 2024-01-01,
 * received in GPS week 2295 at 87 400 s, with its synthetic SPARTN set.
 * [design: a public test's inputs, so that the engine's record checks every
 * value] */
void ssrq_scene(ssrq_inputs *in);

#endif
