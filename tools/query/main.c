/* Host driver of one query (query.h): prints every value, or times the parts.
 *
 *   one_query            every input and result, one "name value" per line,
 *                        values with 17 significant digits (exact doubles)
 *   one_query time N     seconds per query and per epoch, best of 5 runs
 *                        of N each
 *
 * Built by tools/portability.py and tests/test_portability.py with the
 * kernels and -lm. [design: the I/O stays in this file, so that query.c and
 * scene.c build unchanged for a target without a C library's stdio] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "query.h"

static void scalar(const char *name, double v)
{
    printf("%s %.17g\n", name, v);
}

static void integer(const char *name, long v)
{
    printf("%s %ld\n", name, v);
}

static void vector(const char *name, const double *v, int n)
{
    int k;
    for (k = 0; k < n; k++) {
        printf("%s[%d] %.17g\n", name, k, v[k]);
    }
}

static void print_all(const ssrq_inputs *in, const ssrq_receiver_state *r, const ssrq_choice_state *c,
                      const ssrq_epoch_state *e, const ssrq_query_state *q)
{
    int k;
    /* inputs */
    vector("in.receiver", in->receiver, 3);
    integer("in.t_rx.week", in->t_rx.week);
    integer("in.t_rx.second", in->t_rx.second);
    scalar("in.t_rx.fraction", in->t_rx.fraction);
    scalar("in.sqrt_a", in->sqrt_a);
    scalar("in.e", in->e);
    scalar("in.m0", in->m0);
    scalar("in.delta_n", in->delta_n);
    scalar("in.omega", in->omega);
    scalar("in.omega0", in->omega0);
    scalar("in.omega_dot", in->omega_dot);
    scalar("in.i0", in->i0);
    scalar("in.idot", in->idot);
    scalar("in.cuc", in->cuc);
    scalar("in.cus", in->cus);
    scalar("in.crc", in->crc);
    scalar("in.crs", in->crs);
    scalar("in.cic", in->cic);
    scalar("in.cis", in->cis);
    scalar("in.toe", in->toe);
    scalar("in.toc", in->toc);
    scalar("in.af0", in->af0);
    scalar("in.af1", in->af1);
    scalar("in.af2", in->af2);
    vector("in.d_o", in->d_o, 3);
    scalar("in.yaw_deg", in->yaw_deg);
    scalar("in.d_c", in->d_c);
    scalar("in.hydrostatic", in->hydrostatic);
    vector("in.tropo_c", in->tropo_c, 4);
    vector("in.tropo_residuals", in->tropo_residuals, in->lat_nodes * in->lon_nodes);
    vector("in.iono_c", in->iono_c, 4);
    vector("in.iono_residuals", in->iono_residuals, in->lat_nodes * in->lon_nodes);
    /* per receiver, per choice */
    vector("receiver.geodetic", r->geodetic, 3);
    vector("receiver.east", r->east, 3);
    vector("receiver.north", r->north, 3);
    scalar("receiver.latitude_deg", r->latitude_deg);
    scalar("receiver.longitude_deg", r->longitude_deg);
    scalar("choice.m0", c->m0);
    scalar("choice.omega0", c->omega0);
    /* per epoch */
    scalar("epoch.doy", e->doy);
    vector("epoch.sun", e->sun, 3);
    vector("epoch.moon", e->moon, 3);
    scalar("epoch.tt_centuries", e->tt_centuries);
    scalar("epoch.gmst", e->gmst);
    vector("epoch.tide_displacement", e->tide_displacement, 3);
    scalar("epoch.area_distance", e->area_distance);
    scalar("epoch.tropo_polynomial", e->tropo_polynomial);
    scalar("epoch.tropo_grid", e->tropo_grid);
    scalar("epoch.height_factor", e->height_factor);
    scalar("epoch.hydrostatic_zenith", e->hydrostatic_zenith);
    scalar("epoch.wet_zenith", e->wet_zenith);
    scalar("epoch.tropo_sigma", e->tropo_sigma);
    scalar("epoch.iono_sigma", e->iono_sigma);
    /* per query */
    integer("query.iterations", q->iterations);
    integer("query.t_tx.week", q->t_tx.week);
    integer("query.t_tx.second", q->t_tx.second);
    scalar("query.t_tx.fraction", q->t_tx.fraction);
    vector("query.broadcast_position", q->broadcast_position, 3);
    vector("query.velocity", q->velocity, 3);
    scalar("query.broadcast_clock", q->broadcast_clock);
    vector("query.position", q->position, 3);
    scalar("query.range", q->range);
    scalar("query.range_rate", q->range_rate);
    scalar("query.sagnac", q->sagnac);
    vector("query.orbit_ecef", q->orbit_ecef, 3);
    scalar("query.clock_drift", q->clock_drift);
    scalar("query.clock", q->clock);
    scalar("query.ure", q->ure);
    scalar("query.tau", q->tau);
    vector("query.unit", q->unit, 3);
    scalar("query.elevation", q->elevation);
    scalar("query.azimuth", q->azimuth);
    scalar("query.shapiro", q->shapiro);
    scalar("query.relativity_clock", q->relativity_clock);
    scalar("query.solid_tide", q->solid_tide);
    scalar("query.windup", q->windup);
    scalar("query.hydrostatic_mapping", q->hydrostatic_mapping);
    scalar("query.wet_mapping", q->wet_mapping);
    scalar("query.tropo_slant", q->tropo_slant);
    scalar("query.iono_polynomial", q->iono_polynomial);
    scalar("query.iono_grid", q->iono_grid);
    scalar("query.slant_tec", q->slant_tec);
    for (k = 0; k < in->signals; k++) {
        char name[32];
        sprintf(name, "query.code_delay[%d]", k);
        scalar(name, q->code_delay[k]);
        sprintf(name, "query.phase_delay[%d]", k);
        scalar(name, q->phase_delay[k]);
    }
}

/* [design: clock(), the C library's processor time, since C99 has no wall
 * clock; the best of five runs] */
static double best_seconds(const ssrq_inputs *in, const ssrq_receiver_state *r, const ssrq_choice_state *c,
                           ssrq_epoch_state *e, ssrq_query_state *q, long n, int epoch)
{
    double best = -1.0;
    int run;
    for (run = 0; run < 5; run++) {
        long k;
        clock_t start = clock();
        for (k = 0; k < n; k++) {
            if ((epoch ? ssrq_epoch(in, r, e) : ssrq_query(in, r, c, e, q)) != SSRQ_OK) {
                return -1.0;
            }
        }
        {
            double s = (double)(clock() - start) / (double)CLOCKS_PER_SEC / (double)n;
            if (best < 0.0 || s < best) {
                best = s;
            }
        }
    }
    return best;
}

int main(int argc, char **argv)
{
    static ssrq_inputs in;
    static ssrq_receiver_state r;
    static ssrq_choice_state c;
    static ssrq_epoch_state e;
    static ssrq_query_state q;
    int32_t status;

    ssrq_scene(&in);
    status = ssrq_receiver(&in, &r);
    if (status == SSRQ_OK) {
        status = ssrq_choice(&in, &c);
    }
    if (status == SSRQ_OK) {
        status = ssrq_epoch(&in, &r, &e);
    }
    if (status == SSRQ_OK) {
        status = ssrq_query(&in, &r, &c, &e, &q);
    }
    if (status != SSRQ_OK) {
        fprintf(stderr, "one_query: step %d failed\n", (int)status);
        return 1;
    }
    if (argc == 3 && strcmp(argv[1], "time") == 0) {
        long n = strtol(argv[2], NULL, 10);
        printf("query_seconds %.6g\n", best_seconds(&in, &r, &c, &e, &q, n, 0));
        printf("epoch_seconds %.6g\n", best_seconds(&in, &r, &c, &e, &q, n, 1));
        return 0;
    }
    print_all(&in, &r, &c, &e, &q);
    return 0;
}
