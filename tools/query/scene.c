/* The fixed inputs of one query (query.h, ssrq_scene).
 *
 * They are the values the engine holds after its adapters in the public
 * test tests/test_engine.py, test_every_term_of_a_gps_record_is_the_
 * toolkits_value: the 1019 of GPS 5 built from the public broadcast record
 * of 2024-01-01 (tests/test_joins.py) and decoded at RTCM 3 resolution, and
 * the synthetic SPARTN GAD, HPAC and OCB of that test decoded at SPARTN
 * resolution; tests/test_portability.py checks each against the engine's
 * decoded input. Decimal constants are the shortest that give each double
 * exactly. */
#include "query.h"
#include "ssrcore/timescale.h"

/* [design: the scene of a public test, so that the record checks every value] */
void ssrq_scene(ssrq_inputs *in)
{
    static const double tropo_residuals[12] = {
        -0.01999999999999999, -0.016, -0.011999999999999997, -0.007999999999999993,
        -0.0040000000000000036, 0.0, 0.0040000000000000036, 0.008000000000000007,
        0.01200000000000001, 0.016000000000000014, 0.020000000000000018, 0.023999999999999994};
    static const double iono_residuals[12] = {
        -0.20000000000000018, -0.16000000000000014, -0.1200000000000001, -0.08000000000000007,
        -0.040000000000000036, 0.0, 0.040000000000000036, 0.08000000000000007,
        0.1200000000000001, 0.16000000000000014, 0.20000000000000018, 0.2400000000000002};
    int32_t k;

    in->receiver[0] = -3807024.596862735;
    in->receiver[1] = -1185756.1189071152;
    in->receiver[2] = -4961507.1109465;
    in->t_rx.week = 2295;
    in->t_rx.second = 87400;
    in->t_rx.fraction = 0.0;
    in->leap_seconds = 18.0;

    in->system = SSR_SYSTEM_GPS;
    in->toe_week = 2295;
    in->toe = 93600.0;
    in->sqrt_a = 5153.783403396606;
    in->e = 0.005790898343548179;
    in->m0 = -0.7939378307200968;
    in->delta_n = 1.4176748663885519e-09;
    in->omega = 0.39317016769200563;
    in->omega0 = -0.20053031388670206;
    in->omega_dot = -2.596038939373102e-09;
    in->i0 = 0.3085547690279782;
    in->idot = -1.2505552149377763e-11;
    in->cuc = -1.7471611499786377e-06;
    in->cus = 3.3173710107803345e-06;
    in->crc = 319.34375;
    in->crs = -34.125;
    in->cic = 1.0989606380462646e-07;
    in->cis = 4.0978193283081055e-08;
    in->toc_week = 2295;
    in->toc = 93600.0;
    in->af0 = -0.0001562642864882946;
    in->af1 = -1.4779288903810084e-12;
    in->af2 = 0.0;

    in->d_o[0] = 0.5;
    in->d_o[1] = -0.25;
    in->d_o[2] = 0.12399999999999878;
    in->yaw_deg = 60.0;
    in->d_c = 1.2339999999999982;
    in->ure_index = 2u;

    in->ref_lat = -50.0;
    in->ref_lon = -164.0;
    in->lat_nodes = 3;
    in->lon_nodes = 4;
    in->lat_spacing = 1.0;
    in->lon_spacing = 1.0;

    in->hydrostatic = 2.312;
    in->tropo_equation = 2;
    in->tropo_c[0] = 0.1;
    in->tropo_c[1] = 0.009999999999999995;
    in->tropo_c[2] = -0.019999999999999997;
    in->tropo_c[3] = 0.0040000000000000036;
    in->tropo_quality = 2u;
    in->iono_equation = 2;
    in->iono_c[0] = 10.0;
    in->iono_c[1] = 0.3999999999999986;
    in->iono_c[2] = -0.8000000000000007;
    in->iono_c[3] = 0.02000000000000135;
    in->iono_quality = 5u;
    for (k = 0; k < 64; k++) {
        in->tropo_residuals[k] = k < 12 ? tropo_residuals[k] : 0.0;
        in->iono_residuals[k] = k < 12 ? iono_residuals[k] : 0.0;
    }

    in->signals = 2;
    in->frequency[0] = 1575420000.0;     /* GPS L1C */
    in->frequency[1] = 1227600000.0;     /* GPS L2W */
    in->frequency[2] = 0.0;
    in->frequency[3] = 0.0;
}
