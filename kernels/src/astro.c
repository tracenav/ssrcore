/* Positions of the Sun and the Moon in the terrestrial frame. See ssrcore/astro.h.
 *
 * The series are those of ELP 2000-82B (files ELP1 to ELP36) and VSOP87D
 * (file VSOP87D.ear), truncated [design: the tide of TN36 7.1.1 changes by
 * at most 0.33 m per radian of Moon direction, 0.15 m per radian of Sun
 * direction, and 3 times its size per unit relative distance; terms of
 * ELP 2000-82B with |A| >= 0.2 (arcsec, km) and of VSOP87D with |A| >= 1e-6
 * (rad, au) keep the Moon within 4.3 arcsec and 3.3 km and the Sun within
 * 2 arcsec and 820 km of the complete series over 1990-2060, which moves
 * the tide by under 0.02 mm and a wind-up by far less]. The selected terms
 * are copied verbatim, in file order, so that each can be checked against
 * its file. The corrections of ELP 2000-82B section 7 to the main problem
 * amplitudes are not applied [design: over 1990-2060 they change the Moon
 * by at most 0.18 arcsec and 51 m, below the truncation].
 */

#include <math.h>
#include <stdint.h>

#include "ssrcore/astro.h"
#include "ssrcore/geodesy.h"

/* [design: SSR_PI, the double nearest pi (DESIGN 4); the astronomical series
 * are not SPARTN or ICD formulae, which use 3.1415926535898] */
#define PI SSR_PI
#define TWO_PI (2.0 * PI)
/* [spec: definition of the arcsecond] radians per arcsecond. */
#define ARCSEC (PI / 648000.0)
#define DEG (PI / 180.0)

typedef struct {
    int8_t n[4];          /* multipliers of D, l', l, F */
    double a;             /* arcsec (sine) or km (cosine) */
} elp_main_term;

typedef struct {
    int8_t coordinate;    /* 0 longitude, 1 latitude, 2 distance */
    int8_t power;         /* the term is multiplied by t to this power */
    int8_t n[13];         /* multipliers of zeta, D, l', l, F, Me, V, T, Ma, J, S, U, N */
    double phase;         /* degrees */
    double a;             /* arcsec or km */
} elp_perturbation_term;

typedef struct {
    int8_t coordinate;    /* 0 longitude, 1 latitude, 2 radius */
    int8_t alpha;         /* power of the time in thousands of Julian years */
    double a, b, c;       /* A cos(B + C T): rad or au; rad; rad per millennium */
} vsop_term;

/* BEGIN GENERATED TABLES */
/* ELP1: 101 of 1023 terms. */
static const elp_main_term ELP_MAIN_LONGITUDE[101] = {
    {{  0,   0,   0,   2}, -411.60287},
    {{  0,   0,   0,   4}, 0.42034},
    {{  0,   0,   1,  -2}, 39.53393},
    {{  0,   0,   1,   0}, 22639.55},
    {{  0,   0,   1,   2}, -45.10032},
    {{  0,   0,   2,  -2}, -1.37259},
    {{  0,   0,   2,   0}, 769.02326},
    {{  0,   0,   2,   2}, -3.99767},
    {{  0,   0,   3,   0}, 36.12364},
    {{  0,   0,   3,   2}, -0.32866},
    {{  0,   0,   4,   0}, 1.93367},
    {{  0,   1,  -3,   0}, -0.6694},
    {{  0,   1,  -2,   0}, -9.67938},
    {{  0,   1,  -1,  -2}, 0.30157},
    {{  0,   1,  -1,   0}, -147.32654},
    {{  0,   1,   0,   0}, -666.44186},
    {{  0,   1,   0,   2}, 0.41342},
    {{  0,   1,   1,   0}, -109.38419},
    {{  0,   1,   1,   2}, 0.26338},
    {{  0,   1,   2,   0}, -7.63041},
    {{  0,   1,   3,   0}, -0.54594},
    {{  0,   2,  -1,   0}, -2.56813},
    {{  0,   2,   0,   0}, -7.44804},
    {{  0,   2,   1,   0}, -1.16177},
    {{  1,  -1,   0,   0}, -0.55694},
    {{  1,   0,  -2,   0}, -1.75296},
    {{  1,   0,  -1,   0}, -18.58467},
    {{  1,   0,   0,  -2}, -0.58332},
    {{  1,   0,   0,   0}, -124.98806},
    {{  1,   0,   0,   2}, 0.25429},
    {{  1,   0,   1,   0}, -8.45308},
    {{  1,   0,   2,   0}, -0.58399},
    {{  1,   1,  -2,   0}, 0.34966},
    {{  1,   1,  -1,   0}, 1.07773},
    {{  1,   1,   0,   0}, 17.95512},
    {{  1,   1,   1,   0}, 1.26186},
    {{  2,  -3,  -1,   0}, 0.24697},
    {{  2,  -3,   0,   0}, 0.33983},
    {{  2,  -2,  -2,   0}, 0.29422},
    {{  2,  -2,  -1,   0}, 7.37173},
    {{  2,  -2,   0,   0}, 8.05076},
    {{  2,  -2,   1,   0}, 0.75173},
    {{  2,  -1,  -3,   0}, 0.47842},
    {{  2,  -1,  -2,   0}, 8.60582},
    {{  2,  -1,  -1,   0}, 205.44315},
    {{  2,  -1,  -1,   2}, -0.42624},
    {{  2,  -1,   0,  -2}, 2.14619},
    {{  2,  -1,   0,   0}, 164.73458},
    {{  2,  -1,   0,   2}, -0.38215},
    {{  2,  -1,   1,  -2}, -0.37453},
    {{  2,  -1,   1,   0}, 14.53078},
    {{  2,  -1,   2,   0}, 1.17704},
    {{  2,   0,  -4,   0}, 0.94827},
    {{  2,   0,  -3,   0}, 13.194},
    {{  2,   0,  -2,  -2}, -0.56065},
    {{  2,   0,  -2,   0}, 211.65487},
    {{  2,   0,  -2,   2}, -0.53572},
    {{  2,   0,  -1,   0}, 4586.43061},
    {{  2,   0,  -1,   2}, -9.36601},
    {{  2,   0,   0,  -2}, 55.17801},
    {{  2,   0,   0,   0}, 2369.91227},
    {{  2,   0,   0,   2}, -5.7417},
    {{  2,   0,   1,  -2}, -6.38325},
    {{  2,   0,   1,   0}, 191.95575},
    {{  2,   0,   1,   2}, -0.99023},
    {{  2,   0,   2,  -2}, -0.4538},
    {{  2,   0,   2,   0}, 14.37964},
    {{  2,   0,   3,   0}, 1.05949},
    {{  2,   1,  -2,   0}, 2.48897},
    {{  2,   1,  -1,   0}, -28.3981},
    {{  2,   1,   0,  -2}, -1.43724},
    {{  2,   1,   0,   0}, -24.3591},
    {{  2,   1,   1,   0}, -2.91464},
    {{  2,   1,   2,   0}, -0.29023},
    {{  2,   2,  -2,   0}, -0.2499},
    {{  2,   2,  -1,   0}, -2.52138},
    {{  3,  -1,  -1,   0}, -0.23141},
    {{  3,   0,  -2,   0}, -1.22412},
    {{  3,   0,  -1,   0}, -3.20968},
    {{  3,   0,   0,  -2}, -0.25304},
    {{  3,   0,   0,   0}, 0.40423},
    {{  3,   1,  -1,   0}, 0.27377},
    {{  4,  -2,  -1,   0}, 0.30874},
    {{  4,  -1,  -2,   0}, 2.73198},
    {{  4,  -1,  -1,   0}, 4.37416},
    {{  4,  -1,   0,   0}, 1.87083},
    {{  4,  -1,   1,   0}, 0.28251},
    {{  4,   0,  -3,   0}, 1.18682},
    {{  4,   0,  -2,   0}, 30.77247},
    {{  4,   0,  -1,  -2}, 0.30086},
    {{  4,   0,  -1,   0}, 38.42974},
    {{  4,   0,  -1,   2}, -0.20134},
    {{  4,   0,   0,   0}, 13.89903},
    {{  4,   0,   1,   0}, 1.97772},
    {{  4,   0,   2,   0}, 0.21853},
    {{  4,   1,  -2,   0}, -0.35759},
    {{  4,   1,  -1,   0}, -0.63523},
    {{  4,   1,   0,   0}, -0.28911},
    {{  6,   0,  -3,   0}, 0.29255},
    {{  6,   0,  -2,   0}, 0.57156},
    {{  6,   0,  -1,   0}, 0.39451},
};

/* ELP2: 76 of 918 terms. */
static const elp_main_term ELP_MAIN_LATITUDE[76] = {
    {{  0,   0,   0,   1}, 18461.4},
    {{  0,   0,   0,   3}, -6.29664},
    {{  0,   0,   1,  -3}, 2.79871},
    {{  0,   0,   1,  -1}, 999.70079},
    {{  0,   0,   1,   1}, 1010.1743},
    {{  0,   0,   1,   3}, -1.01941},
    {{  0,   0,   2,  -1}, 31.75985},
    {{  0,   0,   2,   1}, 61.91229},
    {{  0,   0,   3,  -1}, 1.58131},
    {{  0,   0,   3,   1}, 3.98407},
    {{  0,   0,   4,   1}, 0.26325},
    {{  0,   1,  -2,  -1}, -0.79322},
    {{  0,   1,  -2,   1}, -0.30129},
    {{  0,   1,  -1,  -1}, -6.73173},
    {{  0,   1,  -1,   1}, -5.6326},
    {{  0,   1,   0,  -1}, -4.83983},
    {{  0,   1,   0,   1}, -6.46036},
    {{  0,   1,   1,  -1}, -5.07614},
    {{  0,   1,   1,   1}, -5.31151},
    {{  0,   1,   2,  -1}, -0.31292},
    {{  0,   1,   2,   1}, -0.63884},
    {{  1,   0,  -1,  -1}, -0.42989},
    {{  1,   0,   0,  -1}, -4.80578},
    {{  1,   0,   0,   1}, -5.36844},
    {{  1,   0,   1,  -1}, -0.58893},
    {{  1,   0,   1,   1}, -0.66741},
    {{  1,   1,   0,  -1}, 0.80426},
    {{  1,   1,   0,   1}, 0.80263},
    {{  2,  -2,  -1,  -1}, 0.26865},
    {{  2,  -2,  -1,   1}, 0.31474},
    {{  2,  -2,   0,  -1}, 1.08587},
    {{  2,  -2,   0,   1}, 0.38353},
    {{  2,  -1,  -2,  -1}, 0.65025},
    {{  2,  -1,  -1,  -1}, 7.43488},
    {{  2,  -1,  -1,   1}, 8.86853},
    {{  2,  -1,   0,  -1}, 29.57794},
    {{  2,  -1,   0,   1}, 7.95891},
    {{  2,  -1,   1,  -1}, 1.76606},
    {{  2,  -1,   1,   1}, 1.13466},
    {{  2,   0,  -3,  -1}, 1.51564},
    {{  2,   0,  -3,   1}, 0.25408},
    {{  2,   0,  -2,  -1}, 15.56635},
    {{  2,   0,  -2,   1}, -1.62443},
    {{  2,   0,  -1,  -3}, 0.32907},
    {{  2,   0,  -1,  -1}, 166.57528},
    {{  2,   0,  -1,   1}, 199.48515},
    {{  2,   0,  -1,   3}, -0.24484},
    {{  2,   0,   0,  -3}, 2.18637},
    {{  2,   0,   0,  -1}, 623.65783},
    {{  2,   0,   0,   1}, 117.26161},
    {{  2,   0,   1,  -3}, -0.29116},
    {{  2,   0,   1,  -1}, 33.35743},
    {{  2,   0,   1,   1}, 15.12165},
    {{  2,   0,   2,  -1}, 2.14618},
    {{  2,   0,   2,   1}, 1.51976},
    {{  2,   1,  -1,  -1}, -0.79105},
    {{  2,   1,  -1,   1}, -1.31788},
    {{  2,   1,   0,  -1}, -12.0947},
    {{  2,   1,   0,   1}, -1.26433},
    {{  2,   1,   1,  -1}, -0.82275},
    {{  2,   1,   1,   1}, -0.23702},
    {{  3,   0,  -1,  -1}, -0.30517},
    {{  3,   0,  -1,   1}, -0.20593},
    {{  3,   0,   0,  -1}, -0.35183},
    {{  4,  -1,  -2,   1}, 0.21376},
    {{  4,  -1,  -1,  -1}, 0.5958},
    {{  4,  -1,  -1,   1}, 0.33882},
    {{  4,  -1,   0,  -1}, 0.41496},
    {{  4,   0,  -2,  -1}, 0.63371},
    {{  4,   0,  -2,   1}, 2.41389},
    {{  4,   0,  -1,  -1}, 6.57962},
    {{  4,   0,  -1,   1}, 2.9985},
    {{  4,   0,   0,  -1}, 3.67449},
    {{  4,   0,   0,   1}, 1.19188},
    {{  4,   0,   1,  -1}, 0.47338},
    {{  4,   0,   1,   1}, 0.21259},
};

/* ELP3: 88 of 704 terms. */
static const elp_main_term ELP_MAIN_DISTANCE[88] = {
    {{  0,   0,   0,   0}, 385000.52719},
    {{  0,   0,   0,   2}, -3.14837},
    {{  0,   0,   1,  -2}, 79.66183},
    {{  0,   0,   1,   0}, -20905.32206},
    {{  0,   0,   2,  -2}, -4.42124},
    {{  0,   0,   2,   0}, -569.92332},
    {{  0,   0,   3,   0}, -23.21032},
    {{  0,   0,   4,   0}, -1.11693},
    {{  0,   1,  -3,   0}, -0.42242},
    {{  0,   1,  -2,   0}, -7.00293},
    {{  0,   1,  -1,   0}, -129.62476},
    {{  0,   1,  -1,   2}, 0.33465},
    {{  0,   1,   0,   0}, 48.8901},
    {{  0,   1,   1,  -2}, -0.24811},
    {{  0,   1,   1,   0}, 104.75896},
    {{  0,   1,   2,   0}, 5.75105},
    {{  0,   1,   3,   0}, 0.35509},
    {{  0,   2,  -1,   0}, -2.11728},
    {{  0,   2,   0,   0}, 1.06575},
    {{  0,   2,   1,   0}, 1.16562},
    {{  1,  -1,   0,   0}, 0.49757},
    {{  1,   0,  -2,   0}, -1.73852},
    {{  1,   0,  -1,   0}, -8.37909},
    {{  1,   0,   0,  -2}, -0.79564},
    {{  1,   0,   0,   0}, 108.74265},
    {{  1,   0,   1,   0}, 6.32199},
    {{  1,   0,   2,   0}, 0.37852},
    {{  1,   1,  -2,   0}, 0.33226},
    {{  1,   1,  -1,   0}, 0.85127},
    {{  1,   1,   0,   0}, -16.67533},
    {{  1,   1,   1,   0}, -0.93335},
    {{  2,  -3,   0,   0}, -0.41076},
    {{  2,  -2,  -2,   0}, 0.34304},
    {{  2,  -2,  -1,   0}, -4.95049},
    {{  2,  -2,   0,   0}, -9.88519},
    {{  2,  -2,   1,   0}, -0.65758},
    {{  2,  -1,  -3,   0}, 0.49506},
    {{  2,  -1,  -2,   0}, 10.05654},
    {{  2,  -1,  -1,  -2}, 0.32336},
    {{  2,  -1,  -1,   0}, -152.14314},
    {{  2,  -1,   0,  -2}, 0.6571},
    {{  2,  -1,   0,   0}, -204.59357},
    {{  2,  -1,   1,  -2}, 0.20942},
    {{  2,  -1,   1,   0}, -12.83185},
    {{  2,  -1,   2,   0}, -0.84883},
    {{  2,   0,  -4,   0}, 0.77854},
    {{  2,   0,  -3,   0}, 14.40262},
    {{  2,   0,  -2,  -2}, 0.47263},
    {{  2,   0,  -2,   0}, 246.15768},
    {{  2,   0,  -2,   2}, 0.77405},
    {{  2,   0,  -1,  -2}, 8.7517},
    {{  2,   0,  -1,   0}, -3699.10468},
    {{  2,   0,  -1,   2}, 0.59633},
    {{  2,   0,   0,  -2}, 10.32129},
    {{  2,   0,   0,   0}, -2955.96651},
    {{  2,   0,   1,  -2}, 4.13118},
    {{  2,   0,   1,   0}, -170.73274},
    {{  2,   0,   2,  -2}, 0.28399},
    {{  2,   0,   2,   0}, -10.44472},
    {{  2,   0,   3,   0}, -0.66968},
    {{  2,   1,  -1,   0}, 24.20935},
    {{  2,   1,   0,   0}, 30.82498},
    {{  2,   1,   1,   0}, 2.6165},
    {{  2,   1,   2,   0}, 0.21252},
    {{  2,   2,  -1,   0}, 2.3538},
    {{  3,  -1,  -1,   0}, 0.25561},
    {{  3,   0,  -2,   0}, 0.86243},
    {{  3,   0,  -1,   0}, 3.25823},
    {{  3,   0,   0,  -2}, 0.20099},
    {{  3,   0,   0,   0}, -1.41893},
    {{  3,   1,  -1,   0}, -0.21259},
    {{  4,  -2,  -1,   0}, -0.27906},
    {{  4,  -1,  -2,   0}, -1.8971},
    {{  4,  -1,  -1,   0}, -3.95812},
    {{  4,  -1,   0,   0}, -1.57145},
    {{  4,  -1,   1,   0}, -0.20286},
    {{  4,   0,  -3,   0}, -0.51423},
    {{  4,   0,  -2,   0}, -21.63627},
    {{  4,   0,  -1,  -2}, -0.32176},
    {{  4,   0,  -1,   0}, -34.78245},
    {{  4,   0,   0,  -2}, -0.50793},
    {{  4,   0,   0,   0}, -11.64993},
    {{  4,   0,   1,   0}, -1.42255},
    {{  4,   1,  -2,   0}, 0.23696},
    {{  4,   1,  -1,   0}, 0.57881},
    {{  4,   1,   0,   0}, 0.24453},
    {{  6,   0,  -2,   0}, -0.4225},
    {{  6,   0,  -1,   0}, -0.28663},
};

/* ELP4 to ELP36: 50 terms. */
static const elp_perturbation_term ELP_PERTURBATIONS[50] = {
    {0, 0, { 1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00127, 0.49141}, /* ELP4 */
    {0, 0, { 1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00094, 7.06304}, /* ELP4 */
    {0, 0, { 1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00071, 0.36061}, /* ELP4 */
    {0, 0, { 1,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00127, 0.49331}, /* ELP4 */
    {1, 0, { 1, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00072, 0.29854}, /* ELP5 */
    {1, 0, { 1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00075, 0.45586}, /* ELP5 */
    {1, 0, { 1,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0}, 0.00086, 0.32622}, /* ELP5 */
    {1, 0, { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 180.00071, 8.04508}, /* ELP5 */
    {1, 0, { 1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 180.00069, 0.41571}, /* ELP5 */
    {2, 0, { 1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, 90.00128, 0.45276}, /* ELP6 */
    {2, 0, { 1,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0}, 270.00126, 0.45648}, /* ELP6 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0}, 318.13776, 0.32455}, /* ELP10 */
    {0, 0, { 0, -2,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0}, 1.74333, 0.21118}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0}, 1.2289, 0.6388}, /* ELP10 */
    {0, 0, { 0,  2,  0, -2,  0,  0,  0,  2,  0, -3,  0,  0,  0}, 10.03835, 0.28281}, /* ELP10 */
    {0, 0, { 0,  2,  0, -1,  0,  0,  0,  2,  0, -3,  0,  0,  0}, 10.07001, 0.44532}, /* ELP10 */
    {0, 0, { 0,  2,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0}, 0.08642, 0.24515}, /* ELP10 */
    {0, 0, { 0,  2,  0, -1,  0,  0,  0,  2,  0, -2,  0,  0,  0}, 180.11977, 1.14307}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0}, 285.98707, 0.90114}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0}, 180.00988, 0.82155}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0}, 269.95393, 0.34355}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0}, 0.20448, 0.30155}, /* ELP10 */
    {0, 0, { 0,  2,  0, -1,  0,  0,  3, -3,  0,  0,  0,  0,  0}, 179.98144, 0.64371}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0,  8, -13,  0,  0,  0,  0,  0}, 235.75013, 0.22653}, /* ELP10 */
    {0, 0, { 0,  0,  0, -1,  0,  0, 10, -3,  0,  0,  0,  0,  0}, 333.30551, 0.56341}, /* ELP10 */
    {0, 0, { 0,  0,  0, -2,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 26.54324, 0.78811}, /* ELP10 */
    {0, 0, { 0,  0,  0, -1,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 26.54261, 14.24883}, /* ELP10 */
    {0, 0, { 0,  0,  0,  0,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 26.5456, 0.7393}, /* ELP10 */
    {1, 0, { 0,  0,  0, -1, -1,  0, 18, -16,  0,  0,  0,  0,  0}, 26.54272, 0.63014}, /* ELP11 */
    {1, 0, { 0,  0,  0, -1,  1,  0, 18, -16,  0,  0,  0,  0,  0}, 26.54287, 0.63037}, /* ELP11 */
    {2, 0, { 0, -2,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0}, 91.79862, 0.20497}, /* ELP12 */
    {2, 0, { 0,  2,  0, -1,  0,  0,  0,  2,  0, -3,  0,  0,  0}, 280.06924, 0.41011}, /* ELP12 */
    {2, 0, { 0,  2,  0, -1,  0,  0,  0,  2,  0, -2,  0,  0,  0}, 90.11969, 1.0587}, /* ELP12 */
    {2, 0, { 0,  2,  0, -1,  0,  0,  3, -3,  0,  0,  0,  0,  0}, 89.98187, 0.59827}, /* ELP12 */
    {2, 0, { 0, -2,  0, -1,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 116.54222, 0.20473}, /* ELP12 */
    {2, 0, { 0,  0,  0, -2,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 116.54311, 0.72783}, /* ELP12 */
    {2, 0, { 0,  0,  0,  0,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 296.54574, 0.68256}, /* ELP12 */
    {2, 0, { 0,  2,  0, -1,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 296.54299, 0.20367}, /* ELP12 */
    {0, 1, { 0,  0,  0, -1,  0,  0, 18, -16,  0,  0,  0,  0,  0}, 114.5655, 0.25425}, /* ELP13 */
    {0, 0, { 0,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0,  0,  0}, 95.13523, 0.28938}, /* ELP16 */
    {1, 0, { 0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0}, 275.13226, 1.37497}, /* ELP17 */
    {0, 1, { 0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 0.0, 0.37115}, /* ELP19 */
    {0, 1, { 0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 0.0, 1.6768}, /* ELP19 */
    {0, 1, { 0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 0.0, 0.2756}, /* ELP19 */
    {0, 1, { 0,  2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 180.0, 0.51642}, /* ELP19 */
    {0, 1, { 0,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 180.0, 0.41383}, /* ELP19 */
    {2, 1, { 0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 90.0, 0.32654}, /* ELP21 */
    {2, 1, { 0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 270.0, 0.26396}, /* ELP21 */
    {2, 1, { 0,  2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 90.0, 0.38245}, /* ELP21 */
    {2, 1, { 0,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}, 90.0, 0.51395}, /* ELP21 */
};

/* VSOP87D Earth: 69 of 2425 terms. */
static const vsop_term VSOP87D_EARTH[69] = {
    {0, 0, 1.75347045673, 0.0, 0.0},
    {0, 0, 0.03341656456, 4.66925680417, 6283.0758499914},
    {0, 0, 0.00034894275, 4.62610241759, 12566.1516999828},
    {0, 0, 3.417571e-05, 2.82886579606, 3.523118349},
    {0, 0, 3.497056e-05, 2.74411800971, 5753.3848848968},
    {0, 0, 3.135896e-05, 3.62767041758, 77713.7714681205},
    {0, 0, 2.676218e-05, 4.41808351397, 7860.4193924392},
    {0, 0, 2.342687e-05, 6.13516237631, 3930.2096962196},
    {0, 0, 1.273166e-05, 2.03709655772, 529.6909650946},
    {0, 0, 1.324292e-05, 0.74246356352, 11506.7697697936},
    {0, 0, 9.01855e-06, 2.04505443513, 26.2983197998},
    {0, 0, 1.199167e-05, 1.10962944315, 1577.3435424478},
    {0, 0, 8.57223e-06, 3.50849156957, 398.1490034082},
    {0, 0, 7.79786e-06, 1.17882652114, 5223.6939198022},
    {0, 0, 9.9025e-06, 5.23268129594, 5884.9268465832},
    {0, 0, 7.53141e-06, 2.53339053818, 5507.5532386674},
    {0, 0, 5.05264e-06, 4.58292563052, 18849.2275499742},
    {0, 0, 4.92379e-06, 4.20506639861, 775.522611324},
    {0, 0, 3.56655e-06, 2.91954116867, 0.0673103028},
    {0, 0, 2.84125e-06, 1.89869034186, 796.2980068164},
    {0, 0, 2.4281e-06, 0.34481140906, 5486.777843175},
    {0, 0, 3.17087e-06, 5.84901952218, 11790.6290886588},
    {0, 0, 2.71039e-06, 0.31488607649, 10977.078804699},
    {0, 0, 2.0616e-06, 4.80646606059, 2544.3144198834},
    {0, 0, 2.05385e-06, 1.86947813692, 5573.1428014331},
    {0, 0, 2.02261e-06, 2.45767795458, 6069.7767545534},
    {0, 0, 1.26184e-06, 1.0830263021, 20.7753954924},
    {0, 0, 1.55516e-06, 0.83306073807, 213.299095438},
    {0, 0, 1.15132e-06, 0.64544911683, 0.9803210682},
    {0, 0, 1.02851e-06, 0.63599846727, 4694.0029547076},
    {0, 0, 1.01724e-06, 4.26679821365, 7.1135470008},
    {0, 0, 1.32212e-06, 3.41118275555, 2942.4634232916},
    {0, 0, 1.01895e-06, 0.97569221824, 15720.8387848784},
    {0, 1, 6283.31966747491, 0.0, 0.0},
    {0, 1, 0.00206058863, 2.67823455584, 6283.0758499914},
    {0, 1, 4.30343e-05, 2.63512650414, 12566.1516999828},
    {0, 1, 4.25264e-06, 1.59046980729, 3.523118349},
    {0, 1, 1.08977e-06, 2.96618001993, 1577.3435424478},
    {0, 1, 1.19261e-06, 5.79557487799, 26.2983197998},
    {0, 2, 0.0005291887, 0.0, 0.0},
    {0, 2, 8.719837e-05, 1.07209665242, 6283.0758499914},
    {0, 2, 3.09125e-06, 0.86728818832, 12566.1516999828},
    {0, 3, 2.89226e-06, 5.84384198723, 6283.0758499914},
    {0, 4, 1.14084e-06, 3.14159265359, 0.0},
    {1, 0, 2.7962e-06, 3.19870156017, 84334.66158130829},
    {1, 0, 1.01643e-06, 5.42248619256, 5507.5532386674},
    {2, 0, 1.00013988799, 0.0, 0.0},
    {2, 0, 0.01670699626, 3.09846350771, 6283.0758499914},
    {2, 0, 0.00013956023, 3.0552460962, 12566.1516999828},
    {2, 0, 3.08372e-05, 5.19846674381, 77713.7714681205},
    {2, 0, 1.628461e-05, 1.17387749012, 5753.3848848968},
    {2, 0, 1.575568e-05, 2.84685245825, 7860.4193924392},
    {2, 0, 9.24799e-06, 5.45292234084, 11506.7697697936},
    {2, 0, 5.42444e-06, 4.56409149777, 3930.2096962196},
    {2, 0, 4.7211e-06, 3.66100022149, 5884.9268465832},
    {2, 0, 3.2878e-06, 5.89983646482, 5223.6939198022},
    {2, 0, 3.45983e-06, 0.96368617687, 5507.5532386674},
    {2, 0, 3.06784e-06, 0.29867139512, 5573.1428014331},
    {2, 0, 1.74844e-06, 3.01193636534, 18849.2275499742},
    {2, 0, 2.43189e-06, 4.27349536153, 11790.6290886588},
    {2, 0, 2.11829e-06, 5.84714540314, 1577.3435424478},
    {2, 0, 1.85752e-06, 5.02194447178, 10977.078804699},
    {2, 0, 1.09835e-06, 5.05510636285, 5486.777843175},
    {2, 1, 0.00103018608, 1.10748969588, 6283.0758499914},
    {2, 1, 1.721238e-05, 1.06442301418, 12566.1516999828},
    {2, 1, 7.02215e-06, 3.14159265359, 0.0},
    {2, 2, 4.359385e-05, 5.78455133738, 6283.0758499914},
    {2, 2, 1.23633e-06, 5.57934722157, 12566.1516999828},
    {2, 3, 1.44595e-06, 4.27319435148, 6283.0758499914},
};
/* END GENERATED TABLES */

#define COUNT(table) ((int32_t)(sizeof(table) / sizeof((table)[0])))

/* [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4; TM-23 3.1.1; TN36 ch.
 * 5 sections 5.5.7, 5.7.1] J2000.0 TT, JD 2451545.0, is 7300.5 days after the
 * GPS epoch (JD 2444244.5) on the TT scale, i.e. at GPS time 630763200 s -
 * 51.184 s. */
double ssr_astro_tt_centuries(int32_t week, double seconds)
{
    double since = (double)((int64_t)week * 604800 - 630763200) + (seconds + SSR_TT_MINUS_GPS);
    return since / 3155760000.0;
}

/* [spec: TN36 ch. 5 eq. (5.15)] ERA = 2 pi (UT1 Julian day fraction +
 * 0.7790572732640 + 0.00273781191135448 Tu); Tu = JD(UT1) - 2451545.0, and
 * the GPS week origin is JD 2444244.5, so Tu counts from GPS second
 * 630763200 of the UT1 scale. The Julian day fraction is that of Tu. */
double ssr_astro_earth_rotation_angle(int32_t week, double seconds, double ut1_minus_gps)
{
    double since = (double)((int64_t)week * 604800 - 630763200) + (seconds + ut1_minus_gps);
    double tu = since / 86400.0;
    double fraction = fmod(since, 86400.0) / 86400.0;
    double era = fmod(fraction + 0.7790572732640 + 0.00273781191135448 * tu, 1.0);
    if (era < 0.0) {
        era += 1.0;
    }
    return TWO_PI * era;
}

/* [spec: TN36 ch. 5 eq. (5.32)] GMST = ERA(UT1) + 0.014506" + 4612.156534" t
 * + 1.3915817" t^2 - 0.00000044" t^3 - 0.000029956" t^4 - 0.0000000368" t^5,
 * t in TT centuries. */
double ssr_astro_gmst(int32_t week, double seconds, double ut1_minus_gps)
{
    double t = ssr_astro_tt_centuries(week, seconds);
    double p = 0.014506 + t * (4612.156534 + t * (1.3915817 + t * (-0.00000044
               + t * (-0.000029956 + t * -0.0000000368))));
    double g = fmod(ssr_astro_earth_rotation_angle(week, seconds, ut1_minus_gps) + p * ARCSEC, TWO_PI);
    return g < 0.0 ? g + TWO_PI : g;
}

/* [spec: TN36 ch. 5 eq. (5.40)] eps_A = 84381.406" - 46.836769" t
 * - 0.0001831" t^2 + 0.0020034" t^3 - 0.000000576" t^4 - 0.0000000434" t^5. */
double ssr_astro_mean_obliquity(double t)
{
    double e = 84381.406 + t * (-46.836769 + t * (-0.0001831 + t * (0.0020034
               + t * (-0.000000576 + t * -0.0000000434))));
    return e * ARCSEC;
}

typedef struct {
    int8_t n[5];          /* multipliers of l, l', F, D, Omega */
    double in_phase, in_phase_rate, out_of_phase, out_of_phase_rate;   /* microarcsec, per century */
} nutation_term;

/* [spec: TN36 ch. 5 Table 5.3a, the extract printed there] Nutation in
 * longitude: A_i (sine) and A''_i (cosine), with A'_i and A'''_i for the
 * three terms whose rates the extract gives (rows 1321 to 1323). */
static const nutation_term NUTATION_LONGITUDE[5] = {
    {{0, 0, 0, 0, 1}, -17206424.18, -17418.82, 3338.60, 2.89},
    {{0, 0, 2, -2, 2}, -1317091.22, -163.84, -1369.60, 1.20},
    {{0, 0, 2, 0, 2}, -227641.81, 0.0, 279.60, 0.0},
    {{0, 0, 0, 0, 2}, 207455.50, 0.0, -69.80, 0.0},
    {{0, 1, 0, 0, 0}, 147587.77, -363.71, 1181.70, -1.50},
};

/* [spec: TN36 ch. 5 Table 5.3b, the extract printed there] Nutation in
 * obliquity: B_i (cosine) and B''_i (sine), with B'_i and B'''_i for the
 * three terms whose rates the extract gives (rows 1038 to 1040). */
static const nutation_term NUTATION_OBLIQUITY[5] = {
    {{0, 0, 0, 0, 1}, 9205233.10, 883.03, 1537.70, 0.20},
    {{0, 0, 2, -2, 2}, 573033.60, -303.09, -458.70, -0.30},
    {{0, 0, 2, 0, 2}, 97846.10, 0.0, 137.40, 0.0},
    {{0, 0, 0, 0, 2}, -89749.20, 0.0, -29.10, 0.0},
    {{0, 1, 2, -2, 2}, 22438.60, -67.76, -17.40, 0.00},
};

/* [spec: TN36 ch. 5 eq. (5.43)] The Delaunay arguments l, l', F, D, Omega
 * (rad). */
static void delaunay(double t, double f[5])
{
    f[0] = 134.96340251 * DEG + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * -0.00024470))) * ARCSEC;
    f[1] = 357.52910918 * DEG + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * -0.00001149))) * ARCSEC;
    f[2] = 93.27209062 * DEG + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * 0.00000417))) * ARCSEC;
    f[3] = 297.85019547 * DEG + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * -0.00003169))) * ARCSEC;
    f[4] = 125.04455501 * DEG + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * -0.00005939))) * ARCSEC;
}

/* [spec: TN36 ch. 5 eq. (5.35), with the extracts of Tables 5.3a and 5.3b;
 * design: the five largest terms of each, which TN36 prints, leave out
 * terms of at most 0.08 arcsec each (0.35 arcsec together in longitude,
 * 0.1 arcsec in obliquity), far below the 4.3 arcsec of the Moon series] */
void ssr_astro_nutation(double t, double *dpsi, double *deps)
{
    double f[5];
    double psi = 0.0;
    double eps = 0.0;
    int32_t i;
    int32_t j;
    delaunay(t, f);
    for (i = 0; i < 5; i++) {
        const nutation_term *p = &NUTATION_LONGITUDE[i];
        const nutation_term *e = &NUTATION_OBLIQUITY[i];
        double ap = 0.0;
        double ae = 0.0;
        for (j = 0; j < 5; j++) {
            ap += p->n[j] * f[j];
            ae += e->n[j] * f[j];
        }
        psi += (p->in_phase + p->in_phase_rate * t) * sin(ap)
             + (p->out_of_phase + p->out_of_phase_rate * t) * cos(ap);
        eps += (e->in_phase + e->in_phase_rate * t) * cos(ae)
             + (e->out_of_phase + e->out_of_phase_rate * t) * sin(ae);
    }
    *dpsi = psi * 1e-6 * ARCSEC;
    *deps = eps * 1e-6 * ARCSEC;
}

/* [spec: ELP 2000-82B section 7] a polynomial in t of arcseconds, to rad. */
static double arcsec_polynomial(const double c[5], double t)
{
    return (c[0] + t * (c[1] + t * (c[2] + t * (c[3] + t * c[4])))) * ARCSEC;
}

/* [spec: ELP 2000-82B section 7] W1, D, l', l, F with the corrections of the
 * fit to DE200/LE200, arcsec (the notice's degrees, minutes and seconds of
 * the constant terms written as seconds). */
static const double ELP_W1[5] = {785939.95571, 1732559343.73604, -5.8883, 0.006604, -0.00003169};
static const double ELP_D[5] = {1072260.73512, 1602961601.4603, -5.8681, 0.006595, -0.00003184};
static const double ELP_LP[5] = {1287104.79306, 129596581.0474, -0.5529, 0.000147, 0.0};
static const double ELP_L[5] = {485868.28096, 1717915923.4728, 32.3893, 0.051651, -0.00024470};
static const double ELP_F[5] = {335779.55755, 1739527263.0983, -12.2505, -0.001021, 0.00000417};
/* [spec: ELP 2000-82B section 7] T, the linear part used in ELP10 to ELP21. */
static const double ELP_T[2] = {361679.22059, 129597742.2758};
/* [spec: ELP 2000-82B section 6 b)] precession constant p in J2000, arcsec per century. */
static const double ELP_P = 5029.0966;
/* [spec: ELP 2000-82B Table F] Me, V, Ma, J, S, U, N: longitude in J2000 and
 * mean motion, arcsec and arcsec per century. */
static const double ELP_PLANETS[7][2] = {
    {908103.25986, 538101628.68898},
    {655127.28305, 210664136.43355},
    {1279559.78866, 68905077.59284},
    {123665.34212, 10925660.42861},
    {180278.89694, 4399609.65932},
    {1130598.01841, 1542481.19393},
    {1095655.19575, 786550.32074},
};
/* [spec: ELP 2000-82B section 8] accumulated precession p_A, arcsec. */
static const double ELP_PA[5] = {0.0, 5029.0966, 1.1120, 0.000077, -0.00002353};

/* [spec: ELP 2000-82B sections 2.2, 3, 4 and 8] V = W1 + the longitude
 * series, U = the latitude series, r = the distance series; the longitude
 * of date is V + p_A. Main problem (ELP1 to ELP3) with the full arguments;
 * the other files with their linear parts and zeta = W1 + p t. */
void ssr_astro_moon_ecliptic(double t, double out[3])
{
    double d = arcsec_polynomial(ELP_D, t);
    double lp = arcsec_polynomial(ELP_LP, t);
    double l = arcsec_polynomial(ELP_L, t);
    double f = arcsec_polynomial(ELP_F, t);
    double linear[13];
    double sum[3] = {0.0, 0.0, 0.0};
    double tpow[3];
    int32_t i;
    int32_t j;

    for (i = 0; i < COUNT(ELP_MAIN_LONGITUDE); i++) {
        const elp_main_term *m = &ELP_MAIN_LONGITUDE[i];
        sum[0] += m->a * sin(m->n[0] * d + m->n[1] * lp + m->n[2] * l + m->n[3] * f);
    }
    for (i = 0; i < COUNT(ELP_MAIN_LATITUDE); i++) {
        const elp_main_term *m = &ELP_MAIN_LATITUDE[i];
        sum[1] += m->a * sin(m->n[0] * d + m->n[1] * lp + m->n[2] * l + m->n[3] * f);
    }
    for (i = 0; i < COUNT(ELP_MAIN_DISTANCE); i++) {
        const elp_main_term *m = &ELP_MAIN_DISTANCE[i];
        sum[2] += m->a * cos(m->n[0] * d + m->n[1] * lp + m->n[2] * l + m->n[3] * f);
    }

    linear[0] = (ELP_W1[0] + (ELP_W1[1] + ELP_P) * t) * ARCSEC;
    linear[1] = (ELP_D[0] + ELP_D[1] * t) * ARCSEC;
    linear[2] = (ELP_LP[0] + ELP_LP[1] * t) * ARCSEC;
    linear[3] = (ELP_L[0] + ELP_L[1] * t) * ARCSEC;
    linear[4] = (ELP_F[0] + ELP_F[1] * t) * ARCSEC;
    linear[5] = (ELP_PLANETS[0][0] + ELP_PLANETS[0][1] * t) * ARCSEC;
    linear[6] = (ELP_PLANETS[1][0] + ELP_PLANETS[1][1] * t) * ARCSEC;
    linear[7] = (ELP_T[0] + ELP_T[1] * t) * ARCSEC;
    for (j = 2; j < 7; j++) {
        linear[6 + j] = (ELP_PLANETS[j][0] + ELP_PLANETS[j][1] * t) * ARCSEC;
    }
    tpow[0] = 1.0;
    tpow[1] = t;
    tpow[2] = t * t;
    for (i = 0; i < COUNT(ELP_PERTURBATIONS); i++) {
        const elp_perturbation_term *p = &ELP_PERTURBATIONS[i];
        double arg = p->phase * DEG;
        for (j = 0; j < 13; j++) {
            arg += p->n[j] * linear[j];
        }
        sum[p->coordinate] += p->a * tpow[p->power] * sin(arg);
    }

    out[0] = arcsec_polynomial(ELP_W1, t) + sum[0] * ARCSEC + arcsec_polynomial(ELP_PA, t);
    out[1] = sum[1] * ARCSEC;
    out[2] = sum[2] * 1000.0;
}

/* [spec: VSOP87 "CONTENTS" form (2) and "COMPUTATION"] each coordinate is
 * the sum over its series of T^alpha A cos(B + C T), T in thousands of
 * Julian years (t / 10). */
void ssr_astro_earth_heliocentric(double t, double out[3])
{
    double tau = t / 10.0;
    double tpow[6];
    int32_t i;
    tpow[0] = 1.0;
    for (i = 1; i < 6; i++) {
        tpow[i] = tpow[i - 1] * tau;
    }
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    for (i = 0; i < COUNT(VSOP87D_EARTH); i++) {
        const vsop_term *v = &VSOP87D_EARTH[i];
        out[v->coordinate] += tpow[v->alpha] * (v->a * cos(v->b + v->c * tau));
    }
}

/* [spec: VSOP87 version D, heliocentric Earth; TN36 ch. 1 Table 1.1, au]
 * the geocentric Sun is the heliocentric Earth reversed: longitude + pi,
 * latitude negated, the same distance. */
void ssr_astro_sun_ecliptic(double t, double out[3])
{
    double earth[3];
    ssr_astro_earth_heliocentric(t, earth);
    out[0] = fmod(earth[0], TWO_PI) + PI;
    out[1] = -earth[1];
    out[2] = earth[2] * SSR_IERS_AU;
}

/* [spec: TN36 ch. 5 section 5.4.5] rectangular ecliptic coordinates with the
 * longitude moved by dpsi, then rotated about x by -(eps_a + deps): the
 * classical nutation from the mean ecliptic of date to the true equator. */
void ssr_astro_ecliptic_to_true_equator(const double lbr[3], double dpsi, double eps_a,
                                  double deps, double out[3])
{
    double lon = lbr[0] + dpsi;
    double eps = eps_a + deps;
    double x = lbr[2] * cos(lbr[1]) * cos(lon);
    double y = lbr[2] * cos(lbr[1]) * sin(lon);
    double z = lbr[2] * sin(lbr[1]);
    out[0] = x;
    out[1] = y * cos(eps) - z * sin(eps);
    out[2] = y * sin(eps) + z * cos(eps);
}

/* [spec: TN36 ch. 5 sections 5.4.3 and 5.4.1, R3(GST); design: polar motion
 * left out, see astro.h] */
void ssr_astro_true_equator_to_ecef(const double v[3], double gst, double out[3])
{
    double c = cos(gst);
    double s = sin(gst);
    out[0] = c * v[0] + s * v[1];
    out[1] = -s * v[0] + c * v[1];
    out[2] = v[2];
}

/* [spec: ELP 2000-82B; VSOP87; TN36 ch. 5 eqs. (5.31), (5.32), (5.35),
 * (5.40), the functions above composed; design: the equation of the
 * equinoxes without its complementary terms, which are below 0.003 arcsec] */
void ssr_astro_sun_moon_ecef(int32_t week, double seconds, double ut1_minus_gps,
                             double sun[3], double moon[3])
{
    double t = ssr_astro_tt_centuries(week, seconds);
    double eps_a = ssr_astro_mean_obliquity(t);
    double dpsi;
    double deps;
    double gst;
    double lbr[3];
    double v[3];
    ssr_astro_nutation(t, &dpsi, &deps);
    gst = ssr_astro_gmst(week, seconds, ut1_minus_gps) + dpsi * cos(eps_a);

    ssr_astro_sun_ecliptic(t, lbr);
    ssr_astro_ecliptic_to_true_equator(lbr, dpsi, eps_a, deps, v);
    ssr_astro_true_equator_to_ecef(v, gst, sun);

    ssr_astro_moon_ecliptic(t, lbr);
    ssr_astro_ecliptic_to_true_equator(lbr, dpsi, eps_a, deps, v);
    ssr_astro_true_equator_to_ecef(v, gst, moon);
}
