/* Conversion of encoded SPARTN field values. See ssrcore/scale.h. */

#include "ssrcore/scale.h"

/* [spec: SPARTN v2.0.3 section 6.3] Vdecoded = (Vencoded * Res) + Rngmin. */
double ssr_scale_decode(uint32_t encoded, double resolution, double range_min)
{
    return (double)encoded * resolution + range_min;
}

/* No special value. [design: every special value of a ranged field in
 * Table 6.2 is an all-ones code of at most 14 bits, so a code no field can
 * have marks its absence] */
#define NONE 0xFFFFFFFFu

/* A ranged field of Table 6.2. top is the encoded value of the upper end of
 * the range, (Rngmax - Rngmin) / Res, written out as the integer the table
 * implies. [design: an integer bound decides "inside the range" exactly,
 * without comparing rounded doubles] */
typedef struct {
    uint32_t field;
    uint32_t bits;
    uint32_t top;
    uint32_t special;
    double resolution;
    double range_min;
} ranged_field;

/* [spec: SPARTN v2.0.3 Table 6.2] Every field with a numeric range and a
 * numeric resolution: identifiers, counts, satellite corrections, biases,
 * yaw, area geometry, troposphere and ionosphere polynomial coefficients and
 * residuals, BPAC VTEC. Enumerations, flags, masks and the key and
 * authentication fields (SF084 to SF092, which are gnsswire's) are codes,
 * not converted. A range "a to b" has Rngmin a; a range "±x" has Rngmin -x. */
static const ranged_field RANGED[] = {
    /* field bits  top    special  Res      Rngmin */
    {5,   9,  511,   NONE,   1.0,     0.0},      /* SIOU, 0 to 511 */
    {18,  8,  255,   NONE,   1.0,     0.0},      /* GPS IODE, 0 to 255 */
    {19,  7,  127,   NONE,   1.0,     0.0},      /* GLONASS IODE, 0 to 127 */
    {20,  14, 16382, NONE,   0.002,   -16.382},  /* satellite corrections, ±16.382 m */
    {21,  6,  59,    0x3Fu,  6.0,     0.0},      /* satellite yaw, 0 to 354 degrees */
    {29,  11, 2046,  NONE,   0.02,    -20.46},   /* code bias, ±20.46 m */
    {30,  5,  31,    NONE,   1.0,     1.0},      /* area count, 1 to 32 */
    {31,  8,  255,   NONE,   1.0,     0.0},      /* area ID, 0 to 255 */
    {32,  11, 1800,  NONE,   0.1,     -90.0},    /* area reference latitude, ±90.00 degrees */
    {33,  12, 3600,  NONE,   0.1,     -180.0},   /* area reference longitude, ±180.00 degrees */
    {34,  3,  7,     NONE,   1.0,     1.0},      /* area latitude node count, 1 to 8 */
    {35,  3,  7,     NONE,   1.0,     1.0},      /* area longitude node count, 1 to 8 */
    {36,  5,  31,    NONE,   0.1,     0.1},      /* area latitude node spacing, 0.1 to 3.2 degrees */
    {37,  5,  31,    NONE,   0.1,     0.1},      /* area longitude node spacing, 0.1 to 3.2 degrees */
    {39,  7,  127,   NONE,   1.0,     0.0},      /* grid points present, 0 to 127 */
    {43,  8,  254,   NONE,   0.004,   -0.508},   /* area average hydrostatic delay, ±0.508 m (from 2.3 m) */
    {45,  7,  126,   NONE,   0.004,   -0.252},   /* small T00, ±0.252 m (offset 0.252 m) */
    {46,  7,  126,   NONE,   0.001,   -0.063},   /* small T10/T01, ±0.063 m/degree */
    {47,  9,  510,   NONE,   0.0002,  -0.0510},  /* small T11, ±0.0510 m/degree^2 */
    {48,  9,  510,   NONE,   0.004,   -1.020},   /* large T00, ±1.020 m (offset 0.252 m) */
    {49,  9,  510,   NONE,   0.001,   -0.255},   /* large T10/T01, ±0.255 m/degree */
    {50,  11, 2046,  NONE,   0.0002,  -0.2046},  /* large T11, ±0.2046 m/degree^2 */
    {52,  6,  62,    0x3Fu,  0.004,   -0.124},   /* small troposphere residual, ±0.124 m */
    {53,  8,  254,   0xFFu,  0.004,   -0.508},   /* large troposphere residual, ±0.508 m */
    {57,  12, 4094,  NONE,   0.04,    -81.88},   /* small C00, ±81.88 TECU */
    {58,  12, 4094,  NONE,   0.008,   -16.376},  /* small C10/C01, ±16.376 TECU/degree */
    {59,  13, 8190,  NONE,   0.002,   -8.190},   /* small C11, ±8.190 TECU/degree^2 */
    {60,  14, 16382, NONE,   0.04,    -327.64},  /* large C00, ±327.64 TECU */
    {61,  14, 16382, NONE,   0.008,   -65.528},  /* large C10/C01, ±65.528 TECU/degree */
    {62,  15, 32766, NONE,   0.002,   -32.766},  /* large C11, ±32.766 TECU/degree^2 */
    {64,  4,  14,    0xFu,   0.04,    -0.28},    /* small ionosphere residual, ±0.28 TECU */
    {65,  7,  126,   0x7Fu,  0.04,    -2.52},    /* medium ionosphere residual, ±2.52 TECU */
    {66,  10, 1022,  0x3FFu, 0.04,    -20.44},   /* large ionosphere residual, ±20.44 TECU */
    {67,  14, 16382, 0x3FFFu, 0.04,   -327.64},  /* extra-large ionosphere residual, ±327.64 TECU */
    {68,  4,  15,    NONE,   1.0,     0.0},      /* AIOU, 0 to 15 */
    {71,  2,  3,     NONE,   1.0,     1.0},      /* BPAC area count, 1 to 4 */
    {72,  2,  3,     NONE,   1.0,     0.0},      /* BPAC area ID, 0 to 3 */
    {73,  8,  170,   NONE,   1.0,     -85.0},    /* BPAC reference latitude, ±85.00 degrees */
    {74,  9,  360,   NONE,   1.0,     -180.0},   /* BPAC reference longitude, ±180.00 degrees */
    {75,  4,  15,    NONE,   1.0,     1.0},      /* BPAC latitude node count, 1 to 16 */
    {76,  4,  15,    NONE,   1.0,     1.0},      /* BPAC longitude node count, 1 to 16 */
    {80,  12, 4094,  NONE,   0.25,    -511.75},  /* area average VTEC, ±511.75 TECU */
    {82,  7,  126,   0x7Fu,  0.25,    -15.75},   /* small VTEC residual, ±15.75 TECU */
    {83,  11, 2046,  0x7FFu, 0.25,    -255.75},  /* large VTEC residual, ±255.75 TECU */
    {99,  10, 1023,  NONE,   1.0,     0.0},      /* Galileo IODnav, 0 to 1023 */
    {100, 8,  255,   NONE,   1.0,     0.0},      /* BDS IODE/IODC, 0 to 255 */
    {101, 8,  255,   NONE,   1.0,     0.0},      /* QZSS IODE, 0 to 255 */
};

/* A field whose special-value column gives the quantity of each code. */
typedef struct {
    uint32_t field;
    uint32_t bits;
    int32_t unknown_at_zero; /* code 0 means "unknown" */
    int32_t above_at_top;    /* the last code means "greater than" the previous step */
    double values[16];       /* by code; the last code holds the bound it exceeds */
} listed_field;

/* [spec: SPARTN v2.0.3 Table 6.2, special-value column of each field] */
static const listed_field LISTED[] = {
    /* continuity indicator, seconds without a discontinuity */
    {15, 3, 0, 0, {0.0, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0, 320.0}},
    /* IODE continuity, seconds with the same IODE */
    {22, 3, 0, 0, {0.0, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0, 320.0}},
    /* user range error, m; 0 unknown, 7 "> 1.0 m" */
    {24, 3, 1, 1, {0.0, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0, 1.0}},
    /* troposphere quality, m at zenith (1 sigma), "<=" each step; 0 unknown, 7 "> 0.320 m" */
    {42, 3, 1, 1, {0.0, 0.010, 0.020, 0.040, 0.080, 0.160, 0.320, 0.320}},
    /* ionosphere quality, TECU (1 sigma), "<=" each step; 0 unknown, 15 "> 143.36 TECU" */
    {55, 4, 1, 1, {0.0, 0.03, 0.05, 0.07, 0.14, 0.28, 0.56, 1.12, 2.24, 4.48, 8.96,
                   17.92, 35.84, 71.68, 143.36, 143.36}},
    /* ionosphere shell height, km */
    {70, 2, 0, 0, {350.0, 400.0, 450.0, 500.0}},
    /* BPAC latitude and longitude node spacing, degrees */
    {77, 2, 0, 0, {2.5, 5.0, 10.0, 15.0}},
    {78, 2, 0, 0, {2.5, 5.0, 10.0, 15.0}},
};

#define COUNT(a) ((int32_t)(sizeof(a) / sizeof((a)[0])))

/* [spec: SPARTN v2.0.3 section 6.3, Table 6.2] See scale.h. The order of
 * checks is width, special value, range: 6.3 requires special values to be
 * checked before conversion, and a code above the range end that is not a
 * special value has no meaning in Table 6.2. [design: such a code is
 * reported as SSR_SCALE_OUTSIDE rather than converted, since the table
 * defines no value for it] */
int32_t ssr_scale_field(uint32_t field, uint32_t encoded, double *value)
{
    int32_t i;
    for (i = 0; i < COUNT(RANGED); ++i) {
        const ranged_field *f = &RANGED[i];
        if (f->field != field) {
            continue;
        }
        if ((encoded >> f->bits) != 0u) {
            return SSR_SCALE_OUTSIDE;
        }
        if (encoded == f->special) {
            return SSR_SCALE_SPECIAL;
        }
        if (encoded > f->top) {
            return SSR_SCALE_OUTSIDE;
        }
        *value = ssr_scale_decode(encoded, f->resolution, f->range_min);
        return SSR_SCALE_OK;
    }
    for (i = 0; i < COUNT(LISTED); ++i) {
        const listed_field *f = &LISTED[i];
        uint32_t top;
        if (f->field != field) {
            continue;
        }
        if ((encoded >> f->bits) != 0u) {
            return SSR_SCALE_OUTSIDE;
        }
        if (f->unknown_at_zero && encoded == 0u) {
            return SSR_SCALE_SPECIAL;
        }
        top = (1u << f->bits) - 1u;
        *value = f->values[encoded];
        return (f->above_at_top && encoded == top) ? SSR_SCALE_ABOVE : SSR_SCALE_OK;
    }
    return SSR_SCALE_NO_FIELD;
}

/* [spec: SPARTN v2.0.3 Table 6.2, notes of SF043 ("area offset from 2.3m"),
 * SF045 and SF048 ("offset by 0.252m"; SF048's offset added in 1.6.2 "to
 * match the offset on the small troposphere coefficient")] See scale.h. */
int32_t ssr_scale_quantity(uint32_t field, uint32_t encoded, double *value)
{
    double v;
    int32_t status = ssr_scale_field(field, encoded, &v);
    if (status != SSR_SCALE_OK) {
        if (status == SSR_SCALE_ABOVE) {
            *value = v;
        }
        return status;
    }
    if (field == 43u) {
        v = v + 2.3;
    } else if (field == 45u || field == 48u) {
        v = v + 0.252;
    }
    *value = v;
    return SSR_SCALE_OK;
}
