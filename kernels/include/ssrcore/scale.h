/* Conversion of encoded SPARTN field values.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 *
 * A field is named by the number of its SPARTN identifier: 20 for SF020.
 */
#ifndef SSRCORE_SCALE_H
#define SSRCORE_SCALE_H

#include <stdint.h>

/* Status of a field conversion. [design: DESIGN 3.1, errors are int32_t
 * status codes and values go through out-parameters; each kind of
 * non-value is told apart so that the caller can report why a term is
 * missing] */
#define SSR_SCALE_OK 0       /* *value is the field's value */
#define SSR_SCALE_SPECIAL 1  /* a special value with no number (invalid, unknown); *value untouched */
#define SSR_SCALE_OUTSIDE 2  /* encoded value outside the field's width or range; *value untouched */
#define SSR_SCALE_ABOVE 3    /* the field says "greater than *value" (the last step of a quality or URE table) */
#define SSR_SCALE_NO_FIELD 4 /* not a field this kernel converts; *value untouched */

/* [spec: SPARTN v2.0.3 section 6.3] Decoded value of an encoded field:
 * encoded * resolution + range_min. The caller checks the field's special
 * values (Table 6.2) before calling, as section 6.3 requires. */
double ssr_scale_decode(uint32_t encoded, double resolution, double range_min);

/* [spec: SPARTN v2.0.3 section 6.3, Table 6.2] Value of field SF<field> from
 * its encoded value, in the unit of Table 6.2 (m, TECU, degrees, s, or a
 * count). Ranged fields: special values are checked first (6.3), then the
 * range, then Vdecoded = Vencoded * Res + Rngmin. Fields whose special-value
 * column lists the quantity of each code (SF015, SF022 in seconds; SF024 in
 * m; SF042 in m; SF055 in TECU; SF070 in m; SF077, SF078 in degrees) are
 * read from that list. */
int32_t ssr_scale_field(uint32_t field, uint32_t encoded, double *value);

/* [spec: SPARTN v2.0.3 Table 6.2, notes of SF043, SF045, SF048] The quantity
 * a field carries: ssr_scale_field plus the offset the notes give (2.3 m
 * for the area average hydrostatic delay SF043; 0.252 m for the troposphere
 * coefficient T00, SF045 and SF048); for every other field the same as
 * ssr_scale_field. */
int32_t ssr_scale_quantity(uint32_t field, uint32_t encoded, double *value);

#endif
