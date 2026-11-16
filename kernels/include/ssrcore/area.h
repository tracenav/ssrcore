/* SPARTN geographic areas (GAD, SM 2-0): extent, inclusion, centroid,
 * polynomial offsets, grid nodes and bilinear weights.
 *
 * An area is given by the decoded fields of its GAD block (SPARTN v2.0.3
 * Table 6.24, Table 6.2 SF032 to SF037), in degrees and node counts:
 *   ref_lat      SF032, the northern-most latitude of the area
 *   ref_lon      SF033, the western-most longitude of the area
 *   lat_nodes    SF034, grid nodes north to south from ref_lat (1 to 8)
 *   lon_nodes    SF035, grid nodes west to east from ref_lon (1 to 8)
 *   lat_spacing  SF036, spacing of the nodes north to south
 *   lon_spacing  SF037, spacing of the nodes west to east
 * Positions are geodetic latitude and longitude in degrees, the unit of the
 * SPARTN area fields and of the polynomial offsets of 8.11 and 8.12
 * (DESIGN 4).
 *
 * Grid nodes are numbered as the figure of SPARTN 8.10 orders them: from the
 * area reference point (north-west corner), west to east along a row, rows
 * from north to south, so node (row, col) has index row * lon_nodes + col.
 *
 * Kernel: numbers in, numbers out; C99, explicit integer widths, double
 * precision, no allocation, no state (DESIGN 3.1, Layers).
 */
#ifndef SSRCORE_AREA_H
#define SSRCORE_AREA_H

#include <stdint.h>

/* Status codes. [design: DESIGN 3.1, errors as int32_t status codes] */
#define SSR_AREA_OK 0
#define SSR_AREA_INVALID 1 /* the definition is not an area */
#define SSR_AREA_OUTSIDE 2 /* the position is not inside the area */
#define SSR_AREA_EQUATION_TYPE 3 /* SF041 or SF054 equation type 3 to 7, "TBD" */
#define SSR_AREA_INVALID_NODE 4  /* a grid node holds the invalid value */

/* [spec: SPARTN v2.0.3 8.10 and Table 6.2 SF032 to SF037] Extent of the
 * area: out = {north, south, west, east} in degrees. north = ref_lat,
 * south = ref_lat - (lat_nodes - 1) lat_spacing, west = ref_lon,
 * east = west + (lon_nodes - 1) lon_spacing, not wrapped (it may exceed 180
 * when the area crosses the antimeridian). Returns SSR_AREA_INVALID, writing
 * nothing, for a node count outside 1 to 8, a spacing not positive, or a
 * latitude extent outside [-90, 90]. */
int32_t ssr_area_extent(double ref_lat, double ref_lon, int32_t lat_nodes,
                        int32_t lon_nodes, double lat_spacing,
                        double lon_spacing, double out[4]);

/* [spec: SPARTN v2.0.3 8.10] Whether the position may use the area's data:
 * *inside = 1 when it lies within the extent (edges included) and its
 * latitude is within +/-85 degrees, else 0. Longitude is compared modulo
 * 360 degrees. */
int32_t ssr_area_contains(double lat, double lon, double ref_lat,
                          double ref_lon, int32_t lat_nodes,
                          int32_t lon_nodes, double lat_spacing,
                          double lon_spacing, int32_t *inside);

/* [spec: SPARTN v2.0.3 8.10 figure, the polynomial expansion point
 * (centroid) at the centre of the node grid] out = {latitude, longitude} of
 * the centroid in degrees, longitude in [-180, 180). */
int32_t ssr_area_centroid(double ref_lat, double ref_lon, int32_t lat_nodes,
                          int32_t lon_nodes, double lat_spacing,
                          double lon_spacing, double out[2]);

/* [spec: SPARTN v2.0.3 8.11, 8.12] Offsets of the position from the area
 * centroid, out = {dlat, dlon} = {lat_u - lat_c, lon_u - lon_c} in degrees,
 * dlon reduced to [-180, 180). The position need not be inside the area. */
int32_t ssr_area_offsets(double lat, double lon, double ref_lat,
                         double ref_lon, int32_t lat_nodes, int32_t lon_nodes,
                         double lat_spacing, double lon_spacing,
                         double out[2]);

/* [spec: SPARTN v2.0.3 8.11, 8.12 bilinear grid interpolation, with the
 * weights of the rectangular case of 8.13] The four nodes around the
 * position and their weights, in the order W1 to W4 of 8.13:
 *   nodes[0] north-east (lambda2, phi2), weights[0] = x y
 *   nodes[1] north-west (lambda1, phi2), weights[1] = (1 - x) y
 *   nodes[2] south-west (lambda1, phi1), weights[2] = (1 - x)(1 - y)
 *   nodes[3] south-east (lambda2, phi1), weights[3] = x (1 - y)
 * with x = (lon - lambda1) / (lambda2 - lambda1) and
 * y = (lat - phi1) / (phi2 - phi1); node indices as in the header comment.
 * Returns SSR_AREA_OUTSIDE, writing nothing, when the position is outside
 * the extent (the +/-85 degree rule is ssr_area_contains's). */
int32_t ssr_area_grid(double lat, double lon, double ref_lat, double ref_lon,
                      int32_t lat_nodes, int32_t lon_nodes,
                      double lat_spacing, double lon_spacing,
                      int32_t nodes[4], double weights[4]);

/* [spec: SPARTN v2.0.3 8.11, 8.12; Table 6.2 SF041, SF054] The area
 * polynomial at offsets (dlat, dlon) from the centroid, degrees:
 *   type 0: C00
 *   type 1: C00 + C01 dlat + C10 dlon
 *   type 2: C00 + C01 dlat + C10 dlon + C11 dlat dlon
 * with c = {C00, C01, C10, C11} in the order of SF041 (troposphere, T00 to
 * T11, m) and SF054 (ionosphere, C00 to C11, TECU); the value is in the
 * unit of C00. Coefficients a type does not use are not read. Types 3 to 7
 * return SSR_AREA_EQUATION_TYPE and write nothing. [design: 8.11 and 8.12
 * give one polynomial for both atmospheres; it is evaluated in one place] */
int32_t ssr_area_polynomial(int32_t equation_type, const double c[4],
                            double dlat, double dlon, double *value);

/* [spec: SPARTN v2.0.3 8.11, 8.12 grid interpolation] The residual at the
 * position, sum over k of W_k R_k, for the four nodes and weights of
 * ssr_area_grid. valid[k] is 0 when node k holds the invalid value (SF052,
 * SF053; SF064 to SF067); then SSR_AREA_INVALID_NODE is returned and
 * nothing written [design: DESIGN rule A7, no extrapolation from three
 * nodes]. */
int32_t ssr_area_interpolate(const double weights[4], const double residuals[4],
                             const int32_t valid[4], double *value);

#endif
