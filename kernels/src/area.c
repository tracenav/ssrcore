/* SPARTN geographic areas. See ssrcore/area.h. */

#include <math.h>

#include "ssrcore/area.h"

/* [design: the full circle in degrees, the unit of the SPARTN area fields] */
#define FULL_CIRCLE 360.0
/* [spec: SPARTN v2.0.3 8.10] "should not be used above +/- 85 degrees
 * latitude" */
#define LATITUDE_LIMIT 85.0
/* [spec: SPARTN v2.0.3 Table 6.2 SF034, SF035] node counts 1 to 8 */
#define NODES_MAX 8

/* [design: angle x reduced to [0, 360) degrees, so that longitudes compare
 * across the antimeridian; SPARTN 8.13 asks that the +/-180 degree
 * discontinuity be considered] */
static double east_of(double x)
{
    return x - FULL_CIRCLE * floor(x / FULL_CIRCLE);
}

/* [design: angle x reduced to [-180, 180) degrees, the range of SF033] */
static double centred(double x)
{
    return x - FULL_CIRCLE * floor((x + 0.5 * FULL_CIRCLE) / FULL_CIRCLE);
}

/* [spec: SPARTN v2.0.3 Table 6.2 SF032 to SF037] Whether the fields define
 * an area: node counts 1 to 8, spacings positive, the latitude extent on
 * the sphere. */
static int32_t valid(double ref_lat, int32_t lat_nodes, int32_t lon_nodes,
                     double lat_spacing, double lon_spacing)
{
    if (lat_nodes < 1 || lat_nodes > NODES_MAX)
        return 0;
    if (lon_nodes < 1 || lon_nodes > NODES_MAX)
        return 0;
    if (!(lat_spacing > 0.0) || !(lon_spacing > 0.0))
        return 0;
    if (!(ref_lat <= 90.0)
        || !(ref_lat - (double)(lat_nodes - 1) * lat_spacing >= -90.0))
        return 0;
    return 1;
}

/* [spec: SPARTN v2.0.3 8.10 and Table 6.2 SF032 to SF037] */
int32_t ssr_area_extent(double ref_lat, double ref_lon, int32_t lat_nodes,
                        int32_t lon_nodes, double lat_spacing,
                        double lon_spacing, double out[4])
{
    if (!valid(ref_lat, lat_nodes, lon_nodes, lat_spacing, lon_spacing))
        return SSR_AREA_INVALID;
    out[0] = ref_lat;
    out[1] = ref_lat - (double)(lat_nodes - 1) * lat_spacing;
    out[2] = ref_lon;
    out[3] = ref_lon + (double)(lon_nodes - 1) * lon_spacing;
    return SSR_AREA_OK;
}

/* [spec: SPARTN v2.0.3 8.10] Inside the extent, edges included
 * [design: the 8.10 figure marks positions strictly inside or outside;
 * an edge is included because the bilinear weights of 8.13 are defined
 * there and adjacent areas then leave no gap], and not above +/-85 degrees
 * latitude. */
int32_t ssr_area_contains(double lat, double lon, double ref_lat,
                          double ref_lon, int32_t lat_nodes,
                          int32_t lon_nodes, double lat_spacing,
                          double lon_spacing, int32_t *inside)
{
    double e[4];
    double x;
    if (ssr_area_extent(ref_lat, ref_lon, lat_nodes, lon_nodes, lat_spacing,
                        lon_spacing, e) != SSR_AREA_OK)
        return SSR_AREA_INVALID;
    x = east_of(lon - e[2]);
    *inside = lat <= e[0] && lat >= e[1] && x <= e[3] - e[2]
              && fabs(lat) <= LATITUDE_LIMIT;
    return SSR_AREA_OK;
}

/* [spec: SPARTN v2.0.3 8.10 figure] The polynomial expansion point
 * (centroid) is drawn at the centre of the node grid: half the latitude
 * extent south of the reference point and half the longitude extent east of
 * it. */
int32_t ssr_area_centroid(double ref_lat, double ref_lon, int32_t lat_nodes,
                          int32_t lon_nodes, double lat_spacing,
                          double lon_spacing, double out[2])
{
    if (!valid(ref_lat, lat_nodes, lon_nodes, lat_spacing, lon_spacing))
        return SSR_AREA_INVALID;
    out[0] = ref_lat - 0.5 * (double)(lat_nodes - 1) * lat_spacing;
    out[1] = centred(ref_lon + 0.5 * (double)(lon_nodes - 1) * lon_spacing);
    return SSR_AREA_OK;
}

/* [spec: SPARTN v2.0.3 8.11, 8.12] dphi = phi_u - phi_c,
 * dlambda = lambda_u - lambda_c [design: dlambda reduced to [-180, 180) so
 * that an area across the antimeridian has small offsets, as 8.13 asks the
 * discontinuity be considered] */
int32_t ssr_area_offsets(double lat, double lon, double ref_lat,
                         double ref_lon, int32_t lat_nodes, int32_t lon_nodes,
                         double lat_spacing, double lon_spacing,
                         double out[2])
{
    double c[2];
    if (ssr_area_centroid(ref_lat, ref_lon, lat_nodes, lon_nodes, lat_spacing,
                          lon_spacing, c) != SSR_AREA_OK)
        return SSR_AREA_INVALID;
    out[0] = lat - c[0];
    out[1] = centred(lon - c[1]);
    return SSR_AREA_OK;
}

/* [spec: SPARTN v2.0.3 8.11, 8.12 "bi-linear grid interpolation"; weights
 * W1 to W4 of 8.13, rectangular case, with delta lambda_pp read as
 * lambda - lambda1 and delta phi_pp as phi - phi1, the reading under which
 * each weight is 1 at its own node] The cell is the one whose north-west
 * node is (floor of the row, floor of the column) of the position, kept
 * inside the grid so that a position on the southern or eastern edge uses
 * the last cell [design: the cell is then always made of grid nodes]. A
 * direction with one node has no cell: its two nodes of the cell are that
 * node, with x = 0 (one column) or y = 1 (one row) [design: the weights
 * then still sum to one and give the node's value on the line it defines].
 */
int32_t ssr_area_grid(double lat, double lon, double ref_lat, double ref_lon,
                      int32_t lat_nodes, int32_t lon_nodes,
                      double lat_spacing, double lon_spacing,
                      int32_t nodes[4], double weights[4])
{
    double e[4];
    double east, row, col, phi2, lambda1, x, y;
    int32_t i0, j0, i1, j1;

    if (ssr_area_extent(ref_lat, ref_lon, lat_nodes, lon_nodes, lat_spacing,
                        lon_spacing, e) != SSR_AREA_OK)
        return SSR_AREA_INVALID;
    east = east_of(lon - e[2]);
    if (!(lat <= e[0] && lat >= e[1] && east <= e[3] - e[2]))
        return SSR_AREA_OUTSIDE;

    row = (e[0] - lat) / lat_spacing;
    col = east / lon_spacing;
    i0 = (int32_t)floor(row);
    j0 = (int32_t)floor(col);
    if (i0 > lat_nodes - 2)
        i0 = lat_nodes - 2;
    if (i0 < 0)
        i0 = 0;
    if (j0 > lon_nodes - 2)
        j0 = lon_nodes - 2;
    if (j0 < 0)
        j0 = 0;

    if (lat_nodes > 1) {
        i1 = i0 + 1;
        phi2 = e[0] - (double)i0 * lat_spacing;
        y = (lat - (phi2 - lat_spacing)) / lat_spacing;
    } else {
        i1 = i0;
        y = 1.0;
    }
    if (lon_nodes > 1) {
        j1 = j0 + 1;
        lambda1 = (double)j0 * lon_spacing;
        x = (east - lambda1) / lon_spacing;
    } else {
        j1 = j0;
        x = 0.0;
    }

    nodes[0] = i0 * lon_nodes + j1; /* north-east (lambda2, phi2) */
    nodes[1] = i0 * lon_nodes + j0; /* north-west (lambda1, phi2) */
    nodes[2] = i1 * lon_nodes + j0; /* south-west (lambda1, phi1) */
    nodes[3] = i1 * lon_nodes + j1; /* south-east (lambda2, phi1) */
    weights[0] = x * y;
    weights[1] = (1.0 - x) * y;
    weights[2] = (1.0 - x) * (1.0 - y);
    weights[3] = x * (1.0 - y);
    return SSR_AREA_OK;
}

/* [spec: SPARTN v2.0.3 8.11, 8.12; Table 6.2 SF041, SF054] See ssrcore/area.h. */
int32_t ssr_area_polynomial(int32_t equation_type, const double c[4],
                            double dlat, double dlon, double *value)
{
    switch (equation_type) {
    case 0:
        *value = c[0];
        return SSR_AREA_OK;
    case 1:
        *value = c[0] + c[1] * dlat + c[2] * dlon;
        return SSR_AREA_OK;
    case 2:
        *value = c[0] + c[1] * dlat + c[2] * dlon + c[3] * dlat * dlon;
        return SSR_AREA_OK;
    default:
        return SSR_AREA_EQUATION_TYPE;
    }
}

/* [spec: SPARTN v2.0.3 8.11, 8.12] sum_{k=1}^{4} W_k R_k
 * [design: DESIGN rule A7 for invalid nodes] */
int32_t ssr_area_interpolate(const double weights[4], const double residuals[4],
                             const int32_t valid[4], double *value)
{
    int32_t k;
    double sum = 0.0;
    for (k = 0; k < 4; k++) {
        if (!valid[k]) {
            return SSR_AREA_INVALID_NODE;
        }
        sum += weights[k] * residuals[k];
    }
    *value = sum;
    return SSR_AREA_OK;
}
