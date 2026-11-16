"""Host binding of the kernel ``kernels/src/area.c``.

Each function returns the kernel's status code first, then its outputs; the
toolkit (``ssrcore.toolkit.area``) turns status codes into exceptions.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/area.h":
    int32_t SSR_AREA_OK
    int32_t SSR_AREA_INVALID
    int32_t SSR_AREA_OUTSIDE
    int32_t SSR_AREA_EQUATION_TYPE
    int32_t SSR_AREA_INVALID_NODE
    int32_t ssr_area_polynomial(int32_t equation_type, const double c[4],
                                double dlat, double dlon, double *value) nogil
    int32_t ssr_area_interpolate(const double weights[4], const double residuals[4],
                                 const int32_t valid[4], double *value) nogil
    int32_t ssr_area_extent(double ref_lat, double ref_lon, int32_t lat_nodes,
                            int32_t lon_nodes, double lat_spacing,
                            double lon_spacing, double out[4]) nogil
    int32_t ssr_area_contains(double lat, double lon, double ref_lat,
                              double ref_lon, int32_t lat_nodes,
                              int32_t lon_nodes, double lat_spacing,
                              double lon_spacing, int32_t *inside) nogil
    int32_t ssr_area_centroid(double ref_lat, double ref_lon, int32_t lat_nodes,
                              int32_t lon_nodes, double lat_spacing,
                              double lon_spacing, double out[2]) nogil
    int32_t ssr_area_offsets(double lat, double lon, double ref_lat,
                             double ref_lon, int32_t lat_nodes, int32_t lon_nodes,
                             double lat_spacing, double lon_spacing,
                             double out[2]) nogil
    int32_t ssr_area_grid(double lat, double lon, double ref_lat, double ref_lon,
                          int32_t lat_nodes, int32_t lon_nodes,
                          double lat_spacing, double lon_spacing,
                          int32_t nodes[4], double weights[4]) nogil


OK = SSR_AREA_OK
INVALID = SSR_AREA_INVALID
OUTSIDE = SSR_AREA_OUTSIDE
EQUATION_TYPE = SSR_AREA_EQUATION_TYPE
INVALID_NODE = SSR_AREA_INVALID_NODE


def extent(double ref_lat, double ref_lon, int32_t lat_nodes, int32_t lon_nodes,
           double lat_spacing, double lon_spacing):
    """[spec: SPARTN v2.0.3 8.10, Table 6.2 SF032 to SF037] (status, north, south, west, east).

    [design: binding] A direct call of the C kernel.
    """
    cdef double out[4]
    out[:] = [0.0, 0.0, 0.0, 0.0]
    cdef int32_t s = ssr_area_extent(ref_lat, ref_lon, lat_nodes, lon_nodes,
                                     lat_spacing, lon_spacing, out)
    return s, out[0], out[1], out[2], out[3]


def contains(double lat, double lon, double ref_lat, double ref_lon,
             int32_t lat_nodes, int32_t lon_nodes, double lat_spacing,
             double lon_spacing):
    """[spec: SPARTN v2.0.3 8.10] (status, inside), inside 1 or 0.

    [design: binding] A direct call of the C kernel.
    """
    cdef int32_t inside = 0
    cdef int32_t s = ssr_area_contains(lat, lon, ref_lat, ref_lon, lat_nodes,
                                       lon_nodes, lat_spacing, lon_spacing, &inside)
    return s, inside


def centroid(double ref_lat, double ref_lon, int32_t lat_nodes, int32_t lon_nodes,
             double lat_spacing, double lon_spacing):
    """[spec: SPARTN v2.0.3 8.10 figure] (status, latitude, longitude) of the centroid.

    [design: binding] A direct call of the C kernel.
    """
    cdef double out[2]
    out[:] = [0.0, 0.0]
    cdef int32_t s = ssr_area_centroid(ref_lat, ref_lon, lat_nodes, lon_nodes,
                                       lat_spacing, lon_spacing, out)
    return s, out[0], out[1]


def offsets(double lat, double lon, double ref_lat, double ref_lon,
            int32_t lat_nodes, int32_t lon_nodes, double lat_spacing,
            double lon_spacing):
    """[spec: SPARTN v2.0.3 8.11, 8.12] (status, dlat, dlon) from the centroid.

    [design: binding] A direct call of the C kernel.
    """
    cdef double out[2]
    out[:] = [0.0, 0.0]
    cdef int32_t s = ssr_area_offsets(lat, lon, ref_lat, ref_lon, lat_nodes,
                                      lon_nodes, lat_spacing, lon_spacing, out)
    return s, out[0], out[1]


def grid(double lat, double lon, double ref_lat, double ref_lon,
         int32_t lat_nodes, int32_t lon_nodes, double lat_spacing,
         double lon_spacing):
    """[spec: SPARTN v2.0.3 8.11, 8.12; 8.13 weights W1 to W4] (status, nodes, weights).

    [design: binding] A direct call of the C kernel.
    """
    cdef int32_t nodes[4]
    cdef double weights[4]
    nodes[:] = [0, 0, 0, 0]
    weights[:] = [0.0, 0.0, 0.0, 0.0]
    cdef int32_t s = ssr_area_grid(lat, lon, ref_lat, ref_lon, lat_nodes,
                                   lon_nodes, lat_spacing, lon_spacing, nodes, weights)
    return (s, (nodes[0], nodes[1], nodes[2], nodes[3]),
            (weights[0], weights[1], weights[2], weights[3]))


def polynomial(int32_t equation_type, c, double dlat, double dlon):
    """[spec: SPARTN v2.0.3 8.11, 8.12; Table 6.2 SF041, SF054] (status, value).

    ``c`` is (C00, C01, C10, C11). [design: binding] A direct call of the C kernel.
    """
    cdef double cc[4]
    cdef double value = 0.0
    cc[0], cc[1], cc[2], cc[3] = c
    cdef int32_t s = ssr_area_polynomial(equation_type, cc, dlat, dlon, &value)
    return s, value


def interpolate(weights, residuals, valid):
    """[spec: SPARTN v2.0.3 8.11, 8.12] (status, sum of W_k R_k).

    [design: binding] A direct call of the C kernel.
    """
    cdef double w[4]
    cdef double r[4]
    cdef int32_t v[4]
    cdef double value = 0.0
    w[0], w[1], w[2], w[3] = weights
    r[0], r[1], r[2], r[3] = residuals
    v[0], v[1], v[2], v[3] = valid
    cdef int32_t s = ssr_area_interpolate(w, r, v, &value)
    return s, value
