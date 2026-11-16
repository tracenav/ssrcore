"""Host binding of the kernel ``kernels/src/ionosphere.c``.

Each function returns the kernel's status code first, then its outputs; the
toolkit (``ssrcore.toolkit.ionosphere``) turns status codes into exceptions.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/ionosphere.h":
    int32_t SSR_IONOSPHERE_OK
    int32_t SSR_IONOSPHERE_BAD_FREQUENCY
    int32_t SSR_IONOSPHERE_BAD_POSITION
    int32_t SSR_IONOSPHERE_BELOW_LAYER
    int32_t SSR_IONOSPHERE_BAD_DEGREE
    int32_t ssr_ionosphere_delay(double tec_polynomial, double tec_grid,
                                 double frequency, double *code, double *phase) nogil
    int32_t ssr_ionosphere_spherical(const double receiver[3], const double unit[3],
                                     double out[5]) nogil
    int32_t ssr_ionosphere_pierce(double latitude, double longitude, double height,
                                  double elevation, double azimuth,
                                  double layer_height, double out[4]) nogil
    int32_t ssr_ionosphere_vtec(int32_t degree, int32_t order, const double *cosine,
                                const double *sine, double latitude, double longitude,
                                double t, double *vtec) nogil


OK = SSR_IONOSPHERE_OK
BAD_FREQUENCY = SSR_IONOSPHERE_BAD_FREQUENCY
BAD_POSITION = SSR_IONOSPHERE_BAD_POSITION
BELOW_LAYER = SSR_IONOSPHERE_BELOW_LAYER
BAD_DEGREE = SSR_IONOSPHERE_BAD_DEGREE

# [spec: IGS SSR v1.00 Table 25, N and M at most 16: (N + 1)(N + 2)/2 cosine
# coefficients at most]
cdef enum:
    _MAX_COEFFICIENTS = 153


def delay(double tec_polynomial, double tec_grid, double frequency):
    """[spec: SPARTN v2.0.3 8.11] (status, code delay m, phase delay m).

    [design: binding] A direct call of the C kernel.
    """
    cdef double code = 0.0
    cdef double phase = 0.0
    cdef int32_t s = ssr_ionosphere_delay(tec_polynomial, tec_grid, frequency,
                                          &code, &phase)
    return s, code, phase


def spherical(receiver, unit):
    """[spec: IGS SSR v1.00 4.5.1] (status, latitude, longitude, height, elevation, azimuth)
    of a receiver and a unit line of sight in the spherical model.

    [design: binding] A direct call of the C kernel.
    """
    cdef double r[3]
    cdef double u[3]
    cdef double out[5]
    r[:] = [receiver[0], receiver[1], receiver[2]]
    u[:] = [unit[0], unit[1], unit[2]]
    out[:] = [0.0, 0.0, 0.0, 0.0, 0.0]
    cdef int32_t s = ssr_ionosphere_spherical(r, u, out)
    return s, out[0], out[1], out[2], out[3], out[4]


def pierce(double latitude, double longitude, double height, double elevation,
           double azimuth, double layer_height):
    """[spec: IGS SSR v1.00 4.5.1] (status, latitude, longitude, psi, slant factor) of a
    layer's pierce point.

    [design: binding] A direct call of the C kernel.
    """
    cdef double out[4]
    out[:] = [0.0, 0.0, 0.0, 0.0]
    cdef int32_t s = ssr_ionosphere_pierce(latitude, longitude, height, elevation,
                                           azimuth, layer_height, out)
    return s, out[0], out[1], out[2], out[3]


def vtec(int32_t degree, int32_t order, cosine, sine, double latitude,
         double longitude, double t):
    """[spec: IGS SSR v1.00 4.5.1, 8.4.1] (status, VTEC TECU) of one layer.

    [design: binding] A direct call of the C kernel; the coefficient lists are
    copied into arrays of the largest size a message can carry, and a list too
    long for them is refused with BAD_DEGREE.
    """
    cdef double c[_MAX_COEFFICIENTS]
    cdef double sn[_MAX_COEFFICIENTS]
    cdef double value = 0.0
    cdef Py_ssize_t k
    if len(cosine) > _MAX_COEFFICIENTS or len(sine) > _MAX_COEFFICIENTS:
        return SSR_IONOSPHERE_BAD_DEGREE, 0.0
    for k in range(_MAX_COEFFICIENTS):
        c[k] = cosine[k] if k < len(cosine) else 0.0
        sn[k] = sine[k] if k < len(sine) else 0.0
    cdef int32_t s = ssr_ionosphere_vtec(degree, order, c, sn, latitude, longitude,
                                         t, &value)
    return s, value
