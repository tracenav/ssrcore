"""Host binding of the kernel ``kernels/src/scale.c``."""

from libc.stdint cimport int32_t, uint32_t


cdef extern from "ssrcore/scale.h":
    int32_t SSR_SCALE_OK
    int32_t SSR_SCALE_SPECIAL
    int32_t SSR_SCALE_OUTSIDE
    int32_t SSR_SCALE_ABOVE
    int32_t SSR_SCALE_NO_FIELD
    double ssr_scale_decode(uint32_t encoded, double resolution, double range_min) nogil
    int32_t ssr_scale_field(uint32_t field, uint32_t encoded, double *value) nogil
    int32_t ssr_scale_quantity(uint32_t field, uint32_t encoded, double *value) nogil


OK = SSR_SCALE_OK
SPECIAL = SSR_SCALE_SPECIAL
OUTSIDE = SSR_SCALE_OUTSIDE
ABOVE = SSR_SCALE_ABOVE
NO_FIELD = SSR_SCALE_NO_FIELD


def decode(uint32_t encoded, double resolution, double range_min):
    """[spec: SPARTN v2.0.3 section 6.3] Decoded value of an encoded field.

    [design: binding] A direct call of the C kernel, so that host and
    microcontroller run the same arithmetic.
    """
    return ssr_scale_decode(encoded, resolution, range_min)


def field(uint32_t field, uint32_t encoded):
    """[spec: SPARTN v2.0.3 section 6.3, Table 6.2] (status, value) of field SF<field>.

    [design: binding] The value is None unless the status is OK or ABOVE.
    """
    cdef double value = 0.0
    cdef int32_t status = ssr_scale_field(field, encoded, &value)
    return status, (value if status == SSR_SCALE_OK or status == SSR_SCALE_ABOVE else None)


def quantity(uint32_t field, uint32_t encoded):
    """[spec: SPARTN v2.0.3 Table 6.2, notes of SF043, SF045, SF048] (status, value) with the notes' offsets.

    [design: binding] The value is None unless the status is OK or ABOVE.
    """
    cdef double value = 0.0
    cdef int32_t status = ssr_scale_quantity(field, encoded, &value)
    return status, (value if status == SSR_SCALE_OK or status == SSR_SCALE_ABOVE else None)
