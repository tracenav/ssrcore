"""Host binding of the kernel ``kernels/src/windup.c``."""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/windup.h":
    int32_t SSR_WINDUP_OK
    int32_t SSR_WINDUP_DEGENERATE
    int32_t ssr_windup_frame_spartn(const double* r, const double* v_ecef,
                                    double yaw, double* x, double* y,
                                    double* z) nogil
    int32_t ssr_windup_frame_sun(const double* r, const double* r_sun,
                                 double* i, double* j, double* k) nogil
    int32_t ssr_windup_cycles(const double* r_sat, const double* r_rcv,
                              const double* a, const double* b,
                              const double* a_s, const double* b_s,
                              double* cycles) nogil
    double ssr_windup_unwrap(double previous, double current) nogil


OK = SSR_WINDUP_OK
DEGENERATE = SSR_WINDUP_DEGENERATE


cdef inline void _load(object v, double* out) except *:
    # [design: binding] A three-element sequence into a C double[3].
    out[0] = v[0]
    out[1] = v[1]
    out[2] = v[2]


def frame_spartn(r, v_ecef, double yaw):
    """[spec: SPARTN v2.0.3 section 8.14, Table 6.2 SF021] Satellite frame with yaw.

    Returns (status, x, y, z). [design: binding] A direct call of the C kernel.
    """
    cdef double cr[3]
    cdef double cv[3]
    cdef double x[3]
    cdef double y[3]
    cdef double z[3]
    _load(r, cr)
    _load(v_ecef, cv)
    status = ssr_windup_frame_spartn(cr, cv, yaw, x, y, z)
    return status, (x[0], x[1], x[2]), (y[0], y[1], y[2]), (z[0], z[1], z[2])


def frame_sun(r, r_sun):
    """[spec: TM-23/1 section 5.6.3, eqs. (5.76)-(5.78)] Nominal frame from the Sun.

    Returns (status, i, j, k). [design: binding] A direct call of the C kernel.
    """
    cdef double cr[3]
    cdef double cs[3]
    cdef double i[3]
    cdef double j[3]
    cdef double k[3]
    _load(r, cr)
    _load(r_sun, cs)
    status = ssr_windup_frame_sun(cr, cs, i, j, k)
    return status, (i[0], i[1], i[2]), (j[0], j[1], j[2]), (k[0], k[1], k[2])


def cycles(r_sat, r_rcv, a, b, a_s, b_s):
    """[spec: TM-23/1 section 5.5, eqs. (5.72), (5.74)] Fractional wind-up in cycles.

    Returns (status, cycles). [design: binding] A direct call of the C kernel.
    """
    cdef double cs[3]
    cdef double cr[3]
    cdef double ca[3]
    cdef double cb[3]
    cdef double cas[3]
    cdef double cbs[3]
    cdef double cycles = 0.0
    _load(r_sat, cs)
    _load(r_rcv, cr)
    _load(a, ca)
    _load(b, cb)
    _load(a_s, cas)
    _load(b_s, cbs)
    status = ssr_windup_cycles(cs, cr, ca, cb, cas, cbs, &cycles)
    return status, cycles


def unwrap(double previous, double current):
    """[spec: TM-23/1 section 5.5, eqs. (5.71), (5.73)] Wind-up continuous with the previous value.

    [design: binding] A direct call of the C kernel.
    """
    return ssr_windup_unwrap(previous, current)
