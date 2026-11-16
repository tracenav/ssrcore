"""Host binding of the kernel ``kernels/src/correction.c``."""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/correction.h":
    double SSR_SPEED_OF_LIGHT
    int32_t SSR_CORRECTION_OK
    int32_t SSR_CORRECTION_DEGENERATE
    int32_t ssr_correction_orbit_frame(const double *r, const double *r_dot,
                                       double *e_radial, double *e_along, double *e_cross) nogil
    int32_t ssr_correction_orbit(const double *r, const double *r_dot,
                                 const double *d_o, double *d_x) nogil
    int32_t ssr_correction_position(const double *r, const double *r_dot,
                                    const double *d_o, double *x) nogil
    double ssr_correction_clock(double t_b, double d_c) nogil


SPEED_OF_LIGHT = SSR_SPEED_OF_LIGHT
OK = SSR_CORRECTION_OK
DEGENERATE = SSR_CORRECTION_DEGENERATE


cdef inline void _load(v, double out[3]) except *:
    out[0] = v[0]
    out[1] = v[1]
    out[2] = v[2]


def orbit_frame(r, r_dot):
    """[spec: SPARTN v2.0.3 section 8.8] (status, e_radial, e_along, e_cross).

    [design: binding] Vectors are 3-tuples; None when the status is not OK.
    """
    cdef double a[3]
    cdef double b[3]
    cdef double er[3]
    cdef double ea[3]
    cdef double ec[3]
    _load(r, a)
    _load(r_dot, b)
    cdef int32_t status = ssr_correction_orbit_frame(a, b, er, ea, ec)
    if status != SSR_CORRECTION_OK:
        return status, None, None, None
    return status, (er[0], er[1], er[2]), (ea[0], ea[1], ea[2]), (ec[0], ec[1], ec[2])


def orbit(r, r_dot, d_o):
    """[spec: SPARTN v2.0.3 section 8.8] (status, dX) for dO = (radial, along, cross).

    [design: binding] dX is a 3-tuple; None when the status is not OK.
    """
    cdef double a[3]
    cdef double b[3]
    cdef double o[3]
    cdef double d[3]
    _load(r, a)
    _load(r_dot, b)
    _load(d_o, o)
    cdef int32_t status = ssr_correction_orbit(a, b, o, d)
    if status != SSR_CORRECTION_OK:
        return status, None
    return status, (d[0], d[1], d[2])


def position(r, r_dot, d_o):
    """[spec: SPARTN v2.0.3 section 8.8] (status, X = X_b - dX).

    [design: binding] X is a 3-tuple; None when the status is not OK.
    """
    cdef double a[3]
    cdef double b[3]
    cdef double o[3]
    cdef double x[3]
    _load(r, a)
    _load(r_dot, b)
    _load(d_o, o)
    cdef int32_t status = ssr_correction_position(a, b, o, x)
    if status != SSR_CORRECTION_OK:
        return status, None
    return status, (x[0], x[1], x[2])


def clock(double t_b, double d_c):
    """[spec: SPARTN v2.0.3 section 8.7] t_s = t_b - dC / c, s.

    [design: binding] A direct call of the C kernel.
    """
    return ssr_correction_clock(t_b, d_c)


cdef extern from "ssrcore/correction.h":
    void ssr_correction_orbit_at(const double *d_o, const double *rate, double dt, double *out) nogil
    double ssr_correction_clock_at(double c0, double c1, double c2, double dt) nogil
    double ssr_correction_clock_rate_at(double c1, double c2, double dt) nogil


def orbit_at(d_o, rate, double dt):
    """[spec: RTCM 10403.3 3.5.13.5 eq. (3.12-6); IGS SSR v1.00 4.1] dO + dO_dot * dt, m.

    [design: binding] Vectors are 3-tuples.
    """
    cdef double o[3]
    cdef double r[3]
    cdef double out[3]
    _load(d_o, o)
    _load(rate, r)
    ssr_correction_orbit_at(o, r, dt, out)
    return (out[0], out[1], out[2])


def clock_at(double c0, double c1, double c2, double dt):
    """[spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8); IGS SSR v1.00 4.2] C0 + C1 dt + C2 dt^2, m.

    [design: binding] A direct call of the C kernel.
    """
    return ssr_correction_clock_at(c0, c1, c2, dt)


def clock_rate_at(double c1, double c2, double dt):
    """[spec: RTCM 10403.3 3.5.13.6 eq. (3.12-8), its derivative] C1 + 2 C2 dt, m/s.

    [design: binding] A direct call of the C kernel.
    """
    return ssr_correction_clock_rate_at(c1, c2, dt)
