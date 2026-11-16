"""Host binding of the kernel ``kernels/src/glonass.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds no arithmetic.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/glonass.h":
    int32_t SSR_GLONASS_OK
    int32_t SSR_GLONASS_TOO_FAR
    double SSR_GLONASS_MAX_SPAN
    double SSR_GLONASS_MAX_STEP

    void ssr_glonass_derivative(const double state[6], const double lunisolar[3],
                                double rate[6]) nogil
    void ssr_glonass_rk4_step(const double state[6], const double lunisolar[3], double h,
                              double out[6]) nogil
    int32_t ssr_glonass_state(const double state_tb[6], const double lunisolar[3],
                              int32_t tb_week, double tb, int32_t week, double sow,
                              double out[6]) nogil
    double ssr_glonass_clock(double tau_n, double gamma_n,
                             int32_t tb_week, double tb, int32_t week, double sow) nogil


OK = SSR_GLONASS_OK
TOO_FAR = SSR_GLONASS_TOO_FAR
MAX_SPAN = SSR_GLONASS_MAX_SPAN
MAX_STEP = SSR_GLONASS_MAX_STEP


cdef void _load(double out[6], position, velocity):
    """[design: binding] Copy two Python 3-sequences into a C state."""
    out[0], out[1], out[2] = position
    out[3], out[4], out[5] = velocity


cdef void _load3(double out[3], vector):
    """[design: binding] Copy a Python 3-sequence into a C vector."""
    out[0], out[1], out[2] = vector


def derivative(position, velocity, lunisolar):
    """[spec: GLONASS ICD 5.1 A.3.1.2 with RTCM 10403.3 3.5.13] (velocity, acceleration)."""
    cdef double s[6]
    cdef double a[3]
    cdef double d[6]
    _load(s, position, velocity)
    _load3(a, lunisolar)
    ssr_glonass_derivative(s, a, d)
    return (d[0], d[1], d[2]), (d[3], d[4], d[5])


def rk4_step(position, velocity, lunisolar, double h):
    """[spec: GLONASS ICD 5.1 A.3.1.1, A.3.1.2] (position, velocity) after one RK4 step of h s."""
    cdef double s[6]
    cdef double a[3]
    cdef double o[6]
    _load(s, position, velocity)
    _load3(a, lunisolar)
    ssr_glonass_rk4_step(s, a, h, o)
    return (o[0], o[1], o[2]), (o[3], o[4], o[5])


def state(position, velocity, lunisolar, int32_t tb_week, double tb, int32_t week, double sow):
    """[spec: GLONASS ICD 5.1 A.3.1.2; RTCM 10403.3 3.5.13] (status, position, velocity) at (week, sow)."""
    cdef double s[6]
    cdef double a[3]
    cdef double o[6]
    _load(s, position, velocity)
    _load3(a, lunisolar)
    o[:] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    cdef int32_t status = ssr_glonass_state(s, a, tb_week, tb, week, sow, o)
    return status, (o[0], o[1], o[2]), (o[3], o[4], o[5])


def clock(double tau_n, double gamma_n, int32_t tb_week, double tb, int32_t week, double sow):
    """[spec: GLONASS ICD 5.1 3.3.3, 4.4] Satellite clock offset (s), -tau_n + gamma_n (t - tb)."""
    return ssr_glonass_clock(tau_n, gamma_n, tb_week, tb, week, sow)
