"""Host binding of the kernel ``kernels/src/relativity.c``."""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/relativity.h":
    int32_t SSR_RELATIVITY_OK
    int32_t SSR_RELATIVITY_THROUGH_CENTRE
    double ssr_relativity_clock(const double position[3], const double velocity[3]) nogil
    double ssr_relativity_clock_kepler(double mu, double e, double sqrt_a, double eccentric_anomaly) nogil
    int32_t ssr_relativity_shapiro(const double satellite[3], const double receiver[3], double *delay) nogil
    int32_t SSR_RELATIVITY_ZERO_RANGE
    int32_t ssr_relativity_sagnac(const double satellite[3], const double receiver[3], double *delta) nogil


OK = SSR_RELATIVITY_OK
THROUGH_CENTRE = SSR_RELATIVITY_THROUGH_CENTRE
ZERO_RANGE = SSR_RELATIVITY_ZERO_RANGE


def clock(double px, double py, double pz, double vx, double vy, double vz):
    """[spec: SPARTN v2.0.3 section 8.7; IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] -2
    r.v / c^2, s.

    [design: binding] A direct call of the C kernel; adds no arithmetic.
    """
    cdef double p[3]
    cdef double v[3]
    p[0], p[1], p[2] = px, py, pz
    v[0], v[1], v[2] = vx, vy, vz
    return ssr_relativity_clock(p, v)


def clock_kepler(double mu, double e, double sqrt_a, double eccentric_anomaly):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4]
    F e sqrt(A) sin E, s.

    [design: binding] A direct call of the C kernel; adds no arithmetic.
    """
    return ssr_relativity_clock_kepler(mu, e, sqrt_a, eccentric_anomaly)


def shapiro(double sx, double sy, double sz, double rx, double ry, double rz):
    """[spec: Ashby 2003 section 10, Eq. (96)] (status, Shapiro delay in m).

    [design: binding] A direct call of the C kernel; adds no arithmetic. The
    delay is NaN unless the status is ``OK``.
    """
    cdef double s[3]
    cdef double r[3]
    cdef double delay = float("nan")
    s[0], s[1], s[2] = sx, sy, sz
    r[0], r[1], r[2] = rx, ry, rz
    cdef int32_t status = ssr_relativity_shapiro(s, r, &delay)
    return status, delay


def sagnac(double sx, double sy, double sz, double rx, double ry, double rz):
    """[spec: Ashby 2003 section 2] (status, Sagnac term in m), exact in omega_E tau.

    [design: binding] A direct call of the C kernel; adds no arithmetic. The
    term is NaN unless the status is ``OK``.
    """
    cdef double s[3]
    cdef double r[3]
    cdef double delta = float("nan")
    s[0], s[1], s[2] = sx, sy, sz
    r[0], r[1], r[2] = rx, ry, rz
    cdef int32_t status = ssr_relativity_sagnac(s, r, &delta)
    return status, delta
