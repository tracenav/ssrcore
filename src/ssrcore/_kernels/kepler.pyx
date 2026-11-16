"""Host binding of the kernel ``kernels/src/kepler.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds no arithmetic.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/kepler.h":
    int32_t SSR_KEPLER_OK
    int32_t SSR_KEPLER_BAD_ELEMENTS
    int32_t SSR_KEPLER_NOT_CONVERGED

    ctypedef struct ssr_kepler_orbit:
        int32_t system
        int32_t geo
        int32_t toe_week
        double toe
        double sqrt_a
        double e
        double m0
        double delta_n
        double omega
        double omega0
        double omega_dot
        double i0
        double idot
        double cuc
        double cus
        double crc
        double crs
        double cic
        double cis

    double ssr_kepler_semicircles(double semicircles) nogil
    int32_t ssr_kepler_beidou_geo(int32_t prn) nogil
    int32_t ssr_kepler_eccentric_anomaly(double m, double e, double *ecc_anomaly) nogil
    int32_t ssr_kepler_state(const ssr_kepler_orbit *orbit, int32_t week, double sow,
                             double position[3], double velocity[3], double *ecc_anomaly) nogil
    double ssr_kepler_clock(double af0, double af1, double af2,
                            int32_t toc_week, double toc, int32_t week, double sow) nogil
    double ssr_kepler_clock_drift(double af1, double af2,
                                  int32_t toc_week, double toc, int32_t week, double sow) nogil


OK = SSR_KEPLER_OK
BAD_ELEMENTS = SSR_KEPLER_BAD_ELEMENTS
NOT_CONVERGED = SSR_KEPLER_NOT_CONVERGED


def semicircles(double value):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2] Radians from semicircles, pi =
    3.1415926535898."""
    return ssr_kepler_semicircles(value)


def beidou_geo(int32_t prn):
    """[spec: BDS-SIS-ICD-B1I-3.0 4.3 Table 4-1] Whether the BeiDou SV ID is a GEO satellite."""
    return ssr_kepler_beidou_geo(prn) != 0


def eccentric_anomaly(double m, double e):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV] (status, E) for mean
    anomaly m, eccentricity e."""
    cdef double ek = 0.0
    cdef int32_t status = ssr_kepler_eccentric_anomaly(m, e, &ek)
    return status, ek


def state(int32_t system, bint geo, int32_t toe_week, double toe, double sqrt_a, double e,
          double m0, double delta_n, double omega, double omega0, double omega_dot,
          double i0, double idot, double cuc, double cus, double crc, double crs,
          double cic, double cis, int32_t week, double sow):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IV; Galileo OS SIS ICD 2.2 5.1.1;
    BDS-SIS-ICD-B1I-3.0 5.2.4.12]

    (status, position, velocity, E) at (week, sow); angles in radians.
    """
    cdef ssr_kepler_orbit k
    k.system = system
    k.geo = 1 if geo else 0
    k.toe_week = toe_week
    k.toe = toe
    k.sqrt_a = sqrt_a
    k.e = e
    k.m0 = m0
    k.delta_n = delta_n
    k.omega = omega
    k.omega0 = omega0
    k.omega_dot = omega_dot
    k.i0 = i0
    k.idot = idot
    k.cuc = cuc
    k.cus = cus
    k.crc = crc
    k.crs = crs
    k.cic = cic
    k.cis = cis
    cdef double p[3]
    cdef double v[3]
    cdef double ek = 0.0
    p[:] = [0.0, 0.0, 0.0]
    v[:] = [0.0, 0.0, 0.0]
    cdef int32_t status = ssr_kepler_state(&k, week, sow, p, v, &ek)
    return status, (p[0], p[1], p[2]), (v[0], v[1], v[2]), ek


def clock(double af0, double af1, double af2, int32_t toc_week, double toc,
          int32_t week, double sow):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4;
    BDS-SIS-ICD-B1I-3.0 5.2.4.9]

    Clock offset (s) of the broadcast polynomial, without the relativistic term.
    """
    return ssr_kepler_clock(af0, af1, af2, toc_week, toc, week, sow)


def clock_drift(double af1, double af2, int32_t toc_week, double toc, int32_t week, double sow):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] Clock drift (s/s) of the
    broadcast polynomial."""
    return ssr_kepler_clock_drift(af1, af2, toc_week, toc, week, sow)
