"""Host binding of the kernel ``kernels/src/tides.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds no arithmetic.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/tides.h":
    double SSR_IERS_GM_SUN
    double SSR_IERS_GM_EARTH
    double SSR_IERS_MOON_EARTH_MASS_RATIO
    double SSR_IERS_SUN_EARTH_MASS_RATIO
    double SSR_IERS_EARTH_RADIUS
    void ssr_tides_degree2(const double site[3], const double body[3], double mass_ratio,
                          double h0, double h2, double l0, double l2, double out[3]) nogil
    void ssr_tides_degree3(const double site[3], const double body[3], double mass_ratio,
                          double out[3]) nogil
    void ssr_tides_l1(const double site[3], const double body[3], double mass_ratio,
                     double l1_diurnal, double l1_semidiurnal, double out[3]) nogil
    void ssr_tides_out_of_phase(const double site[3], const double body[3], double mass_ratio,
                               double hi_diurnal, double li_diurnal,
                               double hi_semidiurnal, double li_semidiurnal, double out[3]) nogil
    void ssr_tides_step2(const double site[3], double t, double gmst, double out[3]) nogil
    void ssr_tides_displacement(const double site[3], const double sun[3], const double moon[3],
                                double t, double gmst, double out[3]) nogil
    void ssr_tides_permanent(const double site[3], double out[3]) nogil


cdef void _vec(v, double out[3]):
    # [design: binding] a Python 3-sequence into a C array; no arithmetic.
    out[0], out[1], out[2] = v


IERS_GM_SUN = SSR_IERS_GM_SUN
IERS_GM_EARTH = SSR_IERS_GM_EARTH
SUN_MASS_RATIO = SSR_IERS_SUN_EARTH_MASS_RATIO
MOON_MASS_RATIO = SSR_IERS_MOON_EARTH_MASS_RATIO
IERS_EARTH_RADIUS = SSR_IERS_EARTH_RADIUS


def degree2(site, body, double mass_ratio, double h0, double h2, double l0, double l2):
    """[spec: TN36 7.1.1 eqs. (5), (2)] Degree 2 in-phase displacement, m."""
    cdef double s[3]
    cdef double b[3]
    cdef double out[3]
    _vec(site, s)
    _vec(body, b)
    ssr_tides_degree2(s, b, mass_ratio, h0, h2, l0, l2, out)
    return out[0], out[1], out[2]


def degree3(site, body, double mass_ratio):
    """[spec: TN36 7.1.1 eq. (6)] Degree 3 displacement, m."""
    cdef double s[3]
    cdef double b[3]
    cdef double out[3]
    _vec(site, s)
    _vec(body, b)
    ssr_tides_degree3(s, b, mass_ratio, out)
    return out[0], out[1], out[2]


def l1(site, body, double mass_ratio, double l1_diurnal, double l1_semidiurnal):
    """[spec: TN36 7.1.1 eqs. (7) to (9)] Transverse displacement from l(1), m."""
    cdef double s[3]
    cdef double b[3]
    cdef double out[3]
    _vec(site, s)
    _vec(body, b)
    ssr_tides_l1(s, b, mass_ratio, l1_diurnal, l1_semidiurnal, out)
    return out[0], out[1], out[2]


def out_of_phase(site, body, double mass_ratio, double hi_diurnal, double li_diurnal,
                 double hi_semidiurnal, double li_semidiurnal):
    """[spec: TN36 7.1.1 eqs. (10), (11)] Out-of-phase displacement, m."""
    cdef double s[3]
    cdef double b[3]
    cdef double out[3]
    _vec(site, s)
    _vec(body, b)
    ssr_tides_out_of_phase(s, b, mass_ratio, hi_diurnal, li_diurnal, hi_semidiurnal,
                          li_semidiurnal, out)
    return out[0], out[1], out[2]


def step2(site, double t, double gmst):
    """[spec: TN36 7.1.1 eqs. (12), (13), Tables 7.3a, 7.3b] Step 2 corrections, m."""
    cdef double s[3]
    cdef double out[3]
    _vec(site, s)
    ssr_tides_step2(s, t, gmst, out)
    return out[0], out[1], out[2]


def displacement(site, sun, moon, double t, double gmst):
    """[spec: TN36 7.1.1 steps 1 and 2] Displacement of the site, ECEF m."""
    cdef double s[3]
    cdef double su[3]
    cdef double mo[3]
    cdef double out[3]
    _vec(site, s)
    _vec(sun, su)
    _vec(moon, mo)
    ssr_tides_displacement(s, su, mo, t, gmst, out)
    return out[0], out[1], out[2]


def permanent(site):
    """[spec: TN36 7.1.1.2 eqs. (14a), (14b)] Permanent deformation, ECEF m."""
    cdef double s[3]
    cdef double out[3]
    _vec(site, s)
    ssr_tides_permanent(s, out)
    return out[0], out[1], out[2]
