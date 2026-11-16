"""Host binding of the kernel ``kernels/src/astro.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds no arithmetic.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/astro.h":
    double SSR_IERS_AU
    double SSR_TT_MINUS_GPS
    double ssr_astro_tt_centuries(int32_t week, double seconds) nogil
    double ssr_astro_earth_rotation_angle(int32_t week, double seconds, double ut1_minus_gps) nogil
    double ssr_astro_gmst(int32_t week, double seconds, double ut1_minus_gps) nogil
    double ssr_astro_mean_obliquity(double t) nogil
    void ssr_astro_nutation(double t, double *dpsi, double *deps) nogil
    void ssr_astro_moon_ecliptic(double t, double out[3]) nogil
    void ssr_astro_earth_heliocentric(double t, double out[3]) nogil
    void ssr_astro_sun_ecliptic(double t, double out[3]) nogil
    void ssr_astro_ecliptic_to_true_equator(const double lbr[3], double dpsi, double eps_a,
                                      double deps, double out[3]) nogil
    void ssr_astro_true_equator_to_ecef(const double v[3], double gst, double out[3]) nogil
    void ssr_astro_sun_moon_ecef(int32_t week, double seconds, double ut1_minus_gps,
                                 double sun[3], double moon[3]) nogil


IERS_AU = SSR_IERS_AU
TT_MINUS_GPS = SSR_TT_MINUS_GPS


def tt_centuries(int32_t week, double seconds):
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4; TM-23 3.1.1; TN36 ch. 5 sections
    5.5.7, 5.7.1] TT centuries from J2000.0."""
    return ssr_astro_tt_centuries(week, seconds)


def earth_rotation_angle(int32_t week, double seconds, double ut1_minus_gps):
    """[spec: TN36 ch. 5 eq. (5.15)] Earth Rotation Angle, rad."""
    return ssr_astro_earth_rotation_angle(week, seconds, ut1_minus_gps)


def gmst(int32_t week, double seconds, double ut1_minus_gps):
    """[spec: TN36 ch. 5 eq. (5.32)] Greenwich mean sidereal time, rad."""
    return ssr_astro_gmst(week, seconds, ut1_minus_gps)


def mean_obliquity(double t):
    """[spec: TN36 ch. 5 eq. (5.40)] Mean obliquity of date, rad."""
    return ssr_astro_mean_obliquity(t)


def nutation(double t):
    """[spec: TN36 ch. 5 eq. (5.35), Tables 5.3a and 5.3b extract] (dpsi, deps), rad."""
    cdef double dpsi, deps
    ssr_astro_nutation(t, &dpsi, &deps)
    return dpsi, deps


def moon_ecliptic(double t):
    """[spec: ELP 2000-82B sections 2 to 4, 8] Moon longitude, latitude (rad, mean ecliptic and equinox of date), distance (m)."""
    cdef double out[3]
    ssr_astro_moon_ecliptic(t, out)
    return out[0], out[1], out[2]


def earth_heliocentric(double t):
    """[spec: VSOP87 version D] Heliocentric Earth longitude, latitude (rad), radius (au)."""
    cdef double out[3]
    ssr_astro_earth_heliocentric(t, out)
    return out[0], out[1], out[2]


def sun_ecliptic(double t):
    """[spec: VSOP87 version D; TN36 ch. 1 Table 1.1] Sun longitude, latitude (rad), distance (m)."""
    cdef double out[3]
    ssr_astro_sun_ecliptic(t, out)
    return out[0], out[1], out[2]


def ecliptic_to_true_equator(lbr, double dpsi, double eps_a, double deps):
    """[spec: TN36 ch. 5 section 5.4.5] Rectangular coordinates on the true equator and equinox of date."""
    cdef double v[3]
    cdef double out[3]
    v[0], v[1], v[2] = lbr
    ssr_astro_ecliptic_to_true_equator(v, dpsi, eps_a, deps, out)
    return out[0], out[1], out[2]


def true_equator_to_ecef(v, double gst):
    """[spec: TN36 ch. 5 section 5.4.3] Rotation by Greenwich apparent sidereal time."""
    cdef double w[3]
    cdef double out[3]
    w[0], w[1], w[2] = v
    ssr_astro_true_equator_to_ecef(w, gst, out)
    return out[0], out[1], out[2]


def sun_moon_ecef(int32_t week, double seconds, double ut1_minus_gps):
    """[spec: ELP 2000-82B; VSOP87; TN36 ch. 5] ((Sun x, y, z), (Moon x, y, z)) in ECEF, m."""
    cdef double sun[3]
    cdef double moon[3]
    ssr_astro_sun_moon_ecef(week, seconds, ut1_minus_gps, sun, moon)
    return (sun[0], sun[1], sun[2]), (moon[0], moon[1], moon[2])
