"""Host binding of the kernel ``kernels/src/geodesy.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; status codes are returned
as they are, and turned into exceptions by :mod:`ssrcore.toolkit.geodesy`.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/geodesy.h":
    double SSR_WGS84_A
    double SSR_WGS84_F
    double SSR_WGS84_GM
    double SSR_WGS84_OMEGA_E
    double SSR_PI
    int32_t SSR_GEODESY_OK
    int32_t SSR_GEODESY_NO_SOLUTION
    int32_t SSR_GEODESY_ZERO_RANGE
    double ssr_geodesy_eccentricity_squared(double f) nogil
    void ssr_geodesy_geodetic_to_ecef(double lat, double lon, double h, double a,
                                      double f, double ecef[3]) nogil
    int32_t ssr_geodesy_ecef_to_geodetic(const double ecef[3], double a, double f,
                                         double geodetic[3]) nogil
    void ssr_geodesy_ecef_to_enu(double lat, double lon, const double d[3],
                                 double enu[3]) nogil
    void ssr_geodesy_enu_to_ecef(double lat, double lon, const double enu[3],
                                 double d[3]) nogil
    int32_t ssr_geodesy_line_of_sight(const double satellite[3],
                                      const double receiver[3], double lat,
                                      double lon, double unit[3],
                                      double *elevation, double *azimuth) nogil
    int32_t ssr_geodesy_range(const double satellite[3], const double velocity[3],
                              const double receiver[3], double *range, double *rate) nogil
    void ssr_geodesy_rotate_z(const double v[3], double angle, double out[3]) nogil


WGS84_A = SSR_WGS84_A
WGS84_F = SSR_WGS84_F
WGS84_GM = SSR_WGS84_GM
WGS84_OMEGA_E = SSR_WGS84_OMEGA_E
PI = SSR_PI
OK = SSR_GEODESY_OK
NO_SOLUTION = SSR_GEODESY_NO_SOLUTION
ZERO_RANGE = SSR_GEODESY_ZERO_RANGE


def eccentricity_squared(double f):
    """[spec: TM-23 B.3] e^2 = 2f - f^2."""
    return ssr_geodesy_eccentricity_squared(f)


def geodetic_to_ecef(double lat, double lon, double h, double a, double f):
    """[spec: TM-23 B.1.1] ECEF (x, y, z) of (lat, lon, h) on (a, f)."""
    cdef double out[3]
    ssr_geodesy_geodetic_to_ecef(lat, lon, h, a, f, out)
    return out[0], out[1], out[2]


def ecef_to_geodetic(double x, double y, double z, double a, double f):
    """[spec: TM-23 B.1.2] (status, (lat, lon, h)) of an ECEF position."""
    cdef double r[3]
    cdef double out[3]
    r[0] = x
    r[1] = y
    r[2] = z
    out[0] = out[1] = out[2] = 0.0
    status = ssr_geodesy_ecef_to_geodetic(r, a, f, out)
    return status, (out[0], out[1], out[2])


def ecef_to_enu(double lat, double lon, double dx, double dy, double dz):
    """[spec: TM-23 B.2.2] (east, north, up) of an ECEF vector at (lat, lon)."""
    cdef double d[3]
    cdef double out[3]
    d[0] = dx
    d[1] = dy
    d[2] = dz
    ssr_geodesy_ecef_to_enu(lat, lon, d, out)
    return out[0], out[1], out[2]


def enu_to_ecef(double lat, double lon, double e, double n, double u):
    """[spec: TM-23 B.2.1] ECEF (dx, dy, dz) of an ENU vector at (lat, lon)."""
    cdef double v[3]
    cdef double out[3]
    v[0] = e
    v[1] = n
    v[2] = u
    ssr_geodesy_enu_to_ecef(lat, lon, v, out)
    return out[0], out[1], out[2]


def line_of_sight(satellite, receiver, double lat, double lon):
    """[spec: TM-23 B.3] (status, unit, elevation, azimuth) from receiver to satellite."""
    cdef double rcv[3]
    cdef double sat[3]
    cdef double unit[3]
    cdef double el = 0.0
    cdef double az = 0.0
    rcv[0], rcv[1], rcv[2] = receiver
    sat[0], sat[1], sat[2] = satellite
    unit[0] = unit[1] = unit[2] = 0.0
    status = ssr_geodesy_line_of_sight(sat, rcv, lat, lon, unit, &el, &az)
    return status, (unit[0], unit[1], unit[2]), el, az


def range_rate(satellite, velocity, receiver):
    """[spec: TM-23 B.3, eq. B.13] (status, range, rate) from receiver to satellite.

    [design: binding] A direct call of the C kernel; adds no arithmetic.
    """
    cdef double sat[3]
    cdef double vel[3]
    cdef double rcv[3]
    cdef double rho = float("nan")
    cdef double rate = float("nan")
    sat[0], sat[1], sat[2] = satellite
    vel[0], vel[1], vel[2] = velocity
    rcv[0], rcv[1], rcv[2] = receiver
    status = ssr_geodesy_range(sat, vel, rcv, &rho, &rate)
    return status, rho, rate


def rotate_z(v, double angle):
    """[spec: Ashby 2003 section 2] v in a frame turned by ``angle`` about z.

    [design: binding] A direct call of the C kernel; adds no arithmetic.
    """
    cdef double a[3]
    cdef double out[3]
    a[0], a[1], a[2] = v
    ssr_geodesy_rotate_z(a, angle, out)
    return out[0], out[1], out[2]
