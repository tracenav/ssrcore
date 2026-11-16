"""Host binding of the kernel ``kernels/src/unb3m.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds none. A status
other than 0 is returned to the caller as it is.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/unb3m.h":
    int32_t SSR_UNB3M_OK
    int32_t SSR_UNB3M_HEIGHT
    int32_t SSR_UNB3M_PARAMETERS
    void ssr_unb3m_parameters(double latitude, double doy, double *pressure, double *humidity,
                              double *lambda_) nogil
    double ssr_unb3m_vapour_pressure(double humidity, double temperature, double pressure) nogil
    int32_t ssr_unb3m_zenith(double latitude, double height, double pressure, double temperature,
                             double vapour_pressure, double lapse_rate, double lambda_,
                             double *hydrostatic, double *wet) nogil


OK = SSR_UNB3M_OK
HEIGHT = SSR_UNB3M_HEIGHT
PARAMETERS = SSR_UNB3M_PARAMETERS


def parameters(double latitude, double doy):
    """[spec: UNB3m eqs. (1) to (3), Tables 1 and 2; TM-23 5.4.2.1 eq. (5.62)] (P0 mbar, RH %, lambda)."""
    cdef double p = 0.0
    cdef double rh = 0.0
    cdef double lam = 0.0
    ssr_unb3m_parameters(latitude, doy, &p, &rh, &lam)
    return p, rh, lam


def vapour_pressure(double humidity, double temperature, double pressure):
    """[spec: UNB3m eqs. (9) to (11)] e0 in mbar."""
    return ssr_unb3m_vapour_pressure(humidity, temperature, pressure)


def zenith(double latitude, double height, double pressure, double temperature,
           double vapour_pressure, double lapse_rate, double lambda_):
    """[spec: UNB3m eqs. (4) to (7)] (status, d_zh m, d_znh m)."""
    cdef double dh = 0.0
    cdef double dw = 0.0
    cdef int32_t status = ssr_unb3m_zenith(latitude, height, pressure, temperature, vapour_pressure,
                                           lapse_rate, lambda_, &dh, &dw)
    return status, dh, dw
