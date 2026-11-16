"""Host binding of the kernel ``kernels/src/troposphere.c``.

[design: binding] Each function is a direct call of the C kernel, so that
host and microcontroller run the same arithmetic; it adds none. A status
other than 0 is returned to the caller as it is.
"""

from libc.stdint cimport int32_t


cdef extern from "ssrcore/troposphere.h":
    int32_t SSR_TROPOSPHERE_OK
    int32_t SSR_TROPOSPHERE_ELEVATION
    int32_t SSR_TROPOSPHERE_HEIGHT
    double ssr_troposphere_marini(double sin_e, double a, double b, double c) nogil
    void ssr_troposphere_niell_coefficients(double latitude, double doy, double hydrostatic[3], double wet[3]) nogil
    int32_t ssr_troposphere_niell_mapping(double latitude, double height, double doy, double elevation,
                                          double *hydrostatic, double *wet) nogil
    double ssr_troposphere_slant(double hydrostatic_zenith, double polynomial_residual,
                                 double grid_residual, double hydrostatic_mapping,
                                 double wet_mapping) nogil
    void ssr_troposphere_unb3m_temperature(double latitude, double doy, double *temperature,
                                           double *lapse_rate) nogil
    int32_t ssr_troposphere_unb3m_height_factor(double latitude, double temperature,
                                                double lapse_rate, double height,
                                                double *factor) nogil


OK = SSR_TROPOSPHERE_OK
ELEVATION = SSR_TROPOSPHERE_ELEVATION
HEIGHT = SSR_TROPOSPHERE_HEIGHT


def marini(double sin_e, double a, double b, double c):
    """[spec: TM-23 5.4.2.2.1 eq. (5.69)] m(E, a, b, c) for sin E."""
    return ssr_troposphere_marini(sin_e, a, b, c)


def niell_coefficients(double latitude, double doy):
    """[spec: TM-23 5.4.2.2.1 eq. (5.70), Tables 5.2, 5.3] ((a_d, b_d, c_d), (a_w, b_w, c_w))."""
    cdef double h[3]
    cdef double w[3]
    ssr_troposphere_niell_coefficients(latitude, doy, h, w)
    return (h[0], h[1], h[2]), (w[0], w[1], w[2])


def niell_mapping(double latitude, double height, double doy, double elevation):
    """[spec: TM-23 5.4.2.2.1 eqs. (5.67) to (5.70); SPARTN v2.0.3 8.12] (status, M_dry, M_wet)."""
    cdef double mh = 0.0
    cdef double mw = 0.0
    cdef int32_t status = ssr_troposphere_niell_mapping(latitude, height, doy, elevation, &mh, &mw)
    return status, mh, mw


def slant(double hydrostatic_zenith, double polynomial_residual, double grid_residual,
          double hydrostatic_mapping, double wet_mapping):
    """[spec: SPARTN v2.0.3 8.12] T_h m_h + (T_P + T_G) m_w."""
    return ssr_troposphere_slant(hydrostatic_zenith, polynomial_residual, grid_residual,
                                 hydrostatic_mapping, wet_mapping)


def unb3m_temperature(double latitude, double doy):
    """[spec: UNB3m eqs. (1) to (3), Table 1; TM-23 5.4.2.1 eq. (5.62)] (T0 K, beta K/m)."""
    cdef double t = 0.0
    cdef double beta = 0.0
    ssr_troposphere_unb3m_temperature(latitude, doy, &t, &beta)
    return t, beta


def unb3m_height_factor(double latitude, double temperature, double lapse_rate,
                        double height):
    """[spec: UNB3m eqs. (4), (6); design: DESIGN Q7] (status, d_zh(h) / d_zh(0))."""
    cdef double factor = 0.0
    cdef int32_t status = ssr_troposphere_unb3m_height_factor(latitude, temperature, lapse_rate,
                                                              height, &factor)
    return status, factor
