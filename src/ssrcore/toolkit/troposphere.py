"""Troposphere: Niell mapping functions, SPARTN slant delay, height reduction.

DESIGN.md section 2.7. Typed wrappers of ``ssrcore._kernels.troposphere``;
no state, no clock. Angles in radians, heights in metres (ellipsoidal),
``doy`` the day of year counted from January 0.0 (1.0 at 00:00 of 1
January), as TM-23 5.4.2.2.1 eq. (5.70) counts it.

The slant delay of a satellite at the user's position is, after SPARTN
v2.0.3 8.12::

    m_h, m_w = niell_mapping(lat, h, doy, elevation)
    t_p = area.polynomial(sf041, (t00, t01, t10, t11), dlat, dlon)
    t_g = area.interpolate(weights, residuals) or 0.0  # 0.0 without a grid, or by A7
    slant = slant(t_h, t_p, t_g, m_h, m_w)

where ``t_h`` is the area's hydrostatic zenith delay SF043 as decoded (with
its 2.3 m offset), optionally multiplied by
``hydrostatic_height_factor(lat, h, doy)`` (DESIGN Q7), and ``area`` is
``ssrcore.toolkit.area`` (offsets in degrees from the centroid, DESIGN 4).
"""

from typing import NamedTuple

from ssrcore._kernels import troposphere as _k
from ssrcore.toolkit import KernelError


class NiellCoefficients(NamedTuple):
    """Niell coefficients at one latitude and day. [spec: TM-23 5.4.2.2.1 eq. (5.70)]"""

    hydrostatic: tuple[float, float, float]
    wet: tuple[float, float, float]


class NiellMapping(NamedTuple):
    """Hydrostatic and wet mapping values. [spec: TM-23 5.4.2.2.1 eqs. (5.67), (5.68)]"""

    hydrostatic: float
    wet: float


def marini(sin_elevation: float, a: float, b: float, c: float) -> float:
    """[spec: TM-23 5.4.2.2.1 eq. (5.69)] Marini continued fraction normalised to 1 at zenith."""
    return _k.marini(sin_elevation, a, b, c)


def niell_coefficients(latitude: float, doy: float) -> NiellCoefficients:
    """[spec: TM-23 5.4.2.2.1 eq. (5.70), Tables 5.2, 5.3] Hydrostatic (a, b, c) and wet (a, b, c).

    [design: southern hemisphere, half a year added to the phase; latitudes
    beyond 15 and 75 degrees take the end values of the tables]
    """
    h, w = _k.niell_coefficients(latitude, doy)
    return NiellCoefficients(h, w)


def niell_mapping(latitude: float, height: float, doy: float, elevation: float) -> NiellMapping:
    """[spec: TM-23 5.4.2.2.1 eqs. (5.67) to (5.70); SPARTN v2.0.3 8.12] Niell m_h(E, h) and m_w(E).

    ``height`` is ellipsoidal, in metres (SPARTN 8.12). Raises KernelError
    for an elevation outside (0, pi/2].
    """
    status, mh, mw = _k.niell_mapping(latitude, height, doy, elevation)
    if status != _k.OK:
        raise KernelError("troposphere", status, f"elevation {elevation!r} rad is not in (0, pi/2]")
    return NiellMapping(mh, mw)


def slant(hydrostatic_zenith: float, polynomial_residual: float, grid_residual: float,
          hydrostatic_mapping: float, wet_mapping: float) -> float:
    """[spec: SPARTN v2.0.3 8.12] Slant delay T_h m_h(E, h) + (T_P + T_G) m_w(E), m."""
    return _k.slant(hydrostatic_zenith, polynomial_residual, grid_residual,
                    hydrostatic_mapping, wet_mapping)


def unb3m_temperature(latitude: float, doy: float) -> tuple[float, float]:
    """[spec: UNB3m eqs. (1) to (3), Table 1; TM-23 5.4.2.1 eq. (5.62)] (T0 in K, beta in K/m)."""
    return _k.unb3m_temperature(latitude, doy)


def unb3m_height_factor(latitude: float, temperature: float, lapse_rate: float,
                        height: float) -> float:
    """[spec: UNB3m eqs. (4), (6); design: DESIGN Q7] d_zh(h) / d_zh(0) for given T0 and beta.

    Raises KernelError where 1 - beta h / T0 <= 0.
    """
    status, factor = _k.unb3m_height_factor(latitude, temperature, lapse_rate, height)
    if status != _k.OK:
        raise KernelError("troposphere", status, f"height {height!r} m is beyond the lapse-rate profile")
    return factor


def hydrostatic_height_factor(latitude: float, height: float, doy: float) -> float:
    """[spec: UNB3m eqs. (1) to (4), (6); design: DESIGN Q7] Factor bringing a hydrostatic
    zenith delay referred to 0 m ellipsoidal height (SPARTN SF043) to ``height``.
    """
    temperature, lapse_rate = unb3m_temperature(latitude, doy)
    return unb3m_height_factor(latitude, temperature, lapse_rate, height)
