"""Ionosphere: the delay of slant TEC on a signal.

DESIGN.md section 2.7; SPARTN v2.0.3 section 8.11. TEC in TEC units (TECU,
10^16 electrons/m^2). The slant TEC of an area is ``toolkit.area``'s::

    i_p = area.polynomial(sf054, (c00, c01, c10, c11), *area.offsets(a, lat, lon))
    nodes, weights = area.grid(a, lat, lon)
    i_g = area.interpolate(weights, [residuals[n] for n in nodes])  # None: A7
    code, phase = delay(i_p, frequency, i_g or 0.0)
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

from ssrcore._kernels import ionosphere as _k
from ssrcore.toolkit import KernelError

Vector = tuple[float, float, float]

# [design: DESIGN 3.1, the toolkit names each status of the kernel]
_MESSAGES = {
    _k.BAD_FREQUENCY: "frequency must be positive (Hz)",
    _k.BAD_POSITION: "no receiver position or no line of sight",
    _k.BELOW_LAYER: "the line of sight does not cross the layer outward",
    _k.BAD_DEGREE: "degree and order are 0 to 16 (IGS SSR Table 25)",
}


def _check(status: int) -> None:
    """[design: DESIGN 3.1, a kernel status other than OK is a KernelError]"""
    if status != _k.OK:
        raise KernelError("ionosphere", status, _MESSAGES.get(status, "kernel error"))


class SphericalLook(NamedTuple):
    """The receiver and the line of sight in the 6370 km sphere of IGS SSR 4.5.1:
    geocentric latitude, longitude (rad), height above the sphere (m), elevation
    and azimuth (rad). [spec: IGS SSR v1.00 4.5.1]"""

    latitude: float
    longitude: float
    height: float
    elevation: float
    azimuth: float


class PiercePoint(NamedTuple):
    """A layer's pierce point: latitude, longitude, central angle ψ (rad) and the
    slant factor 1/sin(E + ψ). [spec: IGS SSR v1.00 4.5.1]"""

    latitude: float
    longitude: float
    psi: float
    slant_factor: float


def spherical(receiver: Vector, unit: Vector) -> SphericalLook:
    """[spec: IGS SSR v1.00 4.5.1] the receiver (ECEF, m) and the unit line of sight to
    the satellite (ECEF, frame of reception) in the spherical Earth model. Raises
    KernelError for a receiver at the centre or a zero line of sight."""
    s, *out = _k.spherical(receiver, unit)
    _check(s)
    return SphericalLook(*out)


def pierce(look: SphericalLook, layer_height: float) -> PiercePoint:
    """[spec: IGS SSR v1.00 4.5.1] the pierce point of a layer ``layer_height`` metres
    above the sphere. Raises KernelError when the line of sight does not cross it."""
    s, *out = _k.pierce(look.latitude, look.longitude, look.height, look.elevation, look.azimuth,
                        layer_height)
    _check(s)
    return PiercePoint(*out)


def vtec(degree: int, order: int, cosine: Sequence[float], sine: Sequence[float], latitude: float,
         longitude: float, seconds_of_day: float) -> float:
    """[spec: IGS SSR v1.00 4.5.1, 8.4.1] one layer's VTEC (TECU) at a pierce point:
    ``cosine`` C(n,m) for m = 0..M, n = m..N and ``sine`` S(n,m) for m = 1..M, n = m..N,
    in TECU, the order of the message; ``seconds_of_day`` the t of the sun-fixed
    longitude. Negative values are returned as they are (the caller replaces them
    by 0, 4.5.1). Raises KernelError for a degree or order outside 0 to 16."""
    m = min(order, degree)
    need_c = (m + 1) * (degree + 1) - m * (m + 1) // 2
    need_s = need_c - (degree + 1)
    if len(cosine) != need_c or len(sine) != need_s:
        raise ValueError(f"degree {degree}, order {order} take {need_c} cosine and {need_s} sine "
                         f"coefficients, not {len(cosine)} and {len(sine)}")
    s, value = _k.vtec(degree, order, cosine, sine, latitude, longitude, seconds_of_day)
    _check(s)
    return value


def delay(tec_polynomial: float, frequency: float, tec_grid: float = 0.0) -> tuple[float, float]:
    """[spec: SPARTN v2.0.3 8.11] (code, phase) slant ionosphere delay of a signal, metres.

    code = +40.3·10¹⁶ / f² · (I_P + I_G) and phase = −code, for frequency
    ``f`` in Hz and TEC in TECU; ``tec_grid`` is zero when no grid is
    present (8.11). Raises :class:`~ssrcore.toolkit.KernelError` for a
    frequency that is not positive.
    """
    s, code, phase = _k.delay(tec_polynomial, tec_grid, frequency)
    _check(s)
    return code, phase
