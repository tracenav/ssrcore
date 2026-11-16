"""SPARTN geographic areas: extent, inclusion, grid nodes and weights, polynomial and grid interpolation.

DESIGN.md section 2.7. An area is the decoded content of a GAD area block
(SPARTN v2.0.3 Table 6.24; Table 6.2 SF032 to SF037), angles in degrees.
Positions are geodetic latitude and longitude in degrees, the unit of the
SPARTN area fields and polynomial offsets (DESIGN section 4).

Grid nodes are numbered as the figure of SPARTN 8.10 orders them: from the
area reference point (north-west corner), west to east along a row, rows
from north to south; node (row, col) is ``row * longitude_nodes + col``.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

from ssrcore._kernels import area as _k
from ssrcore.toolkit import KernelError


class Area(NamedTuple):
    """[spec: SPARTN v2.0.3 Table 6.24, Table 6.2 SF032 to SF037] One GAD area, decoded.

    ``reference_latitude`` (SF032) is the northern-most latitude and
    ``reference_longitude`` (SF033) the western-most longitude, degrees;
    ``latitude_nodes`` (SF034, 1 to 8) count north to south and
    ``longitude_nodes`` (SF035, 1 to 8) west to east; the spacings (SF036,
    SF037) are in degrees.
    """

    reference_latitude: float
    reference_longitude: float
    latitude_nodes: int
    longitude_nodes: int
    latitude_spacing: float
    longitude_spacing: float


def _checked(status: int) -> None:
    """[design: a kernel status becomes an exception at the toolkit, the public edge]"""
    if status == _k.INVALID:
        raise KernelError("area", status, "not a SPARTN area: node counts 1 to 8, positive spacings, "
                          "latitudes within ±90°")
    if status == _k.OUTSIDE:
        raise KernelError("area", status, "position outside the area")
    if status != _k.OK:
        raise KernelError("area", status, "area polynomial equation type undefined (SF041, SF054)")


def extent(area: Area) -> tuple[float, float, float, float]:
    """[spec: SPARTN v2.0.3 8.10, Table 6.2 SF032 to SF037] (north, south, west, east), degrees.

    ``east`` is ``west`` plus the longitude span, not wrapped: it exceeds 180
    for an area across the antimeridian.
    """
    s, north, south, west, east = _k.extent(*area)
    _checked(s)
    return north, south, west, east


def contains(area: Area, latitude: float, longitude: float) -> bool:
    """[spec: SPARTN v2.0.3 8.10] Whether a position may use the area's atmosphere data.

    True when the position is inside the area, edges included, and its
    latitude is within ±85°. Longitudes are compared modulo 360°.
    """
    s, inside = _k.contains(latitude, longitude, *area)
    _checked(s)
    return bool(inside)


def centroid(area: Area) -> tuple[float, float]:
    """[spec: SPARTN v2.0.3 8.10 figure] The polynomial expansion point (latitude, longitude), degrees.

    The centre of the node grid; longitude in [-180, 180).
    """
    s, latitude, longitude = _k.centroid(*area)
    _checked(s)
    return latitude, longitude


def offsets(area: Area, latitude: float, longitude: float) -> tuple[float, float]:
    """[spec: SPARTN v2.0.3 8.11, 8.12] (δφ, δλ) of a position from the area centroid, degrees.

    δλ is reduced to [-180, 180).
    """
    s, dlat, dlon = _k.offsets(latitude, longitude, *area)
    _checked(s)
    return dlat, dlon


def grid(area: Area, latitude: float, longitude: float) -> tuple[tuple[int, int, int, int],
                                                                 tuple[float, float, float, float]]:
    """[spec: SPARTN v2.0.3 8.11, 8.12; weights of 8.13, rectangular case] Nodes and bilinear weights.

    Returns the four nodes of the grid cell around the position and their
    weights, in the order W1 to W4 of 8.13: north-east, north-west,
    south-west, south-east. The weights sum to one. Raises
    :class:`~ssrcore.toolkit.KernelError` when the position is outside the area.
    """
    s, nodes, weights = _k.grid(latitude, longitude, *area)
    _checked(s)
    return nodes, weights


def polynomial(equation_type: int, coefficients: Sequence[float], dlat: float, dlon: float) -> float:
    """[spec: SPARTN v2.0.3 8.11, 8.12; Table 6.2 SF041, SF054] The area polynomial at (δφ, δλ).

    ``coefficients`` are (C00, C01, C10, C11) in the order of SF041
    (troposphere, m) and SF054 (ionosphere, TECU), as many as the equation
    type uses (1, 3 or 4): type 0 is C00, type 1 adds C01 δφ + C10 δλ, type
    2 adds C11 δφ δλ; ``dlat``, ``dlon`` are the offsets of :func:`offsets`,
    degrees. The value is in the unit of C00. Raises
    :class:`~ssrcore.toolkit.KernelError` for types 3 to 7, which SPARTN
    leaves undefined, and ``ValueError`` for too few coefficients.
    """
    needed = {0: 1, 1: 3, 2: 4}.get(equation_type, 0)
    if len(coefficients) < needed:
        raise ValueError(f"equation type {equation_type} needs {needed} coefficients")
    c = (tuple(coefficients) + (0.0, 0.0, 0.0))[:4]
    s, value = _k.polynomial(equation_type, c, dlat, dlon)
    _checked(s)
    return value


def interpolate(weights: Sequence[float], residuals: Sequence[float | None]) -> float | None:
    """[spec: SPARTN v2.0.3 8.11, 8.12] The grid residual Σ W_k R_k at a position.

    ``weights`` and ``residuals`` are for the four nodes :func:`grid`
    returns, in its order (troposphere SF052, SF053 in m; ionosphere SF064
    to SF067 in TECU). A residual of ``None`` is a node holding the invalid
    value; then the result is ``None`` and the polynomial is used alone
    [design: DESIGN rule A7, no extrapolation from three nodes].
    """
    if len(weights) != 4 or len(residuals) != 4:
        raise ValueError("four weights and four residuals")
    valid = tuple(0 if r is None else 1 for r in residuals)
    values = tuple(0.0 if r is None else r for r in residuals)
    s, value = _k.interpolate(tuple(weights), values, valid)
    return None if s == _k.INVALID_NODE else value
