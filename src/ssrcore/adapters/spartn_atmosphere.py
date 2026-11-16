"""gnsswire.spartn HPAC (SM 1-0 to 1-4) and GAD (SM 2-0) messages to atmosphere inputs.

DESIGN.md section 2.2. A decoded ``gnsswire.spartn.Message`` holds the raw
unsigned integers of SPARTN v2.0.3 section 6.5 (Tables 6.13 to 6.24); the
functions here turn them into ``HpacInput`` and ``GadInput`` with every value
in the unit Table 6.2 gives it, through ``ssrcore.toolkit.ssr.field`` and
``quantity`` (section 6.3 with the special values and note offsets of Table
6.2). No bit parsing: that is gnsswire's.

The time tag and the stream (TF008 to TF011) are carried as sent; resolving
the tag to GPS time is the engine's (rules T1 to T3).

A value that the specification does not define (a field outside its range,
an equation type of 3 to 7) raises ``ValueError``: the message is then a
frame that failed to decode (DESIGN 3.5), not a partial input.
"""

from __future__ import annotations

from gnsswire import spartn

from ssrcore.model.types import (
    AreaDefinition,
    GadInput,
    HpacArea,
    HpacInput,
    IonosphereBlock,
    IonosphereSatellite,
    Sat,
    SpartnTag,
    StreamId,
    System,
    TroposphereBlock,
)
from ssrcore.adapters.spartn_ocb import SATELLITE_MASK_BITS, satellites  # noqa: F401  (one definition)
from ssrcore.model.atmosphere import area_of
from ssrcore.toolkit import area as _area
from ssrcore.toolkit import ssr as _ssr

# [spec: SPARTN v2.0.3 Table 6.2 SF041, SF054] coefficients each defined
# equation type uses, in the order (00, 01, 10, 11); types 3 to 7 are TBD.
_COEFFICIENT_COUNT = {0: 1, 1: 3, 2: 4}

# [spec: SPARTN v2.0.3 Tables 6.17, 6.18, 6.21, 6.22] the field of each
# coefficient (00, 01, 10, 11) of the small and large blocks.
_TROPO_FIELDS = {0: (45, 46, 46, 47), 1: (48, 49, 49, 50)}       # SF044
_IONO_FIELDS = {0: (57, 58, 58, 59), 1: (60, 61, 61, 62)}        # SF056

# [spec: SPARTN v2.0.3 Table 6.16 SF051; Table 6.20 SF063] the residual field
# each residual field size indicator selects.
_TROPO_RESIDUAL_FIELD = {0: 52, 1: 53}
_IONO_RESIDUAL_FIELD = {0: 64, 1: 65, 2: 66, 3: 67}


def _value(sf: int, encoded: int, quantity: bool = False) -> float:
    """[spec: SPARTN v2.0.3 6.3 and Table 6.2; design: a value outside the field's range or a
    special value where none is allowed is undefined, so the message fails (DESIGN 3.5)]"""
    v = (_ssr.quantity if quantity else _ssr.field)(sf, encoded)
    if v.status != _ssr.Status.OK:
        raise ValueError(f"SF{sf:03d} = {encoded} has no value ({v.status.name})")
    return v.value


def _residual(sf: int, encoded: int) -> float | None:
    """[spec: SPARTN v2.0.3 Table 6.2 SF052, SF053, SF064 to SF067: the all-ones code is invalid]
    A grid residual, None for the invalid value (rule A7)."""
    v = _ssr.field(sf, encoded)
    if v.status == _ssr.Status.SPECIAL:
        return None
    if v.status != _ssr.Status.OK:
        raise ValueError(f"SF{sf:03d} = {encoded} has no value ({v.status.name})")
    return v.value


def _coefficients(equation_type: int, fields: tuple[int, ...], raws: tuple[int | None, ...],
                  t00_offset: bool) -> tuple[float, ...]:
    """[spec: SPARTN v2.0.3 Tables 6.17, 6.18, 6.21, 6.22; Table 6.2 SF041, SF054, and the notes of
    SF045, SF048 (T00 offset by 0.252 m)] the coefficients of the equation type, in order."""
    count = _COEFFICIENT_COUNT.get(equation_type)
    if count is None:
        raise ValueError(f"equation type {equation_type} is undefined (rule A8)")
    if any(r is None for r in raws[:count]):
        raise ValueError(f"equation type {equation_type} needs {count} coefficients")
    return tuple(_value(fields[k], raws[k], quantity=t00_offset and k == 0) for k in range(count))


def _troposphere(t: spartn.Troposphere) -> TroposphereBlock:
    """[spec: SPARTN v2.0.3 Tables 6.16 to 6.18; Table 6.2 SF041 to SF053] a troposphere block."""
    block = t.small if t.coefficient_size == 0 else t.large
    coefficients = _coefficients(t.equation_type, _TROPO_FIELDS[t.coefficient_size],
                                 (block.t00, block.t01, block.t10, block.t11), True)
    residuals = None
    if t.residuals is not None:
        sf = _TROPO_RESIDUAL_FIELD[t.residual_size]
        residuals = tuple(_residual(sf, r) for r in t.residuals)
    return TroposphereBlock(equation_type=t.equation_type, quality=t.quality,
                            hydrostatic=_value(43, t.hydrostatic, quantity=True),
                            coefficients=coefficients, residuals=residuals)


def _ionosphere(system: System, i: spartn.Ionosphere) -> IonosphereBlock:
    """[spec: SPARTN v2.0.3 Tables 6.19 to 6.22; Table 6.2 SF054 to SF067] an ionosphere block."""
    if i.equation_type not in _COEFFICIENT_COUNT:
        raise ValueError(f"ionosphere equation type {i.equation_type} is undefined (rule A8)")
    sats = satellites(system, i.mask_size, i.mask)
    if len(sats) != len(i.satellites):
        raise ValueError("ionosphere satellite blocks do not match the satellite mask")
    out = []
    for sat, s in zip(sats, i.satellites):
        block = s.small if s.coefficient_size == 0 else s.large
        residuals = None
        if s.residuals is not None:
            sf = _IONO_RESIDUAL_FIELD[s.residual_size]
            residuals = tuple(_residual(sf, r) for r in s.residuals)
        out.append(IonosphereSatellite(
            sat=sat, quality=s.quality,
            coefficients=_coefficients(i.equation_type, _IONO_FIELDS[s.coefficient_size],
                                       (block.c00, block.c01, block.c10, block.c11), False),
            residuals=residuals))
    return IonosphereBlock(equation_type=i.equation_type, satellites=tuple(out))


def _header(message: spartn.Message) -> tuple[StreamId, SpartnTag]:
    """[spec: SPARTN v2.0.3 Table 7.1 TF008 to TF011] the stream and the time tag as sent."""
    d = message.description
    return StreamId(d.solution_id, d.processor_id), SpartnTag(d.time_tag, bool(d.time_tag_type))


def from_hpac(message: spartn.Message) -> HpacInput:
    """An HPAC message, SM 1-0 to 1-4, in physical units.

    [spec: SPARTN v2.0.3 6.5.2, Tables 6.13 to 6.22; Table 6.1 for the
    subtype as the constellation] Raises ``TypeError`` for a message that is
    not HPAC and ``ValueError`` for an undefined value.
    """
    body = message.body
    if message.message_type != 1 or not isinstance(body, spartn.Hpac):
        raise TypeError("not a SPARTN HPAC message (SM 1)")
    system = System(message.description.subtype)
    stream, tag = _header(message)
    areas = []
    for a in body.areas:
        areas.append(HpacArea(
            area_id=a.area_id,
            grid_points=a.grid_points,
            troposphere=None if a.troposphere is None else _troposphere(a.troposphere),
            ionosphere=None if a.ionosphere is None else _ionosphere(system, a.ionosphere)))
    return HpacInput(system=system, stream=stream, tag=tag, siou=body.siou, aiou=body.aiou,
                     areas=tuple(areas))


def from_gad(message: spartn.Message) -> GadInput:
    """A GAD message, SM 2-0, in physical units; its time tag is in GPS time.

    [spec: SPARTN v2.0.3 6.5.3, Tables 6.23, 6.24; Table 6.2 SF031 to SF037,
    each scaled by 6.3: latitude 0.1° from −90°, longitude 0.1° from −180°,
    node counts from 1, spacings 0.1° from 0.1°] Raises ``TypeError`` for a
    message that is not GAD and ``ValueError`` for an undefined value or an
    area whose grid reaches beyond ±90° (``toolkit.area.extent``).
    """
    body = message.body
    if message.message_type != 2 or not isinstance(body, spartn.Gad):
        raise TypeError("not a SPARTN GAD message (SM 2)")
    stream, tag = _header(message)
    areas = []
    for a in body.areas:
        d = AreaDefinition(
            area_id=a.area_id,
            latitude=_value(32, a.latitude),
            longitude=_value(33, a.longitude),
            latitude_nodes=int(_value(34, a.latitude_nodes)),
            longitude_nodes=int(_value(35, a.longitude_nodes)),
            latitude_spacing=_value(36, a.latitude_spacing),
            longitude_spacing=_value(37, a.longitude_spacing))
        _area.extent(area_of(d))      # raises ValueError for a grid beyond the poles
        areas.append(d)
    return GadInput(stream=stream, tag=tag, siou=body.siou, aiou=body.aiou, areas=tuple(areas))
