"""gnsswire.spartn OCB messages (SM 0-0 to 0-4) to correction inputs.

DESIGN.md section 2.2. A decoded ``gnsswire.spartn.Message`` whose body is an
``Ocb`` holds the raw integers of SPARTN 6.5.1 (Tables 6.3 to 6.12); this
adapter gives the ``OcbInput`` of the same message in physical units, with
each satellite named by its mask position (Table 6.2 SF011, SF012, SF093 to
SF095) and each bias by its signal (SF025 to SF028, SF102 to SF107). The
scaling is the toolkit's (``toolkit.ssr.field``, SPARTN 6.3 and Table 6.2);
no bits are parsed here.

The time tag is kept as sent (``SpartnTag``): resolving and converting it is
the time base's (rules T1 to T3).
"""

from __future__ import annotations

from gnsswire import spartn

from ssrcore.model.types import (
    ClockCorrection,
    CodeBias,
    OcbHeaderInput,
    OcbInput,
    OcbSatellite,
    OrbitCorrection,
    PhaseBias,
    Sat,
    Signal,
    SpartnTag,
    StreamId,
    System,
)
from ssrcore.toolkit import ssr

# [spec: SPARTN Table 6.1] message type of OCB.
OCB_MESSAGE_TYPE = 0

# [spec: SPARTN Table 6.2 SF011, SF012, SF093, SF094, SF095] satellite mask
# widths by system and the 2-bit size indicator.
SATELLITE_MASK_BITS = {
    System.GPS: (32, 44, 56, 64), System.GLONASS: (24, 36, 48, 63), System.GALILEO: (36, 45, 54, 64),
    System.BEIDOU: (37, 46, 55, 64), System.QZSS: (10, 40, 48, 64),
}

# [spec: SPARTN Table 6.2 SF025 to SF028, SF102 to SF107] phase and code bias
# mask widths by system and the 1-bit size indicator (the same for both masks).
BIAS_MASK_BITS = {
    System.GPS: (6, 11), System.GLONASS: (5, 9), System.GALILEO: (8, 15),
    System.BEIDOU: (8, 15), System.QZSS: (6, 11),
}

# [spec: SPARTN Table 6.2, the SF numbers of the fields scaled here]
_SF_CORRECTION, _SF_YAW, _SF_IODE_CONTINUITY, _SF_CONTINUITY, _SF_CODE_BIAS = 20, 21, 22, 15, 29


def mask_positions(mask: int, bits: int) -> tuple[int, ...]:
    """The positions set in a mask of ``bits`` bits, position 0 the leftmost (first sent) bit.

    [spec: SPARTN Table 6.2 SF011, "Bitmask position 0 is the leftmost bit of
    the satellite mask"; SF025 to SF028, "Bit #'s left to right"]
    """
    return tuple(k for k in range(bits) if (mask >> (bits - 1 - k)) & 1)


def satellites(system: System, mask_size: int, mask: int) -> tuple[Sat, ...]:
    """The satellites a satellite mask names, in mask order; one function for the
    OCB and HPAC adapters.

    [spec: SPARTN Table 6.2 SF011, SF012, SF093 to SF095: bitmask position 0 is
    the leftmost bit; position + 1 is the PRN, the GLONASS orbit slot, or the
    SPARTN QZSS number]
    """
    width = SATELLITE_MASK_BITS[System(system)][mask_size]
    if mask >> width:
        raise ValueError(f"satellite mask wider than its {width} bits")
    return tuple(Sat(System(system), k + 1) for k in mask_positions(mask, width))


def _value(sf: int, encoded: int) -> float | None:
    # [spec: SPARTN 6.3, special values are checked before the conversion; design:
    # a value that is special or outside its range carries no number, and the
    # block holding it is treated as not sent, so that the other blocks stay usable]
    f = ssr.field(sf, encoded)
    return f.value if f.status == ssr.Status.OK else None


def orbit_correction(orbit: spartn.Orbit, yaw_present: bool, ephemeris_type: int = 0) -> OrbitCorrection | None:
    """An orbit block (SPARTN Table 6.5) in metres and degrees; None if a component is unusable.

    Yaw (SF021) is kept only when the header's SF008 says it is present and
    the value is a valid angle; otherwise it is None, and rule C8 uses zero.
    The header's ephemeris type is kept with the block (rule E2).
    [spec: SPARTN Table 6.2 SF008, SF016 to SF021, SF096 to SF101; 8.14]
    """
    components = [_value(_SF_CORRECTION, x) for x in (orbit.radial, orbit.along, orbit.cross)]
    if any(c is None for c in components):
        return None
    yaw = None
    if yaw_present and orbit.yaw is not None:
        yaw = _value(_SF_YAW, orbit.yaw)
    return OrbitCorrection(orbit.iode, *components, yaw=yaw, ephemeris_type=ephemeris_type)


def clock_correction(clock: spartn.Clock) -> ClockCorrection | None:
    """A clock block (SPARTN Table 6.6): IODE continuity in seconds, δC in metres, URE index.

    None if a value is unusable. [spec: SPARTN Table 6.2 SF020, SF022, SF024]
    """
    continuity = _value(_SF_IODE_CONTINUITY, clock.iode_continuity)
    correction = _value(_SF_CORRECTION, clock.correction)
    if continuity is None or correction is None:
        return None
    return ClockCorrection(continuity, correction, clock.ure)


def biases(system: System, bias: spartn.Bias) -> tuple[tuple[PhaseBias, ...], tuple[CodeBias, ...]]:
    """A bias block (SPARTN Tables 6.7 to 6.12): phase and code biases by signal, in mask order.

    A bias at a spare mask position names no known signal and is left out;
    so is one whose value is unusable. [spec: SPARTN Table 6.2 SF015, SF020,
    SF023, SF025 to SF029, SF102 to SF107; design: DESIGN 3.3 C9, a signal is
    only ever its own RINEX code]
    """
    phase: list[PhaseBias] = []
    width = BIAS_MASK_BITS[system][bias.phase_mask_size]
    for position, b in zip(mask_positions(bias.phase_mask, width), bias.phase):
        signal = Signal.from_spartn(system, position)
        value = _value(_SF_CORRECTION, b.correction)
        continuity = _value(_SF_CONTINUITY, b.continuity)
        if signal is not None and value is not None and continuity is not None:
            phase.append(PhaseBias(signal, value, bool(b.fix), continuity))
    code: list[CodeBias] = []
    width = BIAS_MASK_BITS[system][bias.code_mask_size]
    for position, raw in zip(mask_positions(bias.code_mask, width), bias.code):
        signal = Signal.from_spartn(system, position)
        value = _value(_SF_CODE_BIAS, raw)
        if signal is not None and value is not None:
            code.append(CodeBias(signal, value))
    return tuple(phase), tuple(code)


def ocb_input(message: spartn.Message) -> OcbInput:
    """The ``OcbInput`` of a decoded OCB message, SM 0-0 to 0-4.

    Raises ``TypeError`` for a message that is not OCB (a programming error,
    DESIGN 2.1). [spec: SPARTN 6.5.1 Tables 6.3 to 6.12; Table 7.1 TF007 to
    TF011]
    """
    body = message.body
    if message.message_type != OCB_MESSAGE_TYPE or not isinstance(body, spartn.Ocb):
        raise TypeError(f"not an OCB message: SM {message.message_type}-{message.description.subtype}")
    d = message.description
    system = System(d.subtype)
    blocks = []
    for sat, s in zip(satellites(system, body.mask_size, body.mask), body.satellites):
        if s.dnu:
            # [spec: SPARTN Table 6.2 SF013, no other field of the block is present]
            blocks.append(OcbSatellite(sat, True))
            continue
        phase = code = None
        if s.bias is not None:
            phase, code = biases(system, s.bias)
        blocks.append(OcbSatellite(
            sat, False,
            continuity=_value(_SF_CONTINUITY, s.continuity),
            orbit=None if s.orbit is None else orbit_correction(s.orbit, bool(body.yaw_present),
                                                                 body.ephemeris_type),
            clock=None if s.clock is None else clock_correction(s.clock),
            phase_biases=phase, code_biases=code))
    return OcbInput(system, StreamId(d.solution_id, d.processor_id), SpartnTag(d.time_tag, bool(d.time_tag_type)),
                    body.siou, bool(body.end_of_set), bool(body.yaw_present), body.datum, body.ephemeris_type,
                    tuple(blocks))


def header_input(frame: spartn.Frame, header: spartn.OcbHeader) -> OcbHeaderInput:
    """The ``OcbHeaderInput`` of an OCB frame whose subtype has no table, from its
    header as ``gnsswire.spartn.decode_ocb_header`` reads it.

    [spec: SPARTN 8.4 item 5; Table 6.3 SF005, SF009, SF010; Table 7.1 TF007
    to TF011]
    """
    d = frame.description
    return OcbHeaderInput(d.subtype, StreamId(d.solution_id, d.processor_id),
                          SpartnTag(d.time_tag, bool(d.time_tag_type)), header.siou,
                          bool(header.end_of_set), header.datum)


def dnu_table(message: spartn.Message) -> str | None:
    """Which layout of a do-not-use satellite block an OCB message was read with:
    ``"6.2"`` (SF013 alone), ``"8.2"`` (SF013 with SF014 and SF015), or None when
    no satellite of the message is flagged.

    [spec: SPARTN Table 6.2 SF013 ("no fields in the associated satellite
    block will be present") and Table 8.2 (nominal sizes, 7 bits for a DNU
    satellite); gnsswire ``dnu_reading`` says which reading fitted the payload]
    """
    body = message.body
    if not isinstance(body, spartn.Ocb) or not any(s.dnu for s in body.satellites):
        return None
    return "8.2" if spartn.dnu_reading(body) == spartn.DNU_TABLE_8_2 else "6.2"
