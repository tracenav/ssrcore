"""gnsswire.rtcm3 ephemeris messages (1019, 1020, 1042, 1044, 1045, 1046) to ephemeris inputs.

DESIGN.md section 2.2. Each message's fields are the integers the satellite
broadcast, with the widths and scale factors of the system's interface
document (RTCM 3.5.8, 3.5.9, 3.5.18, 3.5.20); the adapter multiplies each by
its resolution, keeps angles in semicircles (``KeplerEphemeris``), and places
the reference times on the GPS time scale.

Two things a message alone does not fix are passed in: ``near``, a GPS time
close to the message's receipt (the receipt time, or the time base's current
time, rule T5), which resolves the week-number rollover and the week or day
of a reference time; and for GLONASS the leap seconds GPS − UTC (rule T4).
The offsets between the time scales, the week rollover and the GLONASS
calendar are the toolkit's (``toolkit.timescale``).
"""

from __future__ import annotations

from gnsswire import rtcm3

from ssrcore.model.ephemeris import beidou_iodc
from ssrcore.model.types import (
    WEEK_SECONDS,
    GlonassEphemeris,
    GnssTime,
    KeplerEphemeris,
    Sat,
    System,
)
from ssrcore.toolkit import timescale

__all__ = [
    "EphemerisError",
    "LeapSecondsNeeded",
    "ephemeris",
    "from_1019",
    "from_1020",
    "from_1042",
    "from_1044",
    "from_1045",
    "from_1046",
]

_HALF_WEEK = WEEK_SECONDS / 2.0
_DAY = timescale.DAY_SECONDS

# The offsets between GPS time and GST, BDT and GLONASS time are the
# toolkit's (``toolkit.timescale``, DESIGN 2.7), not repeated here.
# [spec: RTCM DF076, DF289, DF489] broadcast week-number widths, bits.
_GPS_WEEK_BITS, _GALILEO_WEEK_BITS, _BEIDOU_WEEK_BITS = 10, 12, 13
# [spec: GLONASS ICD 4.4 Table 4.3 (word P1); RTCM DF106] minutes between
# adjacent values of tb for each value of P1, as seconds.
_GLONASS_P1_SECONDS = (0.0, 30 * 60.0, 45 * 60.0, 60 * 60.0)
# [spec: RTCM DF110, tb from 15 to 1425 min in steps of 15 min; GLONASS ICD 4.4]
_TB_STEP = 15 * 60.0
_KM = 1000.0


class EphemerisError(ValueError):
    """A message whose fields cannot form an ephemeris; the text says which rule.

    [design: DESIGN 2.1, the engine turns it into a rejection with its reason]
    """


class LeapSecondsNeeded(LookupError):
    """A GLONASS ephemeris cannot be placed on the GPS time scale without the leap seconds.

    [spec: DESIGN 3.3 rule T4, GLONASS inputs are held until they are known]
    """


def _scaled(value: int, exponent: int) -> float:
    """value × 2^exponent, exact in binary64. [spec: RTCM 3.5, field resolutions]"""
    return value * 2.0 ** exponent


def _full_week(broadcast: int, bits: int, near_week: int) -> int:
    """The full week of a broadcast week number, its rollover resolved nearest
    ``near_week`` of the same numbering (GPS, GST or BDT weeks).

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1, the transmitted week number is
    the full week modulo 1024 (and RTCM DF289, DF489 modulo 4096, 8192);
    ``toolkit.timescale.resolve_week``]
    """
    return timescale.resolve_week(broadcast, bits, near_week)


def _place(week: int, seconds: float, near: GnssTime) -> GnssTime:
    """A seconds-of-week reference time in the week that puts it within half a week of ``near``.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 Table 20-IV, tk accounts for the
    beginning or end of week crossover: the reference time and the time of use are less
    than half a week apart, so a toe or toc broadcast near the end of a week
    may belong to the next week]
    """
    t = GnssTime.normalised(week, seconds)
    dt = t - near
    if dt > _HALF_WEEK:
        t = t - WEEK_SECONDS
    elif dt < -_HALF_WEEK:
        t = t + WEEK_SECONDS
    return t


def _sat(system: System, number: int) -> Sat:
    """[spec: RTCM DF009, DF038, DF252, DF488; SPARTN satellite numbers (types.Sat)]"""
    try:
        return Sat(system, number)
    except ValueError as e:  # 0 is "unknown" (DF038) and has no satellite
        raise EphemerisError(f"satellite ID {number}: {e}") from None


def from_1019(m: rtcm3.Msg1019, near: GnssTime) -> KeplerEphemeris:
    """A GPS LNAV ephemeris from RTCM 1019.

    [spec: RTCM 3.5.8 Table 3.5-21, DF076 to DF103, DF137 (the resolutions of
    IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-III); IS-GPS-200N with IRN-IS-200N-003 and
    -004 20.3.3.4.1, an IODE that differs from the 8 LSBs of the IODC marks a cutover, and the data
    are not used]
    """
    if m.iode != m.iodc & 0xFF:
        raise EphemerisError(f"IODE {m.iode} differs from the 8 LSBs of IODC {m.iodc} "
                             "(IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.1)")
    week = _full_week(m.week, _GPS_WEEK_BITS, near.week)
    return KeplerEphemeris(
        sat=_sat(System.GPS, m.satellite), message=1019, iod=m.iode,
        toe=_place(week, m.toe * 16.0, near), toc=_place(week, m.toc * 16.0, near),
        sqrt_a=_scaled(m.sqrt_a, -19), e=_scaled(m.e, -33), m0=_scaled(m.m0, -31),
        delta_n=_scaled(m.delta_n, -43), omega=_scaled(m.omega, -31),
        omega0=_scaled(m.omega0, -31), omega_dot=_scaled(m.omega_dot, -43),
        i0=_scaled(m.i0, -31), idot=_scaled(m.idot, -43),
        cuc=_scaled(m.cuc, -29), cus=_scaled(m.cus, -29), crc=_scaled(m.crc, -5),
        crs=_scaled(m.crs, -5), cic=_scaled(m.cic, -29), cis=_scaled(m.cis, -29),
        af0=_scaled(m.af0, -31), af1=_scaled(m.af1, -43), af2=_scaled(m.af2, -55),
        group_delays=(_scaled(m.tgd, -31),), health=m.health, accuracy=m.accuracy,
        fit_interval=m.fit_interval, iodc=m.iodc,
    )


# [spec: RTCM DF454, QZSS SV health: bit 0 (MSb) NAV data, bits 1 to 4 the
# L1 C/A, L2C, L5 and L1C signals, bit 5 (LSb) "reserved for LEX"]
_QZSS_HEALTH_DEFINED = 0b111110


def from_1044(m: rtcm3.Msg1044, near: GnssTime) -> KeplerEphemeris:
    """A QZSS ephemeris from RTCM 1044, on the GPS time scale.

    [spec: RTCM 3.5.19 Table 3.5-112, DF429 to DF457, with the resolutions of
    their rows (the same as GPS's of Table 3.5-21); DF429, satellite IDs 1 to
    10 are PRN 193 to 202, the numbers SPARTN SF095 gives QZSS satellites;
    DF452, the week is the GPS week modulo 1024; DF428 and DESIGN T1, QZSS
    time is GPS time; DF456, the IODE and the 8 LSBs of the IODC of one data
    set are equal, so a message where they differ is refused as for GPS;
    DF454, the reserved LSb of the health is not a health flag]
    """
    if m.iode != m.iodc & 0xFF:
        raise EphemerisError(f"IODE {m.iode} differs from the 8 LSBs of IODC {m.iodc} (RTCM DF456)")
    week = _full_week(m.week, _GPS_WEEK_BITS, near.week)
    return KeplerEphemeris(
        sat=_sat(System.QZSS, m.satellite), message=1044, iod=m.iode,
        toe=_place(week, m.toe * 16.0, near), toc=_place(week, m.toc * 16.0, near),
        sqrt_a=_scaled(m.sqrt_a, -19), e=_scaled(m.e, -33), m0=_scaled(m.m0, -31),
        delta_n=_scaled(m.delta_n, -43), omega=_scaled(m.omega, -31),
        omega0=_scaled(m.omega0, -31), omega_dot=_scaled(m.omega_dot, -43),
        i0=_scaled(m.i0, -31), idot=_scaled(m.idot, -43),
        cuc=_scaled(m.cuc, -29), cus=_scaled(m.cus, -29), crc=_scaled(m.crc, -5),
        crs=_scaled(m.crs, -5), cic=_scaled(m.cic, -29), cis=_scaled(m.cis, -29),
        af0=_scaled(m.af0, -31), af1=_scaled(m.af1, -43), af2=_scaled(m.af2, -55),
        group_delays=(_scaled(m.tgd, -31),), health=m.health & _QZSS_HEALTH_DEFINED,
        accuracy=m.accuracy, fit_interval=m.fit_interval, iodc=m.iodc,
    )


def _galileo_health(*signals: tuple[int, int]) -> int:
    """Signal health status and data validity status of each signal, packed so
    that 0 is nominal: per signal (SHS << 1) | DVS, three bits, first signal highest.

    [spec: Galileo OS SIS ICD 2.2 5.1.9.3, Tables 81 and 84 (0 is "Navigation
    data valid" and "Signal OK"); design: the packing keeps every flag for the
    record and makes E5 one test]
    """
    health = 0
    for shs, dvs in signals:
        health = (health << 3) | (shs << 1) | dvs
    return health


def _galileo(m, message: int, health: int, group_delays: tuple[float, ...],
             near: GnssTime) -> KeplerEphemeris:
    """[spec: RTCM 3.5.18 Tables 3.5-110, 3.5-111, DF286 to DF317; Galileo OS SIS ICD 2.2 5.1.1, 5.1.3;
    5.1.2 for GST, through ``toolkit.timescale.gst_to_gps``]"""
    wn = _full_week(m.week, _GALILEO_WEEK_BITS, timescale.gps_to_gst(near.week, near.seconds)[0])

    def gps(tow: float) -> GnssTime:
        return _place(*timescale.gst_to_gps(wn, tow), near)

    return KeplerEphemeris(
        sat=_sat(System.GALILEO, m.satellite), message=message, iod=m.iodnav,
        toe=gps(m.toe * 60.0), toc=gps(m.toc * 60.0),
        sqrt_a=_scaled(m.sqrt_a, -19), e=_scaled(m.e, -33), m0=_scaled(m.m0, -31),
        delta_n=_scaled(m.delta_n, -43), omega=_scaled(m.omega, -31),
        omega0=_scaled(m.omega0, -31), omega_dot=_scaled(m.omega_dot, -43),
        i0=_scaled(m.i0, -31), idot=_scaled(m.idot, -43),
        cuc=_scaled(m.cuc, -29), cus=_scaled(m.cus, -29), crc=_scaled(m.crc, -5),
        crs=_scaled(m.crs, -5), cic=_scaled(m.cic, -29), cis=_scaled(m.cis, -29),
        af0=_scaled(m.af0, -34), af1=_scaled(m.af1, -46), af2=_scaled(m.af2, -59),
        group_delays=group_delays, health=health, accuracy=m.sisa,
    )


def from_1045(m: rtcm3.Msg1045, near: GnssTime) -> KeplerEphemeris:
    """A Galileo F/NAV ephemeris from RTCM 1045; health from E5a (DF314, DF315).

    [spec: RTCM 3.5.18 Table 3.5-110; Galileo OS SIS ICD 2.2 5.1.9.3]
    """
    return _galileo(m, 1045, _galileo_health((m.e5a_health, m.e5a_validity)),
                    (_scaled(m.bgd_e5a_e1, -32),), near)


def from_1046(m: rtcm3.Msg1046, near: GnssTime) -> KeplerEphemeris:
    """A Galileo I/NAV ephemeris from RTCM 1046; health from E1-B (DF287, DF288)
    and E5b (DF316, DF317), packed E1-B first.

    [spec: RTCM 3.5.18 Table 3.5-111; Galileo OS SIS ICD 2.2 5.1.9.3]
    """
    health = _galileo_health((m.e1b_health, m.e1b_validity), (m.e5b_health, m.e5b_validity))
    return _galileo(m, 1046, health, (_scaled(m.bgd_e5a_e1, -32), _scaled(m.bgd_e5b_e1, -32)), near)


def from_1042(m: rtcm3.Msg1042, near: GnssTime) -> KeplerEphemeris:
    """A BeiDou D1/D2 ephemeris from RTCM 1042, toe and toc moved from BDT to GPS time,
    with the issue of data SPARTN names (IODC from toc, rule E1).

    [spec: RTCM 3.5.20 Table 3.5-113, DF488 to DF515; BDS ICD 5.2.4; SPARTN
    Table 6.2 SF100]
    """
    wn = _full_week(m.week, _BEIDOU_WEEK_BITS, timescale.gps_to_bdt(near.week, near.seconds)[0])
    toc_bdt = m.toc * 8

    def gps(sow_bdt: float) -> GnssTime:
        return _place(*timescale.bdt_to_gps(wn, sow_bdt), near)

    return KeplerEphemeris(
        sat=_sat(System.BEIDOU, m.satellite), message=1042, iod=beidou_iodc(toc_bdt),
        toe=gps(m.toe * 8.0), toc=gps(float(toc_bdt)),
        sqrt_a=_scaled(m.sqrt_a, -19), e=_scaled(m.e, -33), m0=_scaled(m.m0, -31),
        delta_n=_scaled(m.delta_n, -43), omega=_scaled(m.omega, -31),
        omega0=_scaled(m.omega0, -31), omega_dot=_scaled(m.omega_dot, -43),
        i0=_scaled(m.i0, -31), idot=_scaled(m.idot, -43),
        cuc=_scaled(m.cuc, -31), cus=_scaled(m.cus, -31), crc=_scaled(m.crc, -6),
        crs=_scaled(m.crs, -6), cic=_scaled(m.cic, -31), cis=_scaled(m.cis, -31),
        af0=_scaled(m.a0, -33), af1=_scaled(m.a1, -50), af2=_scaled(m.a2, -66),
        group_delays=(m.tgd1 * 1e-10, m.tgd2 * 1e-10), health=m.health, accuracy=m.urai,
    )


def from_1020(m: rtcm3.Msg1020, near: GnssTime, leap_seconds: int | None) -> GlonassEphemeris:
    """A GLONASS ephemeris from RTCM 1020, tb placed on the GPS time scale.

    tb (DF110) counts 15-minute intervals of the GLONASS day (UTC(SU) + 3 h);
    the day is the one that puts tb within half a day of ``near``, and the
    calendar placement and the conversion to GPS time are the toolkit's
    (``timescale.gps_to_glonass``, ``glonass_to_gps``, with GPS − UTC
    ``leap_seconds``). The issue of data is tb itself, the 7
    bits SPARTN SF019 sends. Health packs the Bn MSB (bit 0) and, for a
    GLONASS-M satellite (DF130 = 1), ln of the third string (bit 1).

    [spec: RTCM 3.5.9 Table 3.5-22, DF038 to DF136 (2^-11 km, 2^-20 km/s,
    2^-30 km/s², γn 2^-40, τn and Δτn 2^-30 s); GLONASS ICD 3.3.3, 4.4 and
    Table 4.3; RTCM DF392, the GLONASS IOD is the 7-bit tb (DESIGN Q3);
    design: the day nearest ``near``, since an ephemeris is broadcast within
    an hour of its tb; DF129 NT and DF134 N4 are not needed for that, and
    are not sent by every station (DF134 is 0 when unknown)]
    """
    if leap_seconds is None:
        raise LeapSecondsNeeded("a GLONASS ephemeris needs GPS − UTC (rule T4)")
    sat = _sat(System.GLONASS, m.satellite)
    channel = m.channel - 7
    if not -7 <= channel <= 6:
        raise EphemerisError(f"frequency channel {channel} outside -7 to 6 (GLONASS ICD 3.3.1.1)")
    if not 1 <= m.tb <= 95:
        raise EphemerisError(f"tb {m.tb} outside 1 to 95 (RTCM DF110)")
    days, sod = timescale.gps_to_glonass(near.week, near.seconds, leap_seconds)
    tb_sod = m.tb * _TB_STEP
    day = min((days - 1, days, days + 1), key=lambda d: abs((d - days) * _DAY + tb_sod - sod))
    tb = GnssTime(*timescale.glonass_to_gps(day, tb_sod, leap_seconds))
    health = m.bn_msb | ((m.ln_third << 1) if m.m == 1 else 0)
    return GlonassEphemeris(
        sat=sat, frequency_channel=channel, iod=m.tb, tb=tb,
        position=tuple(_scaled(v, -11) * _KM for v in (m.x, m.y, m.z)),
        velocity=tuple(_scaled(v, -20) * _KM for v in (m.x_dot, m.y_dot, m.z_dot)),
        acceleration=tuple(_scaled(v, -30) * _KM for v in (m.x_ddot, m.y_ddot, m.z_ddot)),
        tau=_scaled(m.tau, -30), gamma=_scaled(m.gamma, -40), delta_tau=_scaled(m.delta_tau, -30),
        health=health, p1=_GLONASS_P1_SECONDS[m.p1],
    )


def ephemeris(m, near: GnssTime, leap_seconds: int | None = None) -> KeplerEphemeris | GlonassEphemeris:
    """The ephemeris of a decoded RTCM ephemeris message (DESIGN 2.2).

    Raises ``EphemerisError`` for a message that cannot form one,
    ``LeapSecondsNeeded`` for a GLONASS message without leap seconds, and
    ``TypeError`` for a message that is not one of the six.
    [spec: RTCM 3.5.8, 3.5.9, 3.5.18, 3.5.19, 3.5.20; design: DESIGN 2.1, a type the
    engine does not know is a programming error]
    """
    if isinstance(m, rtcm3.Msg1019):
        return from_1019(m, near)
    if isinstance(m, rtcm3.Msg1020):
        return from_1020(m, near, leap_seconds)
    if isinstance(m, rtcm3.Msg1042):
        return from_1042(m, near)
    if isinstance(m, rtcm3.Msg1044):
        return from_1044(m, near)
    if isinstance(m, rtcm3.Msg1045):
        return from_1045(m, near)
    if isinstance(m, rtcm3.Msg1046):
        return from_1046(m, near)
    raise TypeError(f"not an RTCM ephemeris message: {type(m).__name__}")
