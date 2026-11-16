"""Time scales: GPS, Galileo, BeiDou, GLONASS and UTC; weeks, days, SPARTN time tags.

DESIGN.md section 2.7, over the kernel ``kernels/src/timescale.c``.

A time of a scale is ``(week, sow)``: whole weeks and seconds since
1980-01-06 00:00:00 read on that scale's clock, 86 400 s to the day. For GPS
time these are the GPS week and time of week (IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4).
Galileo and BeiDou times in their own week numbering go through :func:`gst_to_gps` and
:func:`bdt_to_gps`. A calendar time (UTC, GLONASS time) is ``(days, sod)``:
days since 1980-01-06 and the second of that day, which reaches 86 400 only
inside a leap second; :func:`date_from_days` gives its civil date.

Arguments outside the range of their document raise
:class:`~ssrcore.toolkit.KernelError` with the kernel's status (DESIGN 2.7).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

from ssrcore._kernels import timescale as _k
from ssrcore.toolkit import KernelError

__all__ = [
    "System", "UtcParameters", "WEEK_SECONDS", "DAY_SECONDS",
    "normalize", "add", "diff", "resolve_week",
    "days_from_date", "date_from_days", "day_of_year", "glonass_day", "glonass_n4_nt",
    "gst_to_gps", "gps_to_gst", "ggto", "bdt_to_gps", "gps_to_bdt", "bdt_offset",
    "glonass_to_gps", "gps_to_glonass", "utc_to_gps", "gps_to_utc", "gst_to_utc", "bdt_to_utc",
    "spartn_tag_to_gps", "spartn_tag_from_gps", "spartn_resolve_half_day",
]


#: [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 3.3.4] one week, s; the kernel's constant.
WEEK_SECONDS: float = _k.WEEK_SECONDS
#: [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] one day, s; the kernel's constant.
DAY_SECONDS: float = _k.DAY_SECONDS


class System(IntEnum):
    """A constellation, numbered as the SPARTN message subtype that carries it.

    [spec: SPARTN Table 6.1, subtypes 0 to 4 of OCB and HPAC messages] The
    codes are the kernel's ``SSR_SYSTEM_*``. A SPARTN time tag is in the time
    scale of its message's system, GPS for a message of no or several
    constellations (Table 7.1 TF009). [design: one enumeration of systems for
    every layer; it lives at the lowest layer that needs it, and
    ``ssrcore.model.types`` re-exports it (DESIGN 3.1)]
    """

    GPS = _k.SYSTEM_GPS
    GLONASS = _k.SYSTEM_GLONASS
    GALILEO = _k.SYSTEM_GALILEO
    BEIDOU = _k.SYSTEM_BEIDOU
    QZSS = _k.SYSTEM_QZSS


def _checked(result: tuple, what: str):
    """[design: DESIGN 2.7, a kernel status other than success is a KernelError]"""
    status, value = result
    if status != _k.OK:
        raise KernelError("timescale", status, f"{what}: argument outside the range of its document")
    return value


@dataclass(frozen=True)
class UtcParameters:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-IX; Galileo OS SIS ICD 2.2 Table
    75; BDS ICD Table 5-18] Broadcast UTC parameters in SI units, the weeks as broadcast (8 bits).

    ``dt_ls``, ``dt_lsf`` leap seconds before and after the event (s);
    ``wn_lsf``, ``dn`` week and day of the event; ``a0`` (s), ``a1`` (s/s);
    ``wn_t``, ``t_ot`` reference week and time of ``a1`` (unused by BeiDou).
    """

    dt_ls: int
    dt_lsf: int
    wn_lsf: int
    dn: int
    a0: float = 0.0
    a1: float = 0.0
    wn_t: int = 0
    t_ot: float = 0.0


# ---- weeks and seconds of week ----------------------------------------------

def normalize(week: int, sow: float) -> tuple[int, float]:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] ``(week, sow)`` with sow in
    [0, 604800)."""
    return _k.normalize(week, sow)


def add(week: int, sow: float, dt: float) -> tuple[int, float]:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1] ``(week, sow) + dt`` seconds,
    normalized."""
    return _k.add(week, sow, dt)


def diff(week_a: int, sow_a: float, week_b: int, sow_b: float) -> float:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] ``a - b`` in seconds, weeks
    differenced first."""
    return _k.diff(week_a, sow_a, week_b, sow_b)


def resolve_week(truncated: int, bits: int, reference: int) -> int:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.1, 20.3.3.5.2.4; Galileo OS SIS
    ICD 2.2 5.1.2, 5.1.7, 5.1.8; BDS ICD 5.2.4.4, 5.2.4.18] The full week whose ``bits`` low bits
    are ``truncated``, within half the modulus of ``reference`` (same numbering).

    Broadcast widths: GPS WN 10, WNt and WN_LSF 8; Galileo WN 12, WN0t and
    WN_LSF 8, WN0G 6; BeiDou WN 13, WN_LSF 8.
    """
    return _checked(_k.week_resolve(truncated, bits, reference), "resolve_week")


# ---- calendar -----------------------------------------------------------------

def days_from_date(year: int, month: int, day: int) -> int:
    """[spec: TM-23 A.1.4 (A.4), (A.5)] Days from 1980-01-06 to the civil date,
    1900-03-01 to 2100-02-28."""
    return _checked(_k.days_from_date(year, month, day), "days_from_date")


def date_from_days(days: int) -> tuple[int, int, int]:
    """[spec: TM-23 A.1.4 (A.6), (A.7)] ``(year, month, day)`` of a day count."""
    return _checked(_k.date_from_days(days), "date_from_days")


def day_of_year(days: int) -> int:
    """[design: ordinal day, 1 January = 1] Day of year of a day count."""
    return _checked(_k.day_of_year(days), "day_of_year")


def glonass_day(n4: int, nt: int) -> int:
    """[spec: GLONASS ICD A.3.1.3] Day count of the GLONASS date (N4, NT)."""
    return _checked(_k.glonass_day(n4, nt), "glonass_day")


def glonass_n4_nt(days: int) -> tuple[int, int]:
    """[spec: GLONASS ICD A.3.1.3] ``(N4, NT)`` of a day count, from 1996."""
    return _checked(_k.glonass_n4_nt(days), "glonass_n4_nt")


# ---- scales --------------------------------------------------------------------

def gst_to_gps(wn: int, tow: float) -> tuple[int, float]:
    """[spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8] GPS ``(week, sow)`` of GST ``(WN, TOW)``,
    WN untruncated; the broadcast offset is :func:`ggto`."""
    return _k.gst_to_gps(wn, tow)


def gps_to_gst(week: int, sow: float) -> tuple[int, float]:
    """[spec: Galileo OS SIS ICD 2.2 5.1.2, 5.1.8] GST ``(WN, TOW)`` of GPS ``(week, sow)``."""
    return _k.gps_to_gst(week, sow)


def ggto(wn: int, tow: float, a0g: float, a1g: float, wn0g: int, t0g: float) -> float:
    """[spec: Galileo OS SIS ICD 2.2 5.1.8 Eq. 24] t_Galileo - t_GPS (s); WN0G as broadcast (6 bits)."""
    return _k.ggto(wn, tow, a0g, a1g, resolve_week(wn0g, 6, wn), t0g)


def bdt_to_gps(wn: int, sow_bdt: float) -> tuple[int, float]:
    """[spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] GPS ``(week, sow)`` of BDT
    ``(WN, SOW)``, WN untruncated; the broadcast offset is :func:`bdt_offset`."""
    return _k.bdt_to_gps(wn, sow_bdt)


def gps_to_bdt(week: int, sow: float) -> tuple[int, float]:
    """[spec: BDS ICD 3.3, 5.2.4.4; RTCM DF427] BDT ``(WN, SOW)`` of GPS ``(week, sow)``."""
    return _k.gps_to_bdt(week, sow)


def bdt_offset(a0: float, a1: float, sow_bdt: float) -> float:
    """[spec: BDS ICD 5.2.4.19 to 5.2.4.21] Delta t = A0 + A1 tE (s) of BDT against
    GPS, Galileo or GLONASS time: t_other = tE - Delta t."""
    return _k.bdt_offset(a0, a1, sow_bdt)


def glonass_to_gps(days: int, sod: float, leap_seconds: int) -> tuple[int, float]:
    """[spec: GLONASS ICD 3.3.3] GPS ``(week, sow)`` of a GLONASS ``(days, sod)``,
    given GPS - UTC in seconds."""
    return _k.glonass_to_gps(days, sod, leap_seconds)


def gps_to_glonass(week: int, sow: float, leap_seconds: int) -> tuple[int, float]:
    """[spec: GLONASS ICD 3.3.3] GLONASS ``(days, sod)`` of GPS ``(week, sow)``."""
    return _k.gps_to_glonass(week, sow, leap_seconds)


def utc_to_gps(days: int, sod: float, leap_seconds: int) -> tuple[int, float]:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a, A0 = A1 = 0] GPS ``(week,
    sow)`` of a UTC ``(days, sod)``; 23:59:60 is ``sod`` 86400 with the count before the leap."""
    return _k.utc_to_gps(days, sod, leap_seconds)


def gps_to_utc(week: int, sow: float, p: UtcParameters) -> tuple[int, float, int]:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] UTC ``(days, sod, case)`` of
    GPS ``(week, sow)``; case 1, 2, 3 is a, b (leap second window), c."""
    return _checked(_k.gps_utc(week, sow, p.dt_ls, p.dt_lsf, resolve_week(p.wn_lsf, 8, week),
                               p.dn, p.a0, p.a1, resolve_week(p.wn_t, 8, week), p.t_ot), "gps_to_utc")


def gst_to_utc(wn: int, tow: float, p: UtcParameters) -> tuple[int, float, int]:
    """[spec: Galileo OS SIS ICD 2.2 5.1.7 Eq. 21 to 23] UTC ``(days, sod, case)``
    of GST ``(WN, TOW)``; case a with delta t_LS (Eq. 21), as in IS-GPS-200N with IRN-IS-200N-003
    and -004 20.3.3.5.2.4 a. Issue 2.1 printed delta t_LSF there (its Eq. 20) and 2.2
    corrects it (``ssr_timescale_gst_utc``)."""
    return _checked(_k.gst_utc(wn, tow, p.dt_ls, p.dt_lsf, resolve_week(p.wn_lsf, 8, wn),
                               p.dn, p.a0, p.a1, resolve_week(p.wn_t, 8, wn), p.t_ot), "gst_to_utc")


def bdt_to_utc(wn: int, sow_bdt: float, p: UtcParameters) -> tuple[int, float, int]:
    """[spec: BDS ICD 5.2.4.18] UTC ``(days, sod, case)`` of BDT ``(WN, SOW)``; DN 0 to 6."""
    return _checked(_k.bdt_utc(wn, sow_bdt, p.dt_ls, p.dt_lsf, resolve_week(p.wn_lsf, 8, wn),
                               p.dn, p.a0, p.a1), "bdt_to_utc")


# ---- SPARTN time tags -----------------------------------------------------------

def _leap(scale: System, leap_seconds: int | None) -> int:
    """[design: DESIGN T4, a GLONASS time is not converted without the leap seconds]"""
    if leap_seconds is None:
        if scale == System.GLONASS:
            raise ValueError("a GLONASS time tag needs GPS - UTC leap seconds")
        return 0
    return leap_seconds


def spartn_tag_to_gps(tag: int, scale: System, leap_seconds: int | None = None) -> tuple[int, float]:
    """[spec: SPARTN Table 7.1 TF009] GPS ``(week, sow)`` of a 32-bit time tag,
    seconds since 2010-01-01 00:00:00 of the message's scale. GLONASS needs
    ``leap_seconds`` (GPS - UTC); its epoch is read on the GLONASS clock
    (DESIGN Q2, proposed)."""
    return _checked(_k.spartn_tag32_to_gps(tag, scale, _leap(scale, leap_seconds)), "spartn_tag_to_gps")


def spartn_tag_from_gps(week: int, sow: float, scale: System,
                        leap_seconds: int | None = None) -> int:
    """[spec: SPARTN Table 7.1 TF009] The 32-bit time tag, in the message's scale,
    of a GPS time, whole seconds rounded down."""
    return _checked(_k.spartn_tag32_from_gps(week, sow, scale, _leap(scale, leap_seconds)),
                    "spartn_tag_from_gps")


def spartn_resolve_half_day(tag16: int, reference: int) -> int:
    """[spec: SPARTN Table 7.1 TF009; design: nearest the reference, DESIGN T2]
    The 32-bit tag of a 16-bit half-day tag, ``reference`` a 32-bit tag of the
    same scale (for example of the receipt time, :func:`spartn_tag_from_gps`)."""
    return _checked(_k.spartn_resolve_tag16(tag16, reference), "spartn_resolve_half_day")
