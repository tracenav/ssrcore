"""Broadcast ephemeris store: issues per satellite, validity, lookup by issue of data.

DESIGN.md sections 2.5 (validity of a broadcast ephemeris), 3.2 (the store)
and 3.3, rules E1 to E5. The store keeps, per satellite and ephemeris
message, the three issues with the latest reference times; a SPARTN orbit
block finds its ephemeris by the issue of data it names (E1) and the
ephemeris type its OCB header names (E2), and the store answers either with
the ephemeris or with the reason there is none.

The store holds no clock: every question about validity takes the time it
is asked for (DESIGN 2.1, I/O at the edge).
"""

from __future__ import annotations

from dataclasses import dataclass

from ssrcore.model.types import (
    FeedOutcome,
    FeedReason,
    FeedResult,
    GlonassEphemeris,
    GnssTime,
    KeplerEphemeris,
    Missing,
    Sat,
    System,
)
from ssrcore.toolkit import orbits, timescale

Ephemeris = KeplerEphemeris | GlonassEphemeris

# [spec: RTCM Table 3.5-22] the message a GlonassEphemeris comes from; it has
# no ``message`` field because 1020 is the only GLONASS ephemeris message.
GLONASS_MESSAGE = 1020

# [design: DESIGN 3.2, three issues cover the provider switching IODE up to
# 60 s after the receiver has a new issue (RTCM 3.5.13) and a receiver that
# has a new issue before the provider switches]
DEPTH = 3

# Half-widths of the validity intervals of DESIGN 2.5, seconds.
# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1 (fit interval flag 0 = 4 hours) and
# Table 20-XIII (toe is 2 h after the start of a 4 h fit, its midpoint)]
_GPS_FIT_FLAG_0 = 4 * 3600.0
# [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-XII, curve fit intervals for a fit
# interval flag of 1 by IODC range: 8 h for 240-247, 14 h for 248-255 and 496, 26 h for
# 497-503 and 1021-1023; the IODC values of its Note 2 belong to the 6 h fit]
_GPS_FIT_BY_IODC = (((240, 247), 8 * 3600.0), ((248, 255), 14 * 3600.0), ((496, 496), 14 * 3600.0),
                    ((497, 503), 26 * 3600.0), ((1021, 1023), 26 * 3600.0))
_GPS_FIT_FLAG_1 = 6 * 3600.0


def gps_fit_interval(fit_flag: int, iodc: int | None) -> float:
    """The curve fit interval of a GPS ephemeris, seconds.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1 (flag 0: 4 hours) and Table 20-XII
    (flag 1: 6, 8, 14 or 26 hours by IODC); design: without an IODC, or with one
    the table gives no longer fit for (Note 2, or reserved by Note 6), the
    shortest fit greater than 4 hours (6 h), so that no ephemeris is used
    beyond its fit]
    """
    if not fit_flag:
        return _GPS_FIT_FLAG_0
    if iodc is not None:
        for (low, high), length in _GPS_FIT_BY_IODC:
            if low <= iodc <= high:
                return length
    return _GPS_FIT_FLAG_1

# [spec: RTCM DF457, QZSS fit interval 0: "curve-fit interval is 2 hours", 1:
# "greater than 2 hours"; design: IS-QZSS, which would give the longer fits,
# is not staged (DESIGN Q14), so both flags take the 2 hours the field
# guarantees, and no ephemeris is used beyond its fit]
_QZSS_FIT = 2 * 3600.0
# [design: DESIGN 2.5, the Galileo OS SIS ICD 2.2 gives no fit interval (5.1.1,
# 5.1.9.2: a "validity duration" without a value); 4 h covers the nominal 3 h
# issue with margin]
_GALILEO_HALF = 4 * 3600.0
# [design: DESIGN 2.5, the B1I ICD gives no fit interval; ephemerides are
# issued hourly]
_BEIDOU_HALF = 2 * 3600.0
# [design: DESIGN 2.5, GLONASS ICD 4.4 gives tb and its interval P1, not a
# validity; 15 min plus half of P1 around tb]
_GLONASS_BASE = 15 * 60.0


@dataclass(frozen=True, slots=True)
class StoredEphemeris:
    """One stored issue and the time it was first received (None when the caller
    gave none). [design: DESIGN 3.2, the store keeps the receipt time with each issue]"""

    ephemeris: Ephemeris
    received: GnssTime | None

    @property
    def message(self) -> int:
        """[spec: DESIGN 2.2 table; RTCM Tables 3.5-21, 3.5-22, 3.5-110, 3.5-111, 3.5-113]
        the RTCM message number the issue came from."""
        return message_of(self.ephemeris)


@dataclass(frozen=True, slots=True)
class EphemerisLookup:
    """The answer to a lookup: the issue, or the reason there is none.

    [design: DESIGN 2.3, a term that is not valid is absent with a reason]
    """

    stored: StoredEphemeris | None
    miss: Missing | None = None

    @property
    def ephemeris(self) -> Ephemeris | None:
        """[design: convenience] the ephemeris found, or None."""
        return None if self.stored is None else self.stored.ephemeris


# ---- rules on one ephemeris ----------------------------------------------------------


def message_of(eph: Ephemeris) -> int:
    """The RTCM message an ephemeris came from, which is its type for rule E2.

    [spec: RTCM Tables 3.5-21, 3.5-22, 3.5-110, 3.5-111, 3.5-113]
    """
    return GLONASS_MESSAGE if isinstance(eph, GlonassEphemeris) else eph.message


def beidou_iodc(toc: int) -> int:
    """BeiDou D1/D2 issue of data from toc in BDT seconds of week, as SPARTN names it (rule E1).

    [spec: SPARTN Table 6.2 SF100, IODC = mod(toc/720, 240); the quotient is
    taken whole, since the result is the 8-bit field]
    """
    return (toc // 720) % 240


def ephemeris_message(sat: Sat, ephemeris_type: int) -> int | None:
    """The RTCM message whose ephemeris an OCB of this ephemeris type applies to
    (rule E2), or None when no message ssrcore takes carries that type.

    [spec: SPARTN Table 6.2 SF016 (GPS 0: L1C/A), SF017 (GLONASS 0: L1C/A),
    SF096 (Galileo 0: F/NAV, 1: I/NAV), SF097 (BeiDou 0: D1 on B1I, 1: D2 on
    B1I), SF098 (QZSS 0: LNAV); RTCM Tables 3.5-21, 3.5-22, 3.5-110, 3.5-111,
    3.5-113; BDS ICD 5.1.1, D1 from MEO/IGSO and D2 from GEO satellites.
    design: RTCM 1042 cites the B1I ICD, so the B3I types (SF097 2, 3), the
    CNAV types and every TBD value have no message here]
    """
    system, t = sat.system, ephemeris_type
    if system == System.GPS:
        return 1019 if t == 0 else None
    if system == System.GLONASS:
        return GLONASS_MESSAGE if t == 0 else None
    if system == System.GALILEO:
        return {0: 1045, 1: 1046}.get(t)
    if system == System.BEIDOU:
        geo = orbits.beidou_geo(sat.number)
        return 1042 if (t == 0 and not geo) or (t == 1 and geo) else None
    if system == System.QZSS:
        return 1044 if t == 0 else None
    return None


def reference_time(eph: Ephemeris) -> GnssTime:
    """The time an ephemeris applies from: toe, or tb for GLONASS (DESIGN 2.5).

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3 (toe); GLONASS ICD 4.4 (tb)]
    """
    return eph.tb if isinstance(eph, GlonassEphemeris) else eph.toe


def validity(eph: Ephemeris) -> float:
    """Half-width in seconds of the interval around the reference time in which
    the ephemeris may be used: the table of DESIGN 2.5.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.1, Tables 20-XII and 20-XIII for
    GPS, half the fit interval since toe is its midpoint (``gps_fit_interval``); RTCM
    DF457 for QZSS (``_QZSS_FIT``);
    design: DESIGN 2.5 for Galileo (4 h), BeiDou (2 h) and GLONASS (15 min +
    P1/2) (module constants)]
    """
    if isinstance(eph, GlonassEphemeris):
        return _GLONASS_BASE + eph.p1 / 2.0
    system = eph.sat.system
    if system == System.GPS:
        return gps_fit_interval(eph.fit_interval, eph.iodc) / 2.0
    if system == System.QZSS:
        return _QZSS_FIT / 2.0
    if system == System.GALILEO:
        return _GALILEO_HALF
    return _BEIDOU_HALF


def is_valid(eph: Ephemeris, t: GnssTime) -> bool:
    """Whether the ephemeris may be used at GPS time t: |t − reference| within ``validity``.

    [spec: DESIGN 2.5 table, see ``validity``]
    """
    return abs(t - reference_time(eph)) <= validity(eph)


def healthy(eph: Ephemeris) -> bool:
    """Rule E5: the ephemeris's health is nominal.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.1.4 (all six bits 0: LNAV data and
    signals OK); GLONASS ICD 4.4 (Bn MSB and ln 0); Galileo OS SIS ICD 2.2 5.1.9.3 (signal
    health status 0 and data validity status 0); BDS ICD 5.2.4.6 (SatH1 0). The adapter packs each system's
    flags so that 0 is nominal (``adapters.rtcm3``)]
    """
    return eph.health == 0


# ---- to the orbit kernels -----------------------------------------------------------


def own_scale(system: System, t: GnssTime) -> tuple[int, float]:
    """GPS time ``t`` as (week, seconds of week) of the scale a system's
    ephemeris reference times are counted in.

    [spec: Galileo OS SIS ICD 2.2 5.1.2 (GST) and BDS ICD 3.3 (BDT), through
    ``toolkit.timescale``; IS-GPS-200N with IRN-IS-200N-003 and -004 and RTCM DF428 for GPS and
    QZSS, whose scale is GPS time; DESIGN T4, GLONASS tb is kept on the GPS scale]
    """
    if system == System.GALILEO:
        return timescale.gps_to_gst(t.week, t.seconds)
    if system == System.BEIDOU:
        return timescale.gps_to_bdt(t.week, t.seconds)
    return t.week, t.seconds


def kepler_orbit(eph: KeplerEphemeris) -> orbits.KeplerOrbit:
    """The orbit kernel's elements of a Keplerian ephemeris.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4.3.2, Galileo OS SIS ICD 2.2 Table 68
    and BDS ICD Table 5-11: angles in semicircles are converted with the π of the ICDs
    (``orbits.semicircles``); BDS ICD 4.3 Table 4-1 for GEO
    (``orbits.beidou_geo``); toe in the system's own scale (``own_scale``),
    since the BeiDou and Galileo node terms use Ω̇e·toe of their own week]
    """
    week, toe = own_scale(eph.sat.system, eph.toe)
    s = orbits.semicircles
    return orbits.KeplerOrbit(
        system=eph.sat.system, toe_week=week, toe=toe, sqrt_a=eph.sqrt_a, e=eph.e,
        m0=s(eph.m0), delta_n=s(eph.delta_n), omega=s(eph.omega), omega0=s(eph.omega0),
        omega_dot=s(eph.omega_dot), i0=s(eph.i0), idot=s(eph.idot),
        cuc=eph.cuc, cus=eph.cus, crc=eph.crc, crs=eph.crs, cic=eph.cic, cis=eph.cis,
        geo=eph.sat.system == System.BEIDOU and orbits.beidou_geo(eph.sat.number))


def kepler_clock(eph: KeplerEphemeris) -> orbits.KeplerClock:
    """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.3.3.1; Galileo OS SIS ICD 2.2 5.1.4;
    BDS ICD 5.2.4.9] The clock polynomial, toc in the system's own scale (``own_scale``)."""
    week, toc = own_scale(eph.sat.system, eph.toc)
    return orbits.KeplerClock(toc_week=week, toc=toc, af0=eph.af0, af1=eph.af1, af2=eph.af2)


def glonass_orbit(eph: GlonassEphemeris) -> orbits.GlonassOrbit:
    """[spec: GLONASS ICD 4.4; DESIGN T4, tb on the GPS scale as the adapter
    placed it] The GLONASS kernel's state at tb."""
    return orbits.GlonassOrbit(
        tb_week=eph.tb.week, tb=eph.tb.seconds, position=eph.position,
        velocity=eph.velocity, lunisolar=eph.acceleration, tau_n=eph.tau, gamma_n=eph.gamma)


def broadcast_state(eph: Ephemeris, t: GnssTime) -> tuple[orbits.OrbitState, float]:
    """The satellite's broadcast position and velocity and its clock offset
    (s, without the relativistic term and group delays) at GPS time ``t``.

    [spec: the orbit kernels, ``toolkit.orbits`` (DESIGN 2.7); ``t`` brought
    to the ephemeris's scale by ``own_scale``. design: a kernel failure
    raises ``KernelError``, which the caller reports as a missing term]
    """
    if isinstance(eph, GlonassEphemeris):
        orbit = glonass_orbit(eph)
        return (orbits.glonass_state(orbit, t.week, t.seconds),
                orbits.glonass_clock(orbit, t.week, t.seconds))
    week, sow = own_scale(eph.sat.system, t)
    return (orbits.kepler_state(kepler_orbit(eph), week, sow),
            orbits.kepler_clock(kepler_clock(eph), week, sow))


# ---- the store -----------------------------------------------------------------------


class EphemerisStore:
    """Broadcast ephemerides by satellite and message, at most ``DEPTH`` issues each.

    [design: DESIGN 3.2; issues are kept per message as well as per satellite
    because rule E2 makes a Galileo F/NAV and an I/NAV issue with the same
    IODnav two different ephemerides]
    """

    def __init__(self, depth: int = DEPTH) -> None:
        """[design: DESIGN 3.2, three issues per satellite]"""
        if depth < 1:
            raise ValueError("depth must be at least 1")
        self._depth = depth
        self._issues: dict[tuple[Sat, int], list[StoredEphemeris]] = {}

    def add(self, eph: Ephemeris, received: GnssTime | None = None) -> FeedResult:
        """Store one ephemeris; the result says what changed.

        - An issue equal to a stored one is a duplicate; the first receipt time is kept.
        - An issue whose issue of data is stored with an older reference time
          replaces it: the IOD was reused. [spec: IS-GPS-200N with IRN-IS-200N-003 and -004
          20.3.4.4, an IODE differs from every value of the preceding six hours, so a repeated
          IOD with a newer toe is a new data set]
        - An issue whose issue of data is stored with the same reference time
          and a different body replaces it and is reported as a conflict.
          [design: the IOD names one data set (IS-GPS-200N with IRN-IS-200N-003 and -004
          20.3.3.4.1), so two bodies cannot both be right; the later is what is broadcast now, and
          the conflict is made visible]
        - An issue whose issue of data is stored with a newer reference time
          is not stored (older).
        - Otherwise the issue is added; beyond ``depth`` issues the one with the
          oldest reference time is dropped, or the new one if it is the oldest.
          [design: DESIGN 3.2, the issues kept are the current, previous and next]
        - An unhealthy issue is stored, with the reason (rule E5).
        """
        sat = eph.sat
        key = (sat, message_of(eph))
        issues = self._issues.setdefault(key, [])
        new = StoredEphemeris(eph, received)
        t_new = reference_time(eph)
        stored_reason = None if healthy(eph) else FeedReason.UNHEALTHY
        for k, old in enumerate(issues):
            if old.ephemeris.iod != eph.iod:
                continue
            if old.ephemeris == eph:
                return FeedResult(FeedOutcome.DUPLICATE, sats=(sat,), detail=f"iod {eph.iod}")
            t_old = reference_time(old.ephemeris)
            if t_new < t_old:
                return FeedResult(FeedOutcome.SUPERSEDED, FeedReason.OLDER, (sat,),
                                  detail=f"iod {eph.iod} is stored with a later reference time")
            issues[k] = new
            if t_new == t_old:
                return FeedResult(FeedOutcome.STORED, stored_reason or FeedReason.CONFLICT, (sat,),
                                  detail=f"iod {eph.iod} replaced by a different body with the same "
                                         "reference time")
            return FeedResult(FeedOutcome.STORED, stored_reason, (sat,),
                              detail=f"iod {eph.iod} reused with a later reference time")
        if len(issues) >= self._depth:
            k_oldest = min(range(len(issues)), key=lambda k: reference_time(issues[k].ephemeris))
            if t_new <= reference_time(issues[k_oldest].ephemeris):
                return FeedResult(FeedOutcome.SUPERSEDED, FeedReason.OLDER, (sat,),
                                  detail=f"{self._depth} issues with later reference times are stored")
            del issues[k_oldest]
        issues.append(new)
        issues.sort(key=lambda s: reference_time(s.ephemeris))
        return FeedResult(FeedOutcome.STORED, stored_reason, (sat,), detail=f"iod {eph.iod}")

    def find(self, sat: Sat, iod: int, ephemeris_type: int, t: GnssTime) -> EphemerisLookup:
        """The ephemeris an orbit block of satellite ``sat`` applies to at GPS time ``t``.

        Rule E1: its issue of data equals ``iod`` (SPARTN SF018, SF019, SF099,
        SF100, SF101 as sent). Rule E2: it comes from the message of the OCB
        header's ephemeris type (``ephemeris_message``). Rule E5: it is
        healthy. DESIGN 2.5: it is valid at ``t``.
        [spec: SPARTN 8.3; Table 6.2 SF016 to SF019, SF096 to SF101]
        Reasons are checked in that order, so the reason names the first rule
        that fails. [design: one reason per query, the first rule in that order]
        """
        held = [s for (s_sat, _), issues in self._issues.items() if s_sat == sat for s in issues]
        if not held:
            return EphemerisLookup(None, Missing.NO_EPHEMERIS)
        message = ephemeris_message(sat, ephemeris_type)
        if message is None:
            return EphemerisLookup(None, Missing.EPHEMERIS_TYPE)
        with_iod = [s for s in held if s.ephemeris.iod == iod]
        match = [s for s in with_iod if s.message == message]
        if not match:
            reason = Missing.EPHEMERIS_TYPE if with_iod else Missing.IODE_MISMATCH
            return EphemerisLookup(None, reason)
        found = match[0]  # one issue per (satellite, message, iod): ``add`` replaces
        if not healthy(found.ephemeris):
            return EphemerisLookup(None, Missing.UNHEALTHY)
        if not is_valid(found.ephemeris, t):
            return EphemerisLookup(None, Missing.OUTSIDE_VALIDITY)
        return EphemerisLookup(found)

    def broadcast(self, sat: Sat, t: GnssTime) -> StoredEphemeris | None:
        """The issue to use at ``t`` when there are no corrections (DESIGN 3.4 step 1):
        healthy (E5), valid at ``t`` (2.5), with its reference time nearest ``t``,
        the later reference time on a tie.

        [design: Table 20-XIII of IS-GPS-200N with IRN-IS-200N-003 and -004 puts toe at the midpoint
        of the fit, so the issue whose reference time is nearest is the best fit]
        """
        usable = [s for (s_sat, _), issues in self._issues.items() if s_sat == sat for s in issues
                  if healthy(s.ephemeris) and is_valid(s.ephemeris, t)]
        if not usable:
            return None
        return min(usable, key=lambda s: (abs(t - reference_time(s.ephemeris)),
                                          -(reference_time(s.ephemeris) - t)))

    def issues(self, sat: Sat) -> tuple[StoredEphemeris, ...]:
        """Every stored issue of a satellite, by message, oldest reference time first.

        [design: DESIGN 3.2, read access for selection and status]
        """
        return tuple(s for (s_sat, message), issues in sorted(self._issues.items(),
                                                              key=lambda kv: kv[0][1])
                     if s_sat == sat for s in issues)

    def satellites(self) -> tuple[Sat, ...]:
        """Satellites with at least one stored issue, sorted. [design: DESIGN 3.2, read access]"""
        return tuple(sorted({s for (s, _), issues in self._issues.items() if issues}))

    def glonass_channel(self, sat: Sat) -> int | None:
        """The frequency channel of a GLONASS satellite from its latest stored ephemeris, or None.

        [spec: GLONASS ICD 3.3.1.1, the channel K gives the carrier frequency;
        RTCM DF040; DESIGN 2.6.3, the channel comes from the ephemerides]
        """
        issues = self._issues.get((sat, GLONASS_MESSAGE))
        if not issues:
            return None
        return issues[-1].ephemeris.frequency_channel
