"""gnsswire.ubx messages: NAV-TIMELS to leap seconds; RXM-PMP user data to SPARTN bytes.

DESIGN.md section 2.2. gnsswire gives the raw integers of each payload; their
meaning is the u-blox F9 HPG 1.51 Interface Description (NAV-TIMELS 3.15.26,
RXM-PMP version 1 3.17.3.1), the u-blox X20 HPG 2.11 Interface Description,
which defines both alike (NAV-TIMELS 3.14.22, RXM-PMP version 1 3.16.3.1), and
the u-blox D9 PMP 1.04 Interface Description (RXM-PMP versions 0 and 1,
2.13.1.1 and 2.13.1.2).
"""

from __future__ import annotations

from gnsswire import ubx

from ssrcore.model.types import GnssTime, LeapSeconds

# [spec: u-blox F9 HPG 1.51 3.15.26.1 and u-blox X20 HPG 2.11 3.14.22.1, srcOfCurrLs]
# 0 is "Default (hardcoded in the firmware, can be outdated)" and 255 "Unknown".
# [design: DESIGN T4, a leap second count that may be wrong is not used, since
# a 1 s error in GLONASS time would pass silently]
UNTRUSTED_LEAP_SOURCES = frozenset({0, 255})

# [spec: u-blox F9 HPG 1.51 3.15.26.1 and u-blox X20 HPG 2.11 3.14.22.1,
# srcOfLsChange 4 = BeiDou; dateOfLsGpsDn "GPS and Galileo DN: from 1 = Sun to
# 7 = Sat. BeiDou DN: from 0 = Sun to 6 = Sat"]
_BEIDOU_SOURCE = 4

# [spec: u-blox F9 HPG 1.51 3.15.26.1 and u-blox X20 HPG 2.11 3.14.22.1, valid bit 0
# validCurrLs, bit 1 validTimeToLsEvent]
_VALID_CURRENT = 0x01
_VALID_EVENT = 0x02

_DAY = 86400.0


def leap_seconds(message: ubx.NavTimels) -> LeapSeconds | None:
    """GPS − UTC and the announced change of a NAV-TIMELS message; None when the
    current count is not valid or comes from an untrusted source.

    ``effective``, the GPS time from which the new count applies, is the
    first instant of the UTC day after day DN of week WN: the end of GPS day
    DN, (DN · 86400 s into week WN), plus the new count, since GPS − UTC is
    the new count from that instant on.
    A zero lsChange gives no event: "the event should be ignored" (F9 HPG
    1.51 and X20).
    [spec: u-blox F9 HPG 1.51 3.15.26.1 and u-blox X20 HPG 2.11 3.14.22.1
    (currLs, lsChange, dateOfLsGpsWn, dateOfLsGpsDn, valid); IS-GPS-200N with
    IRN-IS-200N-003 and -004 20.3.3.5.2.4 (the leap second becomes effective near
    the end of day DN of week WN_LSF; "day one" is the first day of the week)]
    """
    if not message.valid & _VALID_CURRENT or message.src_of_curr_ls in UNTRUSTED_LEAP_SOURCES:
        return None
    current = message.curr_ls
    if not message.valid & _VALID_EVENT or message.ls_change == 0:
        return LeapSeconds(current)
    if message.ls_change not in (-1, 1):
        raise ValueError(f"NAV-TIMELS lsChange {message.ls_change} is not -1, 0 or +1")
    day = message.date_of_ls_gps_dn + (1 if message.src_of_ls_change == _BEIDOU_SOURCE else 0)
    if not 1 <= day <= 7:
        raise ValueError(f"NAV-TIMELS dateOfLsGpsDn {message.date_of_ls_gps_dn} is not a day of week")
    after = current + message.ls_change
    # [design: dateOfLsGpsWn is taken as the full week number, as its U2 type and
    #  the description "GPS week number (WN)" give it, without the modulo 256 of
    #  the broadcast WN_LSF; F9 HPG 1.51 3.15.26.1 and X20 HPG 2.11 3.14.22.1
    #  resolve the week of an upcoming event and warn of a 256-week ambiguity only
    #  for a past one more than 256 weeks back, whose change is no longer announced]
    effective = GnssTime.normalised(message.date_of_ls_gps_wn, day * _DAY + after)
    return LeapSeconds(current, after, effective)


def pmp_user_data(message: ubx.RxmPmp | ubx.RxmPmpV0) -> bytes:
    """The user data of an RXM-PMP message: bytes of the L-band SPARTN stream.

    Version 1 carries numBytesUserData bytes; version 0 a fixed block of 504
    with no count, all of which are returned (bytes that belong to no frame
    are skipped by the framing). The data is not a whole number of SPARTN
    frames: "the validity of the frame must be verified by host software",
    which the SPARTN CRCs do when the bytes are framed.
    [spec: u-blox D9 PMP 1.04 2.13.1.1 and 2.13.1.2; u-blox F9 HPG 1.51
    3.17.3.1; u-blox X20 HPG 2.11 3.16.3.1; SPARTN Table 7.1 TF006, TF018]
    """
    return bytes(message.user_data)
