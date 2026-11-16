"""Navigation words from a receiver (UBX-RXM-SFRBX) through ``adapters.navdata`` (DESIGN Q15).

Every unit is synthetic: chosen field values, encoded with gnsswire.nav's
encoders (parity, CRC, BCH and Hamming codes computed) and carried in
UBX-RXM-SFRBX as ``gnsswire.nav.carriers.to_sfrbx`` lays it out. The
equivalent RTCM message is built independently from the same values, so
that the test states what the broadcast and the RTCM message have in common
(the field names and integers of RTCM Tables 3.5-21, 3.5-22, 3.5-111 and
3.5-113) and checks that both sources give the same ephemeris input.
"""

from __future__ import annotations

import dataclasses

import pytest
from gnsswire import rtcm3, ubx
from gnsswire.bits import layout
from gnsswire.nav import beidou, carriers, galileo, glonass, gps

from ssrcore.adapters import rtcm3 as rtcm3_adapter
from ssrcore.adapters.framing import Framer
from ssrcore.adapters.navdata import GLONASS_ROW, SFRBX, STALE_PARTS, NavData
from ssrcore.engine import Engine
from ssrcore.model.types import (
    Carried,
    FailureKind,
    FeedOutcome,
    FeedReason,
    FeedResult,
    FrameFailure,
    GlonassEphemeris,
    GnssTime,
    KeplerEphemeris,
    LeapSeconds,
    Protocol,
    Sat,
    System,
)

# GPS week 2390 (2025-10-26), Tuesday 07:33:20 GPS time.
NEAR = GnssTime(2390, 200_000.0)
LEAP = 18
SFRBX_TYPE = 0x0213   # UBX class 0x02, ID 0x13 (u-blox F9 HPG 1.51 3.17.9, X20 HPG 2.11 3.16.8)


def zero(cls, **fixed):
    """A table of ``cls`` with every field 0 (tuples of zeros) but ``fixed``."""
    values = {}
    for row in layout(cls).rows:
        n = None if row.kind.count is None else row.kind.count.eval(lambda name: 0)
        values[row.name] = 0 if n is None else (0,) * n
    values.update(fixed)
    return cls(**values)


def sfrbx(system, sat, message, **kw):
    return ubx.RxmSfrbx(**carriers.to_sfrbx(system, sat, message, **kw))


def wire(m):
    """The message as an RTCM receiver would get it: encoded and decoded by gnsswire."""
    return rtcm3.decode(rtcm3.parse(rtcm3.encode(m), 0))


def push_all(nd, messages, transport="", near=NEAR):
    out = []
    for m in messages:
        out += nd.push(m, transport, near)
    return out


def only(items, cls):
    found = [x for x in items if isinstance(x, cls)]
    assert len(found) == 1, items
    return found[0]


# ---- GPS LNAV: subframes 1 to 3, RTCM 1019 ---------------------------------------------------

GPS_1019 = dict(
    satellite=7, week=2390 % 1024, accuracy=2, l2_code=1, idot=-1234, iode=0x2A, toc=12_600,
    af2=-3, af1=-17, af0=123_456, iodc=0x12A, crs=-2345, delta_n=12_345, m0=-123_456_789,
    cuc=-1111, e=12_345_678, cus=2222, sqrt_a=2_702_000_000, toe=12_600, cic=-33,
    omega0=987_654_321, cis=44, i0=654_321_987, crc=5555, omega=-111_111_111, omega_dot=-22_222,
    tgd=-5, health=0, l2p_data=0, fit_interval=0)


def gps_subframes(e=GPS_1019, tow=33_333):
    head = dict(preamble=0x8B, anti_spoof=1)
    sf1 = zero(gps.Subframe1, **head, tow=tow, subframe=1, week=e["week"], l2_code=e["l2_code"],
               ura=e["accuracy"], health=e["health"], iodc_msb=e["iodc"] >> 8, l2p_data=e["l2p_data"],
               tgd=e["tgd"], iodc_lsb=e["iodc"] & 0xFF, toc=e["toc"], af2=e["af2"], af1=e["af1"],
               af0=e["af0"])
    sf2 = zero(gps.Subframe2, **head, tow=tow + 1, subframe=2, iode=e["iode"], crs=e["crs"],
               delta_n=e["delta_n"], m0=e["m0"], cuc=e["cuc"], e=e["e"], cus=e["cus"],
               sqrt_a=e["sqrt_a"], toe=e["toe"], fit_interval=e["fit_interval"], aodo=3)
    sf3 = zero(gps.Subframe3, **head, tow=tow + 2, subframe=3, cic=e["cic"], omega0=e["omega0"],
               cis=e["cis"], i0=e["i0"], crc=e["crc"], omega=e["omega"], omega_dot=e["omega_dot"],
               iode=e["iode"], idot=e["idot"])
    return [sfrbx(carriers.GPS, e["satellite"], gps.solve_t(sf)) for sf in (sf1, sf2, sf3)]


def gps_issue(iode, **changes):
    return {**GPS_1019, "iode": iode, "iodc": (GPS_1019["iodc"] & 0x300) | iode, **changes}


def test_gps_gives_the_ephemeris_of_the_equivalent_1019():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.4, Tables 20-I and 20-III; RTCM Table 3.5-21."""
    out = push_all(NavData(), gps_subframes())
    msg = only(out, rtcm3.Msg1019)
    expected = rtcm3.Msg1019(number=1019, **GPS_1019)
    assert msg == expected
    eph = rtcm3_adapter.ephemeris(msg, NEAR)
    assert isinstance(eph, KeplerEphemeris)
    assert eph == rtcm3_adapter.ephemeris(wire(expected), NEAR)
    assert (eph.iod, eph.iodc, eph.toe) == (0x2A, 0x12A, GnssTime(2390, 201_600.0))


def test_nothing_before_the_last_part_of_an_issue():
    nd = NavData()
    units = gps_subframes()
    assert nd.push(units[0], near=NEAR) == []
    assert nd.push(units[1], near=NEAR) == []
    assert isinstance(only(nd.push(units[2], near=NEAR), rtcm3.Msg1019), rtcm3.Msg1019)
    # The same issue again (the next frame) is not new.
    assert push_all(nd, units) == []


# ---- Galileo I/NAV: words 1 to 5, RTCM 1046 --------------------------------------------------

GAL_1046 = dict(
    satellite=11, week=2390 - 1024, iodnav=77, sisa=107, idot=-321, toc=3360, af2=-2, af1=-1234,
    af0=-123_456_789, crs=321, delta_n=4321, m0=1_234_567_890, cuc=-55, e=3_456_789, cus=66,
    sqrt_a=2_852_000_000, toe=3360, cic=-7, omega0=-1_111_111_111, cis=8, i0=699_999_999,
    crc=999, omega=222_222_222, omega_dot=-3333, bgd_e5a_e1=-12, bgd_e5b_e1=-14, e5b_health=0,
    e5b_validity=0, e1b_health=0, e1b_validity=0)


def gal_words(e=GAL_1046):
    i = e["iodnav"]
    return [
        zero(galileo.InavWord1, word_type=1, iodnav=i, toe=e["toe"], m0=e["m0"], e=e["e"],
             sqrt_a=e["sqrt_a"]),
        zero(galileo.InavWord2, word_type=2, iodnav=i, omega0=e["omega0"], i0=e["i0"],
             omega=e["omega"], idot=e["idot"]),
        zero(galileo.InavWord3, word_type=3, iodnav=i, omega_dot=e["omega_dot"],
             delta_n=e["delta_n"], cuc=e["cuc"], cus=e["cus"], crc=e["crc"], crs=e["crs"],
             sisa=e["sisa"]),
        zero(galileo.InavWord4, word_type=4, iodnav=i, svid=e["satellite"], cic=e["cic"],
             cis=e["cis"], toc=e["toc"], af0=e["af0"], af1=e["af1"], af2=e["af2"]),
        zero(galileo.InavWord5, word_type=5, bgd_e5a_e1=e["bgd_e5a_e1"], bgd_e5b_e1=e["bgd_e5b_e1"],
             e5b_hs=e["e5b_health"], e1b_hs=e["e1b_health"], e5b_dvs=e["e5b_validity"],
             e1b_dvs=e["e1b_validity"], wn=e["week"], tow=200_000),
    ]


def gal_units(words, sat=GAL_1046["satellite"], sig_id=1):
    return [sfrbx(carriers.GALILEO_INAV, sat, galileo.inav_page(w), sig_id=sig_id) for w in words]


def test_galileo_gives_the_ephemeris_of_the_equivalent_1046():
    """Galileo OS SIS ICD 2.2 4.3.5, Tables 42 to 46; RTCM Table 3.5-111."""
    out = push_all(NavData(), gal_units(gal_words()))
    msg = only(out, rtcm3.Msg1046)
    expected = rtcm3.Msg1046(number=1046, **GAL_1046, reserved=0)
    assert msg == expected
    eph = rtcm3_adapter.ephemeris(msg, NEAR)
    assert eph == rtcm3_adapter.ephemeris(wire(expected), NEAR)
    assert (eph.message, eph.iod, eph.toe) == (1046, 77, GnssTime(2390, 201_600.0))


def test_galileo_e1_and_e5b_words_of_one_satellite_join():
    """Galileo OS SIS ICD 2.2 4.3.1: E1-B and E5b-I carry the same I/NAV words."""
    words = gal_words()
    units = gal_units(words[:3], sig_id=1) + gal_units(words[3:], sig_id=5)
    assert only(push_all(NavData(), units), rtcm3.Msg1046) == rtcm3.Msg1046(number=1046, **GAL_1046,
                                                                              reserved=0)


def test_galileo_svid_contradicting_the_carrier_is_refused():
    out = push_all(NavData(), gal_units(gal_words(), sat=12))
    refused = only(out, FeedResult)
    assert (refused.outcome, refused.reason) == (FeedOutcome.REJECTED, FeedReason.INCONSISTENT)
    assert not [x for x in out if isinstance(x, rtcm3.Msg1046)]


# ---- BeiDou D1: subframes 1 to 3 of one frame, RTCM 1042 --------------------------------------

BDS_1042 = dict(
    satellite=11, week=2390 - 1356, urai=1, idot=-99, aode=5, toc=25_200, a2=-1, a1=-4321,
    a0=1_234_567, aodc=4, crs=-4444, delta_n=3333, m0=-1_212_121_212, cuc=-1234, e=2_345_678,
    cus=4321, crc=-999, sqrt_a=2_717_000_000, toe=25_200, cic=-77, omega0=1_313_131_313, cis=88,
    i0=707_070_707, omega=-1_414_141_414, omega_dot=-23_456, tgd1=17, tgd2=-23, health=0)


def bds_subframes(e=BDS_1042, sow=201_000):
    head = dict(preamble=beidou.PREAMBLE)
    sf1 = zero(beidou.D1Subframe1, **head, fraid=1, sow=sow, sat_h1=e["health"], aodc=e["aodc"],
               urai=e["urai"], week=e["week"], toc=e["toc"], tgd1=e["tgd1"], tgd2=e["tgd2"],
               a2=e["a2"], a0=e["a0"], a1=e["a1"], aode=e["aode"])
    sf2 = zero(beidou.D1Subframe2, **head, fraid=2, sow=sow + 6, delta_n=e["delta_n"], cuc=e["cuc"],
               m0=e["m0"], e=e["e"], cus=e["cus"], crc=e["crc"], crs=e["crs"], sqrt_a=e["sqrt_a"],
               toe_msb=e["toe"] >> 15)
    sf3 = zero(beidou.D1Subframe3, **head, fraid=3, sow=sow + 12, toe_lsb=e["toe"] & 0x7FFF,
               i0=e["i0"], cic=e["cic"], omega_dot=e["omega_dot"], cis=e["cis"], idot=e["idot"],
               omega0=e["omega0"], omega=e["omega"])
    return [sfrbx(carriers.BEIDOU_D1, e["satellite"], sf) for sf in (sf1, sf2, sf3)]


def test_beidou_gives_the_ephemeris_of_the_equivalent_1042():
    """BDS ICD 5.2.4, Figures 5-8 to 5-10; RTCM Table 3.5-113."""
    out = push_all(NavData(), bds_subframes())
    msg = only(out, rtcm3.Msg1042)
    expected = rtcm3.Msg1042(number=1042, **BDS_1042)
    assert msg == expected
    eph = rtcm3_adapter.ephemeris(msg, NEAR)
    assert eph == rtcm3_adapter.ephemeris(wire(expected), NEAR)
    # SPARTN SF100: IODC = mod(toc/720, 240) of the BDT toc 201 600 s.
    assert (eph.iod, eph.toe) == (40, GnssTime(2390, 201_614.0))


def test_beidou_subframes_of_two_frames_do_not_join():
    """BDS ICD 5.2.4.10, 5.2.4.12: parameters change at a superframe start, so one
    frame's subframes 1 to 3 are one issue; gnsswire ties them by SOW."""
    first, second = bds_subframes(sow=201_000), bds_subframes(BDS_1042 | {"aode": 6}, sow=201_030)
    out = push_all(NavData(), [first[0], first[1], second[2]])
    assert out == []


# ---- GLONASS: strings 1 to 5 of one frame, RTCM 1020 ------------------------------------------

GLO_1020 = dict(
    satellite=9, channel=-2 + 7, almanac_health=0, almanac_health_available=0, p1=1, tk=0x5A0,
    bn_msb=0, p2=1, tb=42, x_dot=-123_456, x=-12_345_678, x_ddot=-3, y_dot=234_567, y=23_456_789,
    y_ddot=2, z_dot=-345_678, z=9_876_543, z_ddot=-1, p3=1, gamma=-17, p=3, ln_third=0,
    tau=-98_765, delta_tau=-3, en=0, p4=0, ft=2, nt=1234, m=1, additional=1, na=1234,
    tau_c=-1_234_567, n4=8, tau_gps=-4321, ln_fifth=0, reserved=0)


def glo_strings(e=GLO_1020):
    s = [
        zero(glonass.String1, m=1, p1=e["p1"], tk=e["tk"], x_dot=e["x_dot"], x_ddot=e["x_ddot"], x=e["x"]),
        zero(glonass.String2, m=2, bn=e["bn_msb"] << 2, p2=e["p2"], tb=e["tb"], y_dot=e["y_dot"],
             y_ddot=e["y_ddot"], y=e["y"]),
        zero(glonass.String3, m=3, p3=e["p3"], gamma=e["gamma"], p=e["p"], ln=e["ln_third"],
             z_dot=e["z_dot"], z_ddot=e["z_ddot"], z=e["z"]),
        zero(glonass.String4, m=4, tau=e["tau"], delta_tau=e["delta_tau"], en=e["en"], p4=e["p4"],
             ft=e["ft"], nt=e["nt"], n=e["satellite"], sat_type=e["m"]),
        zero(glonass.String5, m=5, na=e["na"], tau_c=e["tau_c"], n4=e["n4"], tau_gps=e["tau_gps"],
             ln=e["ln_fifth"]),
    ]
    return [glonass.with_hamming(x) for x in s]


def glo_units(strings, slot=GLO_1020["satellite"], freq_id=GLO_1020["channel"]):
    return [sfrbx(carriers.GLONASS, slot, s, freq_id=freq_id) for s in strings]


def almanac_pair(na, channel, cn=0, m=6):
    h = channel + 32 if channel < 0 else channel     # GLONASS ICD Table 4.10
    return [glonass.with_hamming(zero(glonass.AlmanacEven, m=m, na=na, cn=cn)),
            glonass.with_hamming(zero(glonass.AlmanacOdd, m=m + 1, h=h))]


def test_glonass_gives_the_ephemeris_of_the_equivalent_1020():
    """GLONASS ICD 4.4, Tables 4.5, 4.6; RTCM Table 3.5-22; the channel is the carrier's
    (u-blox F9 HPG 1.51 3.17.9 and X20 HPG 2.11 3.16.8 freqId, RTCM DF040)."""
    out = push_all(NavData(), glo_units(glo_strings()))
    msg = only(out, rtcm3.Msg1020)
    expected = rtcm3.Msg1020(number=1020, **GLO_1020)
    assert msg == expected
    eph = rtcm3_adapter.ephemeris(msg, NEAR, LEAP)
    assert isinstance(eph, GlonassEphemeris)
    assert eph == rtcm3_adapter.ephemeris(wire(expected), NEAR, LEAP)
    assert (eph.frequency_channel, eph.iod) == (-2, 42)
    # tb 42 (10:30 Moscow time) of the day nearest NEAR, the RTCM adapter's rule.
    assert eph.tb - NEAR == pytest.approx(42 * 900 - (200_000 - 2 * 86_400 - LEAP + 3 * 3600))


def test_glonass_without_a_carrier_channel_takes_the_almanac():
    """GLONASS ICD 4.5 Tables 4.9 to 4.11: HnA of satellite nA, negative channels
    coded 25 to 31; a carrier that gives no channel (as SkyTraq's) is another Carrier."""
    no_channel = dataclasses.replace(SFRBX, channel=lambda m: None)
    nd = NavData((no_channel,))
    out = push_all(nd, glo_units(glo_strings()))
    refused = only(out, FeedResult)
    assert refused.reason == FeedReason.INCONSISTENT and "DF040" in refused.detail
    # Satellite 9's almanac, broadcast by another satellite on another transport: not this one's.
    push_all(nd, glo_units(almanac_pair(na=9, channel=-2, cn=1), slot=3), "other")
    assert only(push_all(nd, glo_units(glo_strings())), FeedResult).reason == FeedReason.INCONSISTENT
    # ... and on the same transport: the same issue, at its next frame, is given.
    push_all(nd, glo_units(almanac_pair(na=9, channel=-2, cn=1), slot=3))
    msg = only(push_all(nd, glo_units(glo_strings())), rtcm3.Msg1020)
    assert msg == rtcm3.Msg1020(number=1020, **GLO_1020 | {"almanac_health": 1,
                                                           "almanac_health_available": 1})
    assert rtcm3_adapter.ephemeris(msg, NEAR, LEAP).frequency_channel == -2


def test_glonass_almanac_strings_14_and_15_need_the_frame_number():
    """GLONASS ICD 4.5: in frame 5 strings 14 and 15 are not almanac; the UBX carrier
    gives the frame number only when the receiver knows it (UBX-13003221 10.3)."""
    no_channel = dataclasses.replace(SFRBX, channel=lambda m: None)
    nd = NavData((no_channel,))
    push_all(nd, glo_units(almanac_pair(na=9, channel=-2, m=14), slot=3))   # frame unknown: 0
    assert only(push_all(nd, glo_units(glo_strings())), FeedResult).reason == FeedReason.INCONSISTENT


def test_glonass_n_contradicting_the_carrier_slot_is_refused():
    out = push_all(NavData(), glo_units(glo_strings(), slot=10))
    assert only(out, FeedResult).reason == FeedReason.INCONSISTENT


def test_glonass_n_zero_takes_the_carrier_slot():
    """GLONASS ICD Table 4.5 note (1): n is a GLONASS-M word."""
    msg = only(push_all(NavData(), glo_units(glo_strings(GLO_1020 | {"satellite": 0}))), rtcm3.Msg1020)
    assert msg.satellite == 9


# ---- issue changes ---------------------------------------------------------------------------

def test_gps_issue_change_mid_set_is_never_mixed():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.4.4: IODE equals the 8 LSBs of IODC within one
    data set."""
    old, new = gps_subframes(gps_issue(0x2A)), gps_subframes(gps_issue(0x2B, toe=12_825, toc=12_825))
    nd = NavData()
    # Subframes 1 and 2 of the old issue, then the cutover: 3, 1, 2 of the new.
    out = push_all(nd, [old[0], old[1], new[2], new[0]])
    assert out == []
    msg = only(nd.push(new[1], near=NEAR), rtcm3.Msg1019)
    assert msg == rtcm3.Msg1019(number=1019, **gps_issue(0x2B, toe=12_825, toc=12_825))


def test_galileo_issue_change_mid_set_is_never_mixed():
    """Galileo OS SIS ICD 2.2 5.1.9.2: words 1 to 4 of one batch share IODnav."""
    old, new = gal_words(), gal_words(GAL_1046 | {"iodnav": 78, "toe": 3370, "toc": 3370})
    nd = NavData()
    assert push_all(nd, gal_units([old[0], old[1], new[2], new[3], old[4]])) == []
    msg = only(push_all(nd, gal_units([new[0], new[1]])), rtcm3.Msg1046)
    assert (msg.iodnav, msg.toe, msg.toc) == (78, 3370, 3370)


def test_glonass_tb_change_between_frames_is_never_mixed():
    """GLONASS ICD 4.3.2: strings 1 to 5 open each 30 s frame, 2 s each. A receiver
    passes on only the strings it decoded, so strings 1 to 3 of one frame may be
    followed directly by 4 and 5 of the next: the gap in time breaks the row."""
    old, new = glo_strings(), glo_strings(GLO_1020 | {"tb": 43, "tau": 1})
    nd = NavData()
    for k, s in enumerate(glo_units(old[:3])):
        assert nd.push(s, near=NEAR + 2.0 * k) == []
    for k, s in enumerate(glo_units(new[3:])):                    # 4, 5 of the next frame
        assert nd.push(s, near=NEAR + 30.0 + 2.0 * (k + 3)) == []
    out = []
    for k, s in enumerate(glo_units(new)):
        out += nd.push(s, near=NEAR + 60.0 + 2.0 * k)
    msg = only(out, rtcm3.Msg1020)
    assert (msg.tb, msg.tau) == (43, 1)


def test_glonass_row_within_one_frame_is_taken():
    nd = NavData()
    out = [x for k, s in enumerate(glo_units(glo_strings())) for x in nd.push(s, near=NEAR + 2.0 * k)]
    assert only(out, rtcm3.Msg1020).tb == 42
    assert GLONASS_ROW == 10.0


# ---- corrupted words ---------------------------------------------------------------------------

def corrupt(msg, word=4, bit=5):
    dwrd = list(msg.dwrd)
    dwrd[word] ^= 1 << bit
    return dataclasses.replace(msg, dwrd=tuple(dwrd))


@pytest.mark.parametrize("units, expected", [
    (gps_subframes(), rtcm3.Msg1019(number=1019, **GPS_1019)),
    (gal_units(gal_words()), rtcm3.Msg1046(number=1046, **GAL_1046, reserved=0)),
    (bds_subframes(), rtcm3.Msg1042(number=1042, **BDS_1042)),
    (glo_units(glo_strings()), rtcm3.Msg1020(number=1020, **GLO_1020)),
], ids=["gps", "galileo", "beidou", "glonass"])
def test_a_corrupted_word_is_a_decode_failure_and_joins_nothing(units, expected):
    """IS-GPS-200N with IRN-IS-200N-003 and -004 Table 20-XIV; Galileo OS SIS ICD 2.2 5.1.9.4; BDS
    ICD 5.1.3; GLONASS ICD 4.7. DESIGN 3.5: the carrier frame passed its checksum, its content
    cannot be read."""
    nd = NavData()
    bad = corrupt(units[1], word=1 if len(units[1].dwrd) == 4 else 4)   # GLONASS: word 1 is all string
    out = push_all(nd, [units[0], bad] + units[2:], near=NEAR)
    failure = only(out, FrameFailure)
    assert failure == FrameFailure(FailureKind.DECODE, Protocol.UBX, SFRBX_TYPE, None, None)
    assert not [x for x in out if not isinstance(x, FrameFailure)]
    # The good unit completes the issue (GLONASS: a new row, its strings in order).
    rest = units if len(units) == 5 else [units[1]]
    assert only(push_all(nd, rest), type(expected)) == expected


def test_beidou_all_zero_parity_fields_are_taken_unchecked():
    """gnsswire.nav: receivers send the BeiDou parity fields as zero (UBX-13003221 10.4)."""
    units = [dataclasses.replace(u, dwrd=tuple(w & ~0xFF if i else w & ~0xF for i, w in enumerate(u.dwrd)))
             for u in bds_subframes()]
    assert only(push_all(NavData(), units), rtcm3.Msg1042) == rtcm3.Msg1042(number=1042, **BDS_1042)


def test_failure_carries_the_receipt_time():
    t = GnssTime(2390, 1.0)
    out = NavData().push(corrupt(gps_subframes()[0]), "a", NEAR, received=t)
    assert out == [FrameFailure(FailureKind.DECODE, Protocol.UBX, SFRBX_TYPE, None, t)]


# ---- transports, replay, memory ----------------------------------------------------------------

def test_two_transports_do_not_mix_words():
    """DESIGN 2.2: per named transport, as the Framer keeps its buffers."""
    a, b = gps_subframes(gps_issue(0x2A)), gps_subframes(gps_issue(0x2A, af0=7))
    nd = NavData()
    interleaved = [(a[0], "rx1"), (b[1], "rx2"), (a[1], "rx1"), (b[2], "rx2"), (b[0], "rx2")]
    out = {name: [] for name in ("rx1", "rx2")}
    for m, name in interleaved:
        out[name] += nd.push(m, name, NEAR)
    assert out["rx1"] == []
    assert only(out["rx2"], rtcm3.Msg1019).af0 == 7
    assert only(nd.push(a[2], "rx1", NEAR), rtcm3.Msg1019).af0 == GPS_1019["af0"]
    assert nd.transports == ("rx1", "rx2")
    nd.end("rx1")
    assert nd.transports == ("rx2",)


def scenario():
    gp, ga, bd, gl = gps_subframes(), gal_units(gal_words()), bds_subframes(), glo_units(glo_strings())
    page18 = sfrbx(carriers.GPS, 7, gps.solve_t(zero(gps.Subframe4Page18, preamble=0x8B, subframe=4,
                                                      data_id=1, sv_id=56, dt_ls=18, dt_lsf=18)))
    seq = [gp[0], ga[0], gl[0], bd[0], gp[1], ga[1], gl[1], corrupt(bd[1]), ga[2], gl[2], gp[2],
           page18, ga[3], gl[3], bd[1], bd[2], ga[4], gl[4]] + bds_subframes(sow=201_030)
    # Receiver 1 tracks GPS and GLONASS, receiver 2 Galileo and BeiDou; GLONASS strings
    # also reach receiver 2, alone, so that they join nothing there.
    receiver = {0: "rx1", 6: "rx1", 2: "rx2", 3: "rx2"}
    return [(m, receiver[m.gnss_id]) for m in seq] + [(gl[0], "rx2"), (gl[1], "rx2")]


def run(seq):
    nd = NavData()
    return [nd.push(m, name, NEAR + i) for i, (m, name) in enumerate(seq)]


def test_replay_gives_the_same_outputs():
    """DESIGN 2.1, I/O at the edge: every run can be replayed from recorded inputs."""
    seq = scenario()
    first = run(seq)
    assert first == run(seq)
    kinds = sorted(type(x).__name__ for out in first for x in out)
    # The BeiDou issue joins at the next frame, after the corrupted subframe 2.
    assert kinds == ["FrameFailure", "LeapSeconds", "Msg1019", "Msg1020", "Msg1042", "Msg1046"]


def test_assembly_state_is_bounded_per_satellite():
    nd = NavData()
    for k in range(50):
        push_all(nd, gps_subframes(gps_issue(k), tow=k * 5), "rx", NEAR + 30.0 * k)
    state = nd._transports["rx"]
    assert len(state.assemblies) == len(state.seen) == 1


def test_parts_older_than_the_stale_limit_are_dropped():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.4.4: an IODE may come back after six hours, so
    an old part must not meet a new issue of the same IODE (design: STALE_PARTS)."""
    old, new = gps_subframes(gps_issue(0x2A)), gps_subframes(gps_issue(0x2A, af0=99, toe=13_725,
                                                                         toc=13_725))
    nd = NavData()
    push_all(nd, old[1:])                                    # 2 and 3 of the old issue
    later = NEAR + (STALE_PARTS + 1.0)
    assert nd.push(new[0], near=later) == []                  # 1 of the new: the old parts are gone
    assert only(push_all(nd, new[1:], near=later), rtcm3.Msg1019).af0 == 99


def test_a_stale_part_is_dropped_while_other_parts_arrive():
    """design: STALE_PARTS ages each kind of part, not the satellite's last unit."""
    old, new = gps_subframes(gps_issue(0x2A)), gps_subframes(gps_issue(0x2A, af0=99, toe=13_725,
                                                                         toc=13_725))
    nd = NavData()
    assert nd.push(old[1], near=NEAR) == []                            # subframe 2, then lost
    assert nd.push(new[0], near=NEAR + 250.0) == []
    assert nd.push(new[2], near=NEAR + (STALE_PARTS + 1.0)) == []      # not joined to the old 2
    # The assembly started over: the next frame's subframes 1 and 2 complete it.
    assert nd.push(new[0], near=NEAR + (STALE_PARTS + 30.0)) == []
    assert only(nd.push(new[1], near=NEAR + (STALE_PARTS + 31.0)), rtcm3.Msg1019).af0 == 99


def test_without_a_time_the_parts_are_kept():
    nd = NavData()
    out = push_all(nd, gps_subframes(), near=None)
    assert only(out, rtcm3.Msg1019) == rtcm3.Msg1019(number=1019, **GPS_1019)


# ---- leap seconds ------------------------------------------------------------------------------

def gps_page18(**f):
    return sfrbx(carriers.GPS, 7, gps.solve_t(zero(gps.Subframe4Page18, preamble=0x8B, subframe=4,
                                                    data_id=1, sv_id=56, **f)))


def test_gps_page_18_gives_the_leap_seconds_once():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4, its NOTE: with ΔtLS = ΔtLSF no
    effectivity time is needed."""
    nd = NavData()
    assert nd.push(gps_page18(dt_ls=18, dt_lsf=18, wn_lsf=0, dn=0)) == [LeapSeconds(18)]
    assert nd.push(gps_page18(dt_ls=18, dt_lsf=18, wn_lsf=0, dn=0)) == []
    assert nd.push(gps_page18(dt_ls=18, dt_lsf=18), "other") == [LeapSeconds(18)]


def test_gps_announced_leap_second():
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4: WN_LSF modulo 256 and DN, day one the
    first of the week; the new count applies from the end of day DN (as adapters.ubx)."""
    nd = NavData()
    leap = only(nd.push(gps_page18(dt_ls=18, dt_lsf=19, wn_lsf=(2391) % 256, dn=7), near=NEAR), LeapSeconds)
    assert leap == LeapSeconds(18, 19, GnssTime(2392, 19.0))
    # Without a time the week cannot be resolved: the current count alone.
    assert NavData().push(gps_page18(dt_ls=18, dt_lsf=19, wn_lsf=2391 % 256, dn=7)) == [LeapSeconds(18)]


@pytest.mark.parametrize("data_id", [0b00, 0b10, 0b11])
def test_gps_page_18_of_another_data_id_is_not_read(data_id):
    """IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.1.1 and 6.4.6.4: IRN-IS-200N-003
    allows future data structures with a Data ID other than 2 (binary 01), and the user
    determines the Data ID before using LNAV data; such a page gives no leap seconds."""
    nd = NavData()
    page = sfrbx(carriers.GPS, 7, gps.solve_t(zero(gps.Subframe4Page18, preamble=0x8B, subframe=4,
                                                    data_id=data_id, sv_id=56, dt_ls=18, dt_lsf=18)))
    assert nd.push(page, near=NEAR) == []
    assert nd.push(gps_page18(dt_ls=18, dt_lsf=18), near=NEAR) == [LeapSeconds(18)]


def test_galileo_word_6_gives_the_leap_seconds():
    """Galileo OS SIS ICD 2.2 5.1.7 Table 75 (DN 1 = Sunday), 5.1.2 (GST − UTC is GPS − UTC)."""
    word = zero(galileo.InavWord6, word_type=6, dt_ls=18, dt_lsf=19, wn_lsf=(2391 - 1024) % 256, dn=1)
    leap = only(NavData().push(gal_units([word])[0], near=NEAR), LeapSeconds)
    assert leap == LeapSeconds(18, 19, GnssTime(2391, 86_400.0 + 19.0))


def test_beidou_utc_page_gives_the_leap_seconds():
    """BDS ICD 5.2.4.18 (BDT − UTC, DN 0 = Sunday), 3.3 (BDT = GPS − 14 s)."""
    page = zero(beidou.Utc, preamble=beidou.PREAMBLE, fraid=5, pnum=10, dt_ls=4, dt_lsf=5,
                wn_lsf=(2391 - 1356) % 256, dn=0)
    leap = only(NavData().push(sfrbx(carriers.BEIDOU_D1, 11, page), near=NEAR), LeapSeconds)
    assert leap == LeapSeconds(18, 19, GnssTime(2391, 86_400.0 + 19.0))


# ---- what is not taken -----------------------------------------------------------------------------

def test_a_signal_without_a_layout_is_refused():
    """gnsswire.nav: QZSS and GPS L2C words raise Unsupported (no staged layout)."""
    qzss = dataclasses.replace(gps_subframes()[0], gnss_id=5)
    out = NavData().push(qzss)
    assert out == [FeedResult(FeedOutcome.REJECTED, FeedReason.UNSUPPORTED_SYSTEM, detail=out[0].detail)]


@pytest.mark.parametrize("gnss_id, sig_id", [
    (2, 8), (2, 9), (2, 10),          # Galileo E6 B, E6 C, E6 A
    (3, 4), (3, 10),                  # BeiDou B3I D1, B3I D2
    (3, 5), (3, 6), (3, 7), (3, 8),   # BeiDou B1 Cp, B1 Cd, B2 ap, B2 ad
    (5, 12),                          # QZSS L1C/B (X20 only)
])
def test_the_signals_f9_hpg_1_51_and_the_x20_add_are_refused(gnss_id, sig_id):
    """u-blox F9 HPG 1.51 and X20 HPG 2.11 1.5.4 (the sigIds they add, and B1C and B2a
    split into pilot, 5 and 7, and data, 6 and 8; QZSS L1C/B is the X20's alone) and
    F9 HPG 1.51 3.17.9.1, X20 HPG 2.11 3.16.8.1 (RXM-SFRBX of up to 16 words):
    gnsswire-nav lays none of these out, so each is unsupported_system, and the GPS words
    that follow are taken as before."""
    unit = dataclasses.replace(gps_subframes()[0], gnss_id=gnss_id, sv_id=11, sig_id=sig_id,
                               num_words=16, dwrd=(0,) * 16)
    nd = NavData()
    out = nd.push(unit, near=NEAR)
    assert out == [FeedResult(FeedOutcome.REJECTED, FeedReason.UNSUPPORTED_SYSTEM, detail=out[0].detail)]
    assert only(push_all(nd, gps_subframes()), rtcm3.Msg1019) == rtcm3.Msg1019(number=1019, **GPS_1019)


def test_a_message_no_carrier_takes_is_a_type_error():
    with pytest.raises(TypeError):
        NavData().push(rtcm3.Msg1019(number=1019, **GPS_1019))


# ---- through the Framer and the engine (DESIGN 2.2, Q15) -------------------------------------------

def test_the_framer_gives_navigation_words_with_their_transport():
    """DESIGN 2.2: RXM-SFRBX leaves the framer as ``Carried``, other messages as they are."""
    unit = gps_subframes()[0]
    wire_bytes = ubx.encode(unit)
    framer = Framer()
    out = framer.push("rx1", wire_bytes[:9]) + framer.push("rx1", wire_bytes[9:])
    assert out == [Carried("rx1", unit)]
    msg1019 = rtcm3.Msg1019(number=1019, **GPS_1019)
    assert framer.push("rx2", rtcm3.encode(msg1019)) == [wire(msg1019)]


def test_the_engine_stores_the_ephemeris_the_words_give_as_the_rtcm_message_would():
    """DESIGN 2.2 and Q15: the words of one issue give the ephemeris of the equivalent 1019,
    stored by the same code (rules E1 to E5)."""
    words = Engine()
    results = [words.feed(Carried("rx1", m), NEAR) for m in gps_subframes()]
    assert [r.detail for r in results[:2]] == ["navigation unit taken into its assembly"] * 2
    assert results[2].outcome == FeedOutcome.STORED
    rtcm = Engine()
    rtcm.feed(wire(rtcm3.Msg1019(number=1019, **GPS_1019)), NEAR)
    sat = Sat(System.GPS, GPS_1019["satellite"])
    assert words.ephemerides.satellites() == rtcm.ephemerides.satellites() == (sat,)
    assert words.query(sat, NEAR, (6378137.0, 0.0, 0.0)).satellite == \
        rtcm.query(sat, NEAR, (6378137.0, 0.0, 0.0)).satellite


def test_the_engine_keeps_the_words_of_two_transports_apart():
    """DESIGN 2.2: subframe 1 of one receiver does not complete the issue of another."""
    e = Engine()
    sf1, sf2, sf3 = gps_subframes()
    e.feed(Carried("rx1", sf1), NEAR)
    e.feed(Carried("rx2", sf2), NEAR)
    r = e.feed(Carried("rx2", sf3), NEAR)
    assert r.detail == "navigation unit taken into its assembly" and not e.ephemerides.satellites()
    e.feed(sf1, NEAR)          # a bare RXM-SFRBX is of the unnamed transport ""
    assert not e.ephemerides.satellites()
    e.feed(Carried("rx2", sf1), NEAR)
    assert e.ephemerides.satellites() == (Sat(System.GPS, GPS_1019["satellite"]),)


def test_the_engine_holds_an_assembled_ephemeris_until_it_knows_the_time():
    """Rule T5: the assembly needs no time; the conversion waits in the engine's hold."""
    e = Engine()
    results = [e.feed(Carried("rx1", m)) for m in gps_subframes()]
    assert (results[2].outcome, results[2].reason) == (FeedOutcome.HELD, FeedReason.NO_TIME)
    r = e.feed(LeapSeconds(LEAP), NEAR)
    assert "released" in r.detail and e.ephemerides.satellites()


def test_the_engine_counts_a_word_that_fails_its_code_and_takes_broadcast_leap_seconds():
    """DESIGN 3.5 (rule S7) for the failure; rule T4 for GPS subframe 4 page 18."""
    e = Engine()
    r = e.feed(Carried("rx1", corrupt(gps_subframes()[0])), NEAR)
    assert (r.outcome, r.reason) == (FeedOutcome.COUNTED, FeedReason.FRAME_FAILURE)
    r = e.feed(Carried("rx1", gps_page18(dt_ls=18, dt_lsf=18, wn_lsf=0, dn=0)), NEAR)
    assert r.outcome == FeedOutcome.STORED and e.timebase.leap_seconds == LeapSeconds(18)
    with pytest.raises(TypeError):
        e.feed(Carried("rx1", rtcm3.Msg1019(number=1019, **GPS_1019)))
