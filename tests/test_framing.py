"""The framing helper (DESIGN 2.2): interleaved transports, frames split across
pieces, bytes that belong to no frame, corrupted frames and frames that fail to
decode. Every frame is synthetic, encoded with gnsswire, or a worked example of
the u-blox F9 HPG 1.51 Interface Description, which the u-blox X20 HPG 2.11
one repeats.
"""

from __future__ import annotations

import random

import pytest
from gnsswire import nmea, rtcm3, spartn, ubx

from ssrcore.adapters.framing import Framer
from ssrcore.model.types import (FailureKind, FrameFailure, GnssTime, OcbHeaderInput, Protocol, SpartnTag,
                                 StreamId)


def description(subtype=0, tag=500_000_000, eaf=0):
    enc = dict(encryption_id=1, sequence=3, auth_indicator=0, auth_length=0) if eaf else {}
    return spartn.Description(subtype=subtype, time_tag_type=1, time_tag=tag, solution_id=5,
                              processor_id=1, **{"encryption_id": None, "sequence": None,
                                                 "auth_indicator": None, "auth_length": None,
                                                 **enc})


def gad(siou=3, tag=500_000_000) -> spartn.Message:
    area = spartn.GadArea(area_id=7, latitude=1300, longitude=1850, latitude_nodes=4,
                          longitude_nodes=5, latitude_spacing=4, longitude_spacing=9)
    body = spartn.Gad(siou=siou, aiou=2, reserved=0, area_count=0, areas=(area,))
    return spartn.Message(2, 0, 2, description(tag=tag), body)


def gad_bytes(**kw) -> bytes:
    return spartn.encode(gad(**kw))


def decoded_gad(**kw):
    return spartn.decode(next(spartn.scan(gad_bytes(**kw))))


# RTCM 10403.3 4.2: the 1005 message of the worked example (reference station 2003,
# ARP 1114104.5999, -4850729.7108, 3975521.4643 m).
MSG1005 = rtcm3.Msg1005(number=1005, station=2003, itrf_year=0, gps=1, glonass=0, galileo=0,
                        reference_station=0, x=11141045999, single_oscillator=0, reserved=0,
                        y=-48507297108, quarter_cycle=0, z=39755214643)

ZDA = b"$GPZDA,082710.00,16,09,2002,00,00*64\r\n"            # u-blox F9 2.7.21.1, X20 2.7.20.1
GGA = b"$GPGGA,092725.00,4717.11399,N,00833.91590,E,1,08,1.01,499.6,M,48.0,M,,*5B\r\n"  # 2.7.5.1


def timels() -> bytes:
    return ubx.encode(ubx.NavTimels(itow=1000, version=0, reserved0=(0, 0, 0), src_of_curr_ls=2,
                                    curr_ls=18, src_of_ls_change=0, ls_change=0,
                                    time_to_ls_event=0, date_of_ls_gps_wn=0, date_of_ls_gps_dn=0,
                                    reserved1=(0, 0, 0), valid=1))


def pmp(user: bytes) -> bytes:
    return ubx.encode(ubx.RxmPmp(version=1, reserved0=0, num_bytes_user_data=len(user),
                                 time_tag=0, unique_word=(0, 0), service_identifier=0, spare=0,
                                 unique_word_bit_errors=0, fec_bits=0, ebno=80, reserved1=0,
                                 user_data=user))


def pieces(data: bytes, rng: random.Random, largest: int = 9) -> list[bytes]:
    out, i = [], 0
    while i < len(data):
        n = rng.randint(1, largest)
        out.append(data[i:i + n])
        i += n
    return out


# ---- one transport ------------------------------------------------------------------------


def test_each_codec_decodes():
    f = Framer()
    out = f.push("port", gad_bytes() + rtcm3.encode(MSG1005) + timels() + ZDA + GGA)
    assert [type(m).__name__ for m in out] == ["Message", "Msg1005", "NavTimels", "Zda", "Gga"]
    assert out[0] == decoded_gad()
    assert out[1] == MSG1005
    assert f.pending("port") == 0


def test_garbage_between_frames_is_skipped():
    """Bytes that belong to no frame, including the start bytes of every codec."""
    junk = bytes([0x73, 0x00, 0xD3, 0xFF, 0xB5, 0x00, 0x24, 0x2C, 0x55, 0xB5, 0x62])
    data = junk + gad_bytes() + b"noise$$**" + ZDA + junk[:6] + rtcm3.encode(MSG1005)
    out = Framer().push_all([("ip", data)])
    assert out == [decoded_gad(), nmea.decode(next(nmea.scan(ZDA))), MSG1005]


def test_a_tag_block_before_a_sentence_is_skipped():
    """u-blox F9 HPG 1.51 and X20 HPG 2.11 2.8.1.1: the secondary output GGA, its worked
    example, is the sentence of 2.7.5.1 behind an NMEA TAG block, whose bytes belong to no
    frame."""
    secondary = b"\\s:1*78\\" + GGA
    assert Framer().push("port", secondary + ZDA) == Framer().push("port", GGA + ZDA)
    assert [type(m).__name__ for m in Framer().push("port", secondary)] == ["Gga"]


def test_corrupted_frames_give_nothing_and_the_next_frame_is_found():
    """DESIGN 3.5: a CRC failure reaches the engine as nothing."""
    for frame in (gad_bytes(), rtcm3.encode(MSG1005), timels(), ZDA):
        for at in (len(frame) // 2, len(frame) - 3):
            bad = bytearray(frame)
            bad[at] ^= 0x10
            out = Framer().push_all([("x", bytes(bad) + gad_bytes(siou=9))])
            assert out == [decoded_gad(siou=9)], (frame[:2], at)


@pytest.mark.parametrize("seed", range(5))
def test_frames_split_across_pieces(seed):
    rng = random.Random(seed)
    data = b"\x00\x73" + gad_bytes() + timels() + GGA + rtcm3.encode(MSG1005) + gad_bytes(siou=4)
    whole = Framer().push_all([("t", data)])
    f = Framer()
    split = []
    for piece in pieces(data, rng):
        split += f.push("t", piece)
    assert split + f.flush("t") == whole
    assert len(whole) == 5


def test_frame_cut_off_by_the_end_gives_nothing():
    """DESIGN 3.5: bytes that stop inside a frame give nothing; flush clears them."""
    f = Framer()
    frame = gad_bytes()
    assert f.push("t", frame[:-4]) == []
    assert f.pending("t") == len(frame) - 4
    assert f.flush("t") == []
    assert f.pending("t") == 0
    assert f.push("t", frame) == [decoded_gad()]


def test_false_start_waits_then_resolves():
    """A cut-off frame followed by whole ones: the whole ones come out once the false
    start is shown invalid by later bytes, or at the end of the transport."""
    big = spartn.encode_frame(0, 0, 2, description(subtype=9), bytes(300))
    cut = big[:12]                          # a valid frame start whose end never comes
    tail = ZDA + gad_bytes(siou=2)          # whole frames, shorter than the cut-off one claims
    expected = [nmea.decode(next(nmea.scan(ZDA))), decoded_gad(siou=2)]
    f = Framer()
    assert f.push("t", cut + tail) == []    # still possibly the cut-off frame
    assert f.pending("t") == len(cut + tail)
    assert f.flush("t") == expected         # the end of the transport gives it up
    # Without an end: enough later bytes show the start invalid, and the frames come out.
    f = Framer()
    assert f.push("t", cut + tail) == []
    assert f.push("t", bytes(len(big))) == expected


def test_a_transport_of_one_protocol_is_not_held_by_another_protocols_start():
    """DESIGN 2.2, *Delay*: in a SPARTN stream, a stray RTCM preamble whose length
    field claims more than follows holds the next SPARTN frame when every protocol
    is looked for; with the transport declared SPARTN-only it is not a start."""
    stray = bytes([rtcm3.PREAMBLE, 0x00, 0xFF]) + bytes(40)   # claims 255 + 6 bytes
    data = stray + gad_bytes(siou=5)
    assert Framer().push("ip", data) == []
    f = Framer(protocols={"ip": (Protocol.SPARTN,)})
    assert f.push("ip", data) == [decoded_gad(siou=5)]
    assert f.pending("ip") == 0
    # a transport it does not name still carries every protocol
    assert f.push("port", rtcm3.encode(MSG1005)) == [MSG1005]
    with pytest.raises(ValueError, match="no protocol"):
        Framer(protocols={"ip": ()})


# ---- interleaved transports -------------------------------------------------------------


@pytest.mark.parametrize("seed", range(5))
def test_interleaved_transports_keep_their_own_buffers(seed):
    """An IP correction stream, an L-band stream inside RXM-PMP and a receiver port,
    each cut into pieces and interleaved; each transport's frames come out whole and in
    order, and the messages of the three appear in the order they were completed."""
    rng = random.Random(100 + seed)
    ip = gad_bytes(siou=1) + rtcm3.encode(MSG1005) + gad_bytes(siou=2)
    lband_spartn = gad_bytes(siou=11) + gad_bytes(siou=12)
    lband = b"".join(pmp(p) for p in pieces(lband_spartn, rng, largest=40))
    port = ZDA + timels() + GGA + b"$GPGSV,1,1,00*79\r\n"

    streams = {"ip": pieces(ip, rng), "lband": pieces(lband, rng), "port": pieces(port, rng)}
    order = [name for name, ps in streams.items() for _ in ps]
    rng.shuffle(order)
    f = Framer()
    out_by_transport = {name: [] for name in streams}
    sequence = []
    for name in order:
        for m in f.push(name, streams[name].pop(0)):
            out_by_transport[name].append(m)
            sequence.append(m)
    for name in streams:
        assert f.flush(name) == []

    assert out_by_transport["ip"] == [decoded_gad(siou=1), MSG1005, decoded_gad(siou=2)]
    assert out_by_transport["lband"] == [decoded_gad(siou=11), decoded_gad(siou=12)]
    assert [type(m).__name__ for m in out_by_transport["port"]] == ["Zda", "NavTimels", "Gga"]
    assert len(sequence) == 8
    assert set(f.transports) <= {"ip", "lband", "port"}


def test_pmp_user_data_is_a_spartn_stream():
    """u-blox D9 PMP 1.04 2.13.1: the user data of successive RXM-PMP messages carry
    the SPARTN frames, cut anywhere; only SPARTN is looked for in it."""
    frame = gad_bytes(siou=6)
    f = Framer()
    assert f.push("lband", pmp(b"$GPZDA" + frame[:10])) == []
    assert f.push("lband", pmp(frame[10:] + ZDA)) == [decoded_gad(siou=6)]
    # Version 0: a fixed 504-byte block, all of it framed.
    v0 = ubx.RxmPmpV0(version=0, reserved0=(0, 0, 0), time_tag=0, unique_word=(0, 0),
                      service_identifier=0, spare=0, unique_word_bit_errors=0,
                      user_data=frame + bytes(504 - len(frame)), fec_bits=0, ebno=0, reserved1=0)
    assert Framer().push("lband", ubx.encode(v0)) == [decoded_gad(siou=6)]


# ---- frames that fail to decode ----------------------------------------------------------


def test_spartn_failures_are_reported():
    received = GnssTime(2300, 10.0)
    unsupported = spartn.encode_frame(1, 0, 2, description(subtype=9), bytes(20))
    undecodable = spartn.encode_frame(2, 0, 2, description(subtype=0), bytes(1))
    encrypted = spartn.encode_frame(2, 1, 2, description(subtype=0, eaf=1), bytes(30))
    short_header = spartn.encode_frame(0, 0, 2, description(subtype=9), bytes(1))
    out = Framer().push("ip", unsupported + undecodable + encrypted + short_header, received)
    assert out == [FrameFailure(FailureKind.UNSUPPORTED, Protocol.SPARTN, 1, 9, received),
                   FrameFailure(FailureKind.DECODE, Protocol.SPARTN, 2, 0, received),
                   FrameFailure(FailureKind.ENCRYPTED, Protocol.SPARTN, 2, 0, received),
                   FrameFailure(FailureKind.UNSUPPORTED, Protocol.SPARTN, 0, 9, received)]


def test_ocb_header_of_a_subtype_without_a_table_is_read():
    """SPARTN 8.4 item 5: SF010 is read in every SM 0 message (DESIGN Q19)."""
    header = spartn.OcbHeader(siou=5, end_of_set=1, reserved=0, yaw_present=0, datum=1,
                              ephemeris_type=None, mask_size=None, mask=None)
    from gnsswire.spartn import pack
    payload = pack(header, subtype=9) + bytes(8)
    frame = spartn.encode_frame(0, 0, 2, description(subtype=9), payload)
    assert Framer().push("ip", frame) == [
        OcbHeaderInput(9, StreamId(5, 1), SpartnTag(500_000_000, True), 5, True, 1)]


def test_data_after_a_table_is_kept_not_refused():
    """SPARTN 8.4 item 4 and RTCM 3.5.17: appended data reaches the engine in the
    message's ``extension`` (DESIGN Q20)."""
    longer_gad = spartn.encode_frame(2, 0, 2, description(), spartn.parse(gad_bytes()).payload + b"\xab")
    longer_1005 = rtcm3.encode_frame(rtcm3.parse(rtcm3.encode(MSG1005)).payload + b"\x12\x34")
    gad_out, rtcm_out = Framer().push("ip", longer_gad + longer_1005)
    assert gad_out.body.extension is not None and gad_out.body.areas == decoded_gad().body.areas
    assert rtcm_out.extension.bits == 16 and rtcm_out.station == 2003


def test_rtcm_failures_are_reported():
    no_table = rtcm3.encode_frame((1043 << 4).to_bytes(2, "big") + bytes(20))
    short = rtcm3.encode_frame((1005 << 4).to_bytes(2, "big") + bytes(3))
    assert Framer().push("ntrip", no_table + short) == [
        FrameFailure(FailureKind.UNSUPPORTED, Protocol.RTCM3, 1043),
        FrameFailure(FailureKind.DECODE, Protocol.RTCM3, 1005)]


def test_every_msm_level_and_1044_decode():
    """gnsswire decodes MSM1 to MSM7 and 1044; the framer passes them on as messages."""
    msm7 = rtcm3.encode_frame((1077 << 4).to_bytes(2, "big") + bytes(20))
    [out] = Framer().push("ntrip", msm7)
    assert isinstance(out, rtcm3.Msm) and out.number == 1077


def test_receiver_messages_without_a_table_are_skipped():
    """NAV-PVT and GSV have no gnsswire table and are not inputs: the framer counts
    them and passes nothing on (DESIGN 3.5); a known message that does not fit its
    table is a failure of its protocol, a UBX one with its class and ID. A
    sentence with fewer fields than its table decodes, its missing fields None."""
    nav_pvt = ubx.encode_frame(0x01, 0x07, bytes(92))
    bad_timels = ubx.encode_frame(0x01, 0x26, bytes(10))
    gsv = b"$GPGSV,1,1,00*79\r\n"
    short_zda = nmea.encode_frame(b"GPZDA,082710.00,16,09")
    framer = Framer()
    out = framer.push("port", nav_pvt + gsv + bad_timels + short_zda + ZDA + nav_pvt)
    assert out[0] == FrameFailure(FailureKind.DECODE, Protocol.UBX, 0x0126)
    assert isinstance(out[1], nmea.Zda) and out[1].year is None
    assert out[2] == nmea.decode(next(nmea.scan(ZDA)))
    assert len(out) == 3
    assert framer.skipped == {("ubx", 0x0107): 2, ("nmea", "GPGSV"): 1}


def test_transport_names_are_strings():
    with pytest.raises(TypeError):
        Framer().push(3, b"")


# ---- the delay a false start can add (DESIGN 2.2, framing helper) ------------------------

# A UBX start whose length field claims 65 535 bytes, as a stray B5 62 in other data can.
FALSE_UBX = bytes((0xB5, 0x62, 0x01, 0x07, 0xFF, 0xFF)) + bytes(10)


def test_false_start_holds_the_stream_without_bounds():
    framer = Framer()
    assert framer.push("port", FALSE_UBX + rtcm3.encode(MSG1005)) == []
    assert framer.pending("port") == len(FALSE_UBX) + len(rtcm3.encode(MSG1005))
    assert framer.flush("port") == [MSG1005]           # the end of the transport releases it


def test_max_frame_releases_frames_behind_a_false_start_at_once():
    framer = Framer(max_frame=2048)
    assert framer.push("port", FALSE_UBX + rtcm3.encode(MSG1005)) == [MSG1005]
    assert framer.pending("port") == 0
    # A frame that is long but within the bound is still waited for.
    whole = timels()
    assert framer.push("port", whole[:10]) == [] and framer.push("port", whole[10:]) != []


def test_max_wait_bounds_the_hold_in_receipt_time_and_replays_alike():
    t0 = GnssTime(2300, 1000.0)
    pieces = [("port", FALSE_UBX, t0), ("port", rtcm3.encode(MSG1005), t0 + 0.5),
              ("port", gad_bytes(), t0 + 1.5)]

    def run():
        framer, out = Framer(max_wait=1.0), []
        for transport, piece, received in pieces:
            out.append(framer.push(transport, piece, received))
        return out

    first = run()
    # Held while the false start is younger than 1 s, released once it is older.
    assert first == [[], [], [MSG1005, decoded_gad()]]
    assert run() == first                               # replay with the same times: the same
    # Without receipt times the time bound cannot act (and nothing reads a clock).
    framer = Framer(max_wait=1.0)
    assert [framer.push(t, p) for t, p, _ in pieces] == [[], [], []]


def test_bounds_are_checked():
    with pytest.raises(ValueError):
        Framer(max_frame=0)
    with pytest.raises(ValueError):
        Framer(max_wait=-1.0)
