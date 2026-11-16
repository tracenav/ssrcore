"""Pieces of byte streams of named transports to typed messages, through gnsswire.

DESIGN.md section 2.2, "Framing helper". At the edge: bytes in, messages out,
no I/O and no clock.

A recording or a live link delivers pieces of byte streams: a frame may start
in one piece and end in a later one, a piece may hold several frames and
bytes that belong to none, and several transports (an IP and an L-band
correction stream, a receiver port) are interleaved. ``Framer`` keeps one
buffer per named transport, asks gnsswire's ``parse`` of each codec whether a
frame starts at each position, and returns the decoded messages, and a
``FrameFailure`` for each frame that could not be decoded, in the order their
last byte arrived.

The four codecs are told apart by the first byte of a frame: SPARTN 0x73
(SPARTN Table 7.1 TF001), RTCM 3 0xD3 (RTCM 4, Table 4-1), UBX 0xB5 0x62
(u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 3.2) and NMEA ``$`` (both 2.1).
Every transport may carry any of them, as a receiver port does.

**Delay.** A frame is released when its last byte and its check have
arrived; the framer holds bytes only while a frame that starts at them is
still incomplete. A false start (a start byte inside other data, with a
length field that claims more than follows) holds the bytes after it until
enough have arrived to fail its check: up to 65 543 bytes for UBX. For live
use the caller bounds this with ``max_frame`` (the bytes a frame may claim)
and ``max_wait`` (the seconds of receipt time a start may wait); a start
beyond either is taken as false. Both act only on the bytes and on the
receipt times the caller passes, so a replay of the same pieces with the
same times gives the same messages (DESIGN 2.2, framing helper).
"""

from __future__ import annotations

import re
from collections import deque
from collections.abc import Iterable, Mapping
from typing import Any

from gnsswire import nmea, rtcm3, spartn, ubx

from ssrcore.adapters.spartn_ocb import header_input
from ssrcore.adapters.ubx import pmp_user_data
from ssrcore.model.types import Carried, FailureKind, FrameFailure, GnssTime, Protocol

# [spec: SPARTN Table 7.1 TF001; RTCM 4 Table 4-1; u-blox F9 HPG 1.51 and u-blox X20
#  HPG 2.11 3.2 and 2.1]
# the first byte of a frame of each codec.
_CODECS: dict[int, Any] = {spartn.PREAMBLE: spartn, rtcm3.PREAMBLE: rtcm3,
                           ubx.SYNC[0]: ubx, ord("$"): nmea}
_ANY_START = re.compile(b"[" + b"".join(re.escape(bytes((b,))) for b in _CODECS) + b"]")
# [spec: as _CODECS] the first byte of each protocol's frames.
_FIRST_BYTE = {Protocol.SPARTN: spartn.PREAMBLE, Protocol.RTCM3: rtcm3.PREAMBLE,
               Protocol.UBX: ubx.SYNC[0], Protocol.NMEA: ord("$")}
# [design: DESIGN 2.2, the L-band service delivers SPARTN in RXM-PMP user data;
#  only SPARTN is looked for there, so that other codecs' start bytes in it are
#  not taken for frames]
_SPARTN_ONLY = re.compile(re.escape(bytes((spartn.PREAMBLE,))))

# The inner stream of a transport: the SPARTN bytes of its RXM-PMP messages.
_PMP = "RXM-PMP"


class Framer:
    """Frames and decodes pieces of byte streams of named transports.

    ``push(transport, piece)`` returns what the piece completed; ``flush``
    ends a transport. Nothing here reads a clock: a receipt time given with a
    piece is copied into the failures it produces, and the caller passes it
    with the messages to ``Engine.feed``.
    [design: DESIGN 2.2, a shared framer makes every caller frame the same way]
    """

    def __init__(self, max_frame: int | None = None, max_wait: float | None = None,
                 protocols: Mapping[str, Iterable[Protocol]] | None = None) -> None:
        """[design: DESIGN 2.2, one buffer per named transport; ``max_frame`` and
        ``max_wait`` bound the delay a false start can add (None: no bound beyond
        each protocol's own maximum length); ``protocols`` names, for a
        transport the caller knows to carry only some protocols, those
        protocols: a start byte of another one there is never a frame, so it
        cannot hold the stream (DESIGN 2.2, *Delay*); a transport it does not
        name carries every protocol]"""
        if max_frame is not None and max_frame < 1:
            raise ValueError("max_frame must be at least 1 byte")
        if max_wait is not None and not max_wait >= 0.0:
            raise ValueError("max_wait must be zero or more seconds")
        self.max_frame = max_frame
        self.max_wait = max_wait
        self._starts: dict[str, re.Pattern] = {}
        for transport, carried in (protocols or {}).items():
            first = sorted({_FIRST_BYTE[Protocol(p)] for p in carried})
            if not first:
                raise ValueError(f"transport {transport!r} carries no protocol")
            self._starts[transport] = re.compile(b"[" + b"".join(re.escape(bytes((b,))) for b in first) + b"]")
        self._buffers: dict[tuple[str, str], bytearray] = {}
        # receipt times of the kept bytes: [end offset in the buffer, received] per piece
        self._arrivals: dict[tuple[str, str], deque[list]] = {}
        self._skipped: dict[tuple[str, int | str], int] = {}

    @property
    def skipped(self) -> dict[tuple[str, int | str], int]:
        """UBX and NMEA messages without a gnsswire table, counted per (protocol,
        message): a UBX class × 256 + ID, an NMEA address such as ``GPGSV``.

        [design: DESIGN 3.5, a receiver port sends many messages that are not
        inputs; they are counted here and not passed to the engine]
        """
        return dict(self._skipped)

    @property
    def transports(self) -> tuple[str, ...]:
        """[design: DESIGN 2.2] the names of the transports seen, in order of first piece."""
        return tuple(dict.fromkeys(name for name, _ in self._buffers))

    def pending(self, transport: str) -> int:
        """[design: DESIGN 2.2] bytes of ``transport`` kept because a frame may still complete."""
        return sum(len(b) for (name, _), b in self._buffers.items() if name == transport)

    def push(self, transport: str, piece: bytes, received: GnssTime | None = None) -> list[Any]:
        """The messages and failures completed by ``piece`` of ``transport``.

        Bytes at which no frame starts are skipped one at a time; bytes at
        which a frame may start but is not complete are kept for the next
        piece. A CRC or checksum failure is not a frame: its bytes are
        skipped like any others.
        [spec: gnsswire ``parse`` (Frame, NeedMore, Invalid); SPARTN Table 7.1
        TF006, TF018; RTCM 4; u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 3.4 and
        2.1; design: DESIGN 3.5, a failed CRC reaches the engine as nothing]
        """
        if not isinstance(transport, str):
            raise TypeError("a transport is named by a string")
        out: list[Any] = []
        self._extend((transport, ""), bytes(piece), received, out, final=False)
        return out

    def flush(self, transport: str, received: GnssTime | None = None) -> list[Any]:
        """The end of ``transport``: what its kept bytes still hold, taken as complete.

        A frame cut off by the end is skipped like any other bytes, and the
        frames after a false start that was waiting for more are found.
        [spec: gnsswire ``scan``, "the buffer is taken as complete"; DESIGN
        3.5, bytes that stop inside a frame give nothing]
        """
        out: list[Any] = []
        for key in ((transport, ""), (transport, _PMP)):   # the outer stream may add to the inner
            if key in self._buffers:
                self._extend(key, b"", received, out, final=True)
        return out

    def push_all(self, pieces: Iterable[tuple[str, bytes]]) -> list[Any]:
        """Every message of a sequence of (transport, piece), then of each transport's end.

        [design: DESIGN 2.2, a recording often holds pieces without receipt
        times, and the time base then takes its current time from the data
        (rule T5)]
        """
        out: list[Any] = []
        for transport, piece in pieces:
            out += self.push(transport, piece)
        for transport in self.transports:
            out += self.flush(transport)
        return out

    # ---- framing ---------------------------------------------------------------------

    def _extend(self, key: tuple[str, str], data: bytes, received: GnssTime | None,
                out: list[Any], final: bool) -> None:
        buf = self._buffers.setdefault(key, bytearray())
        buf += data
        arrivals = self._arrivals.setdefault(key, deque())
        if data:
            arrivals.append([len(buf), received])
        starts = _SPARTN_ONLY if key[1] == _PMP else self._starts.get(key[0], _ANY_START)
        pos, n = 0, len(buf)
        while pos < n:
            m = starts.search(buf, pos)
            if m is None:
                pos = n
                break
            pos = m.start()
            codec = _CODECS[buf[pos]]
            result = codec.parse(buf, pos)
            if isinstance(result, codec.NeedMore):
                if not final and not self._too_long(result.needed, key, pos, received):
                    break
                pos += 1
            elif isinstance(result, codec.Invalid):
                pos += 1
            else:
                pos += result.length
                self._decode(key[0], codec, result, received, out)
        del buf[:pos]
        while arrivals and arrivals[0][0] <= pos:
            arrivals.popleft()
        for a in arrivals:
            a[0] -= pos
        if final and not buf:
            del self._buffers[key]
            del self._arrivals[key]

    def _too_long(self, needed: int, key: tuple[str, str], pos: int, received: GnssTime | None) -> bool:
        """Whether an incomplete start at ``pos`` is past the caller's bounds.

        [design: DESIGN 2.2, framing helper: a start that claims more than
        ``max_frame`` bytes, or whose first byte arrived more than ``max_wait``
        seconds of receipt time before ``received``, is taken as false]
        """
        if self.max_frame is not None and needed > self.max_frame:
            return True
        if self.max_wait is None or received is None:
            return False
        for end, when in self._arrivals[key]:
            if pos < end:
                return when is not None and received - when > self.max_wait
        return False

    def _decode(self, transport: str, codec: Any, frame: Any, received: GnssTime | None,
                out: list[Any]) -> None:
        if codec is spartn:
            out.extend(_spartn(frame, received))
        elif codec is rtcm3:
            out.extend(_rtcm3(frame, received))
        elif codec is ubx:
            message = _ubx(frame, received)
            if isinstance(message, (ubx.RxmPmp, ubx.RxmPmpV0)):
                self._extend((transport, _PMP), pmp_user_data(message), received, out, final=False)
            elif isinstance(message, ubx.RxmSfrbx):
                # [design: DESIGN 2.2, navigation words are joined per transport
                #  (adapters.navdata), so they leave the framer with their transport]
                out.append(Carried(transport, message))
            elif message is not None:
                out.append(message)
            else:
                self._skip(Protocol.UBX, (frame.msg_class << 8) | frame.msg_id)
        else:
            message = _nmea(frame, received)
            if message is not None:
                out.append(message)
            else:
                self._skip(Protocol.NMEA, bytes(frame.payload).split(b",", 1)[0].decode("ascii", "replace"))

    def _skip(self, protocol: Protocol, message: int | str) -> None:
        key = (protocol.value, message)
        self._skipped[key] = self._skipped.get(key, 0) + 1


def _spartn(frame: spartn.Frame, received: GnssTime | None) -> list[Any]:
    """A SPARTN message, the header of an OCB message whose subtype has no table, or a failure.

    [spec: SPARTN 8.15 (encrypted payloads); 8.4 item 4, data after the table
    is kept by gnsswire's tolerant ``decode`` (``extension``), not refused;
    8.4 item 5, decoders read SF010 of every OCB message "even though they
    may not understand the data", which ``decode_ocb_header`` does; design:
    DESIGN 3.5, every other frame that cannot be read is a failure of its
    protocol, caught as the codec's one error base]
    """
    where = dict(protocol=Protocol.SPARTN, message_type=frame.message_type,
                 subtype=frame.description.subtype, received=received)
    try:
        return [spartn.decode(frame)]
    except spartn.Encrypted:
        return [FrameFailure(FailureKind.ENCRYPTED, **where)]
    except spartn.Unsupported:
        if frame.message_type == spartn.Ocb.MESSAGE_TYPE:
            try:
                return [header_input(frame, spartn.decode_ocb_header(frame))]
            except spartn.Error:
                pass
        return [FrameFailure(FailureKind.UNSUPPORTED, **where)]
    except spartn.Error:
        return [FrameFailure(FailureKind.DECODE, **where)]


def _rtcm3(frame: rtcm3.Frame, received: GnssTime | None) -> list[Any]:
    """[spec: RTCM Table 3.4-1, DF002 names the message; 3.5.17, data appended after a
    message's table is kept by gnsswire's tolerant ``decode`` (``extension``); IGS
    SSR v1.00 Table 7 IDF002, the sub-type of message 4076; design: DESIGN 3.5,
    caught as the codec's one error base, a 4076 failure counted by its sub-type]"""
    number = rtcm3.message_number(frame)
    subtype = rtcm3.igs_subtype(bytes(frame.payload)) if number == rtcm3.IGS_SSR else None
    try:
        return [rtcm3.decode(frame)]
    except rtcm3.Unsupported:
        return [FrameFailure(FailureKind.UNSUPPORTED, Protocol.RTCM3, number, subtype, received)]
    except rtcm3.Error:
        return [FrameFailure(FailureKind.DECODE, Protocol.RTCM3, number, subtype, received)]


def _ubx(frame: ubx.Frame, received: GnssTime | None) -> Any:
    """A UBX message; None for a message without a gnsswire table, a failure for one
    that does not fit its table.

    [design: DESIGN 3.5, a receiver port sends many UBX messages that are not
    inputs: they are counted by the framer (``Framer.skipped``), not fed; a
    message with a table that it does not fit is a failure, its type class ×
    256 + ID]
    """
    try:
        return ubx.decode(frame)
    except ubx.Unsupported:
        return None
    except ubx.Error:
        return FrameFailure(FailureKind.DECODE, Protocol.UBX, (frame.msg_class << 8) | frame.msg_id,
                            None, received)


def _nmea(frame: nmea.Frame, received: GnssTime | None) -> Any:
    """An NMEA sentence; None for a sentence without a gnsswire table, a failure for
    one that does not fit its table.

    [design: as for UBX, a receiver sends many sentences the engine does not use;
    ``FrameFailure`` carries numbers only, so a sentence's address is not in it;
    gnsswire's tolerant ``decode`` gives a sentence with fewer or more fields
    than its table, the missing ones ``None``, which the NMEA adapter reads as
    empty]
    """
    try:
        return nmea.decode(frame)
    except nmea.Unsupported:
        return None
    except nmea.Error:
        return FrameFailure(FailureKind.DECODE, Protocol.NMEA, None, None, received)
