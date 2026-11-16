"""Leap seconds, the engine's current time, and the resolution of SPARTN time tags.

DESIGN.md section 3.3, rules T1 to T5. The time base holds what the data have
said about time: GPS − UTC with its announced change (from NAV-TIMELS or a
``LeapSeconds``), and the current time, from receipt times when the caller
gives them and otherwise from unambiguous times in the data (T5). With these
it turns SPARTN time tags into GPS time (T1 to T4) and keeps aside the inputs
that cannot be placed in time yet.

The conversions between time scales, the calendar and the placement of a
16-bit tag in its half day are the toolkit's (``toolkit.timescale``, DESIGN
2.7); the time base decides only which leap second count applies at
which instant, which is what rule T4 needs.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Any

from ssrcore.model.types import (
    FeedReason,
    GnssTime,
    LeapSeconds,
    RoverPosition,
    SpartnTag,
    System,
    UtcTime,
)
from ssrcore.toolkit import timescale

_DAY = timescale.DAY_SECONDS


def utc_to_gps(year: int, month: int, day: int, seconds: float, leap: int) -> GnssTime:
    """The GPS time of a UTC date and second of the day under GPS − UTC ``leap``.

    [spec: TM-23 A.1.4 and IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a with A0 = A1 =
    0, through ``timescale.days_from_date`` and ``timescale.utc_to_gps``]
    """
    return GnssTime(*timescale.utc_to_gps(timescale.days_from_date(year, month, day), seconds, leap))


def gps_to_utc(t: GnssTime, leap: int) -> tuple[int, int, int, float]:
    """The UTC date and second of the day (year, month, day, seconds) at GPS time
    ``t`` under GPS − UTC ``leap``.

    [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a with A0 = A1 = 0 (UTC = GPS −
    leap); TM-23 A.1.4 through ``timescale.date_from_days``]
    """
    week, sow = timescale.add(t.week, t.seconds, -float(leap))
    day_of_week = int(sow // _DAY)
    year, month, day = timescale.date_from_days(week * 7 + day_of_week)
    return year, month, day, sow - day_of_week * _DAY


@dataclass(frozen=True, slots=True)
class _Held:
    item: Any
    system: System | None
    tag: SpartnTag | None
    received: GnssTime | None
    needs_time: bool = False


class TimeBase:
    """Leap seconds, the current time and the held inputs of one engine (DESIGN 3.2).

    [design: DESIGN 3.2, the time base store; one per engine, changed only
    by ``feed``]
    """

    def __init__(self, hold_limit: int = 256) -> None:
        """[design: DESIGN 3.3 rule T5, ``hold_limit`` from ``Config``]"""
        if hold_limit < 1:
            raise ValueError("hold_limit must be at least 1")
        self._hold_limit = hold_limit
        self._leap: LeapSeconds | None = None
        self._current: GnssTime | None = None
        self._held: deque[_Held] = deque()
        self._dropped = 0

    # ---- state ----------------------------------------------------------------------

    @property
    def leap_seconds(self) -> LeapSeconds | None:
        """[spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4] GPS − UTC and its
        announced change, when known."""
        return self._leap

    @property
    def current(self) -> GnssTime | None:
        """The current GPS time, None until the data or the caller have given one.

        [design: DESIGN 3.3 rule T5]
        """
        return self._current

    @property
    def dropped(self) -> int:
        """[design: DESIGN 3.3 rule T5] how many held inputs were dropped for want of room."""
        return self._dropped

    @property
    def held(self) -> int:
        """[design: DESIGN 3.3 rule T5] how many inputs are held."""
        return len(self._held)

    # ---- leap seconds (T4) ----------------------------------------------------------

    def set_leap_seconds(self, leap: LeapSeconds) -> None:
        """Take GPS − UTC and its announced change as the latest statement of them.

        [spec: u-blox F9 HPG 1.51 3.15.26 and u-blox X20 HPG 2.11 3.14.22
        (NAV-TIMELS); design: the latest statement replaces the earlier one,
        since a source reports its current knowledge, which includes every
        change it has seen]
        """
        if not isinstance(leap, LeapSeconds):
            raise TypeError(f"expected LeapSeconds, not {type(leap).__name__}")
        self._leap = leap

    def leap_at(self, t: GnssTime) -> int | None:
        """GPS − UTC at GPS time ``t``, or None without leap seconds.

        [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a and c; ``LeapSeconds.at``]
        """
        return None if self._leap is None else self._leap.at(t)

    def leap_at_utc(self, utc: GnssTime) -> int | None:
        """GPS − UTC for a UTC instant given on the GPS week grid without leap seconds.

        ``utc`` is the UTC reading as if it were GPS time (``utc_to_gps`` with
        a count of 0). The new count applies from the UTC instant at which it
        takes effect, the first second of the next UTC day: ``utc + next``
        reaches ``effective`` exactly then. A UTC second that does not exist
        (after a negative leap second) gets the old count.
        [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 a to c]
        """
        ls = self._leap
        if ls is None:
            return None
        if ls.next is not None and ls.effective is not None and utc + ls.next >= ls.effective:
            return ls.next
        return ls.current

    # ---- the current time (T5) ------------------------------------------------------

    def advance(self, t: GnssTime) -> bool:
        """Make ``t`` the current time unless the current time is later; True if it moved.

        [design: DESIGN 3.3 rule T5, the current time never moves backwards]
        """
        if self._current is None or t > self._current:
            self._current = t
            return True
        return False

    def observe_received(self, received: GnssTime) -> bool:
        """A receipt time the caller gave with an input. [design: DESIGN 3.3 rule T5]"""
        return self.advance(received)

    def observe_tag(self, tag: SpartnTag, system: System = System.GPS) -> GnssTime | None:
        """A 32-bit SPARTN tag seen in the data: its GPS time, which becomes the
        current time if later; None when it cannot be converted (no leap
        seconds for a GLONASS tag) or the tag is a 16-bit one.

        [spec: SPARTN Table 7.1 TF009, the 32-bit tag is unambiguous; design:
        DESIGN 3.3 rule T5]
        """
        if not tag.full:
            return None
        t = self.resolve(tag, system)
        if not isinstance(t, GnssTime):
            return None
        self.advance(t)
        return t

    def observe_utc(self, utc: UtcTime) -> GnssTime | None:
        """A UTC date and time (NMEA ZDA): its GPS time, which becomes the current
        time if later; None without leap seconds.

        [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4; design: DESIGN 3.3 rules T4
        and T5, no time is made from a guessed leap second count]
        """
        seconds = utc.hour * 3600.0 + utc.minute * 60.0 + utc.second
        t = self._from_utc(utc.year, utc.month, utc.day, seconds)
        if t is not None:
            self.advance(t)
        return t

    def observe_rover(self, rover: RoverPosition) -> GnssTime | None:
        """The UTC time of day of a GGA fix, placed in the UTC day that puts it
        nearest the current time; it becomes the current time if later. None
        without a current time or without leap seconds.

        [spec: u-blox F9 HPG 1.51 and u-blox X20 HPG 2.11 2.7.5.1, GGA time is UTC
        of the day; design: DESIGN 3.3 rule T5, a time of day is placed in the day the time base
        already knows; nearest is unambiguous while the time base is within
        12 h of the fix]
        """
        now = self._current
        leap = None if now is None else self.leap_at(now)
        if now is None or leap is None:
            return None
        best = None
        for shift in (-_DAY, 0.0, _DAY):
            year, month, day, _ = gps_to_utc(now + shift, leap)
            t = self._from_utc(year, month, day, rover.utc_seconds)
            if t is not None and (best is None or abs(t - now) < abs(best - now)):
                best = t
        if best is not None:
            self.advance(best)
        return best

    def _from_utc(self, year: int, month: int, day: int, seconds: float) -> GnssTime | None:
        # The leap second itself, 23:59:60, belongs to the day it ends and takes
        # that day's count [spec: IS-GPS-200N with IRN-IS-200N-003 and -004 20.3.3.5.2.4 b].
        grid = utc_to_gps(year, month, day, min(seconds, _DAY - 1.0), 0)
        leap = self.leap_at_utc(grid)
        return None if leap is None else utc_to_gps(year, month, day, seconds, leap)

    # ---- SPARTN time tags (T1 to T4) ------------------------------------------------

    def missing(self, tag: SpartnTag | None, system: System | None,
                received: GnssTime | None = None, needs_time: bool = False) -> FeedReason | None:
        """What the time base lacks to place an input in GPS time, or None.

        GLONASS inputs need leap seconds (T4); a 16-bit tag, and an input that
        ``needs_time`` (an RTCM ephemeris, whose week is resolved near its
        receipt, DESIGN 2.2), need the receipt time or the current time (T2, T5).
        [spec: GLONASS ICD 3.3.3; SPARTN Table 7.1 TF009; design: DESIGN 3.3
        rules T4 and T5]
        """
        if system == System.GLONASS and self._leap is None:
            return FeedReason.NO_LEAP_SECONDS
        if (needs_time or (tag is not None and not tag.full)) and received is None and self._current is None:
            return FeedReason.NO_TIME
        return None

    def resolve(self, tag: SpartnTag, system: System = System.GPS,
                received: GnssTime | None = None) -> GnssTime | FeedReason:
        """The GPS time of a SPARTN time tag, or the reason it cannot be had yet.

        ``system`` is the constellation whose time scale the tag is in: the
        message's for OCB and HPAC, GPS for GAD and messages of several
        constellations (T1). A 32-bit tag counts seconds from 2010-01-01
        00:00:00 of that scale (T3). A 16-bit tag counts seconds of the half
        day of that scale and is placed in the half day that puts it nearest
        ``received``, or the current time when there is no receipt time (T2,
        T5), by ``timescale.spartn_resolve_half_day`` (the earlier half day on
        a tie). GLONASS time is UTC + 3 h, so a GLONASS tag needs the leap
        seconds (T4). [spec: SPARTN Table 7.1 TF009; GLONASS ICD 3.3.3;
        design: DESIGN 3.3 rule T2; the reference is the whole second of the
        reference time, as a 32-bit tag counts whole seconds]
        """
        reason = self.missing(tag, system, received)
        if reason is not None:
            return reason
        value = tag.value
        if not tag.full:
            reference = received if received is not None else self._current
            assert reference is not None   # missing() said NO_TIME otherwise
            at = timescale.spartn_tag_from_gps(reference.week, reference.seconds, system,
                                               self.leap_at(reference) or 0)
            value = timescale.spartn_resolve_half_day(value, at)
        return self._spartn_to_gps(system, value)

    def spartn_to_gps(self, system: System, seconds: int) -> GnssTime | FeedReason:
        """GPS time of a 32-bit tag, ``seconds`` since 2010-01-01 00:00:00 of the system's scale (T1, T3, T4).

        [spec: SPARTN Table 7.1 TF009; GLONASS ICD 3.3.3, GLONASS time is
        UTC(SU) + 3 h and follows the UTC leap seconds; design: Q2, the 2010
        epoch is read in the scale of the tag, and for GLONASS the seconds
        are those of its clock reading, 86400 to the day, so that the 16-bit
        tag is the second of the GLONASS half day]
        """
        if system == System.GLONASS and self._leap is None:
            return FeedReason.NO_LEAP_SECONDS
        return self._spartn_to_gps(system, seconds)

    def _spartn_to_gps(self, system: System, tag: int) -> GnssTime:
        leap = 0
        if system == System.GLONASS:
            # GLONASS time is a UTC reading (plus 3 h): the count is chosen as for UTC.
            leap = self.leap_at_utc(GnssTime(*timescale.spartn_tag_to_gps(tag, system, 0))) or 0
        return GnssTime(*timescale.spartn_tag_to_gps(tag, system, leap))

    # ---- inputs held until they can be placed in time (T4, T5) -----------------------

    def hold(self, item: Any, system: System | None = None, tag: SpartnTag | None = None,
             received: GnssTime | None = None, needs_time: bool = False) -> Any | None:
        """Keep ``item`` until ``missing`` has nothing to say for it; return the
        oldest held item if room had to be made for this one, else None.

        [design: DESIGN 3.3 rules T4 and T5, at most ``hold_limit`` inputs
        held, the oldest dropped and counted]
        """
        dropped = None
        if len(self._held) >= self._hold_limit:
            dropped = self._held.popleft().item
            self._dropped += 1
        self._held.append(_Held(item, system, tag, received, needs_time))
        return dropped

    def release(self) -> list[Any]:
        """The held items that can now be placed in time, in the order they arrived;
        the others stay held.

        [design: DESIGN 3.3 rule T5, released in arrival order when the time
        becomes known]
        """
        ready, waiting = [], deque()
        for h in self._held:
            (waiting if self.missing(h.tag, h.system, h.received, h.needs_time) else ready).append(h)
        self._held = waiting
        return [h.item for h in ready]
