"""ssrcore: GNSS state-space corrections engine.

Layout and rules: docs/DESIGN.md. Public interface (DESIGN 2.1): the names
below, from :mod:`ssrcore.engine`, :mod:`ssrcore.record` and the input types
of :mod:`ssrcore.model.types`; the framing helper
:class:`~ssrcore.adapters.framing.Framer`; the adapters of
:mod:`ssrcore.adapters`; the stateless :mod:`ssrcore.toolkit`; and
:mod:`ssrcore.vrs`. The atmosphere sources of ``Config.atmosphere`` and the
provider interface (rule A9) are here too; a caller's provider returns the
value types of :mod:`ssrcore.model.providers` (``Zenith``, ``SlantTec``, and
``Absent`` of :mod:`ssrcore.model.atmosphere`).

    from ssrcore import Engine, Framer, GnssTime, Sat, Signal
"""

from ssrcore.adapters.framing import Framer
from ssrcore.engine import Engine, EngineStatus
from ssrcore.model.providers import AtmosphereProvider, Kind, Origin
from ssrcore.model.types import (
    CallerSource,
    Config,
    FeedOutcome,
    FeedReason,
    FeedResult,
    FrameFailure,
    GnssTime,
    LeapSeconds,
    MaxAge,
    Missing,
    ModelSource,
    RoverPosition,
    Sat,
    Signal,
    SpartnSource,
    StreamId,
    System,
    Term,
    UtcTime,
)
from ssrcore.record import LineOfSight

__version__ = "0.1.0"

__all__ = [
    "AtmosphereProvider", "CallerSource", "Config", "Engine", "EngineStatus", "FeedOutcome", "FeedReason",
    "FeedResult", "FrameFailure", "Framer", "GnssTime", "Kind", "LeapSeconds", "LineOfSight", "MaxAge", "Missing",
    "ModelSource", "Origin", "RoverPosition", "Sat", "Signal", "SpartnSource", "StreamId", "System", "Term",
    "UtcTime", "__version__",
]
