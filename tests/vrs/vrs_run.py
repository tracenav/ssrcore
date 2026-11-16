"""The synthetic input of examples/vrs_stream.py, run once per scenario and decoded.

The tests of the observations read the run through the wire: every message is
encoded with gnsswire, scanned and decoded back, and each decoded cell is
turned into ranges by the restoration formulas of RTCM 10403.3 3.5.16.3.1.
"""

from __future__ import annotations

import functools
import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

from gnsswire import rtcm3

from ssrcore.model.types import GnssTime, Sat, Signal, System
from ssrcore.vrs import msm

EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "vrs_stream.py"


def example():
    """The example module, loaded once."""
    if "vrs_stream" not in sys.modules:
        spec = importlib.util.spec_from_file_location("vrs_stream", EXAMPLE)
        module = importlib.util.module_from_spec(spec)
        sys.modules["vrs_stream"] = module
        spec.loader.exec_module(module)
    return sys.modules["vrs_stream"]


@functools.lru_cache(maxsize=None)
def run(**events):
    """A run of the example's scenario with ``events`` changed, with the model of
    every observation as the engine gives it at that epoch: ``result.models`` maps
    (t, sat, signal) to (pseudorange_model, phase_model or None, wavelength,
    cycles of the source's observation)."""
    ex = example()
    models = {}

    def record(t, engine, station, messages):
        position = station.station.position
        for o in station.source.last:
            los = engine.query(o.sat, t, position)
            phase = los.phase_model(o.signal) if o.phase is not None else None
            models[(t, o.sat, o.signal)] = (los.pseudorange_model(o.signal), phase, o.wavelength, o.cycles)

    result = ex.run(ex.Scenario(**events), record)
    result.models = models
    return result


SYSTEM_OF = {prefix: system for system, prefix in msm.PREFIX.items()}
CODES = {system: {i: code for code, i in ids.items()} for system, ids in msm.SIGNAL_ID.items()}


@dataclass(frozen=True)
class Cell:
    """One decoded cell: the ranges restored from the wire, and the fields."""

    sat: Sat
    signal: Signal
    pseudorange: float
    phase: float | None
    lock: int
    rate: float | None
    channel: int | None
    cnr: int
    half_cycle: int


def cells(message: rtcm3.Msm) -> list[Cell]:
    """[spec: RTCM 10403.3 3.5.16.3.1, 3.5.16.3.4] the cells of a decoded MSM."""
    system = SYSTEM_OF[message.number // 10]
    level = message.number % 10
    sat_ids = [i for i in range(1, 65) if message.sat_mask >> (64 - i) & 1]
    sig_ids = [i for i in range(1, 33) if message.sig_mask >> (32 - i) & 1]
    width = len(sat_ids) * len(sig_ids)
    out, k = [], 0
    for a, sat_id in enumerate(sat_ids):
        for b, sig_id in enumerate(sig_ids):
            if not message.cell_mask >> (width - 1 - (a * len(sig_ids) + b)) & 1:
                continue
            rough = (message.rough_ms[a], message.rough_mod[a])
            p = msm.pseudorange_metres(*rough, message.fine_range[k], level)
            invalid = message.fine_phase[k] == msm.invalid_phaserange(level)
            phase = None if invalid else msm.phaserange_metres(*rough, message.fine_phase[k], level)
            rate = None
            if message.fine_rate is not None:
                rate = msm.rate_metres(message.rough_rate[a], message.fine_rate[k])
            channel = None
            if system == System.GLONASS and message.extended_info is not None:
                channel = message.extended_info[a] - 7
            out.append(Cell(Sat(system, sat_id), Signal(system, CODES[system][sig_id]), p, phase,
                            message.lock[k], rate, channel, message.cnr[k], message.half_cycle[k]))
            k += 1
    return out


def decoded_epochs(result) -> list[tuple[GnssTime, object, list]]:
    """(t, station, decoded messages) per epoch, every message through encode and decode."""
    out = []
    for t, _, station, messages in result.epochs:
        wire = b"".join(rtcm3.encode(m) for m in messages)
        out.append((t, station, [rtcm3.decode(f) for f in rtcm3.scan(wire)]))
    return out


def msm_cells(messages) -> dict[tuple[Sat, Signal], Cell]:
    """Every cell of an epoch's MSMs by satellite and signal."""
    return {(c.sat, c.signal): c for m in messages if isinstance(m, rtcm3.Msm) for c in cells(m)}
