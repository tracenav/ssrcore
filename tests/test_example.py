"""examples/one_satellite.py runs and gives a complete record; the package's public names exist."""

from __future__ import annotations

import importlib.util
import runpy
from pathlib import Path

import ssrcore

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "one_satellite.py"


def load():
    spec = importlib.util.spec_from_file_location("one_satellite", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_example_gives_a_complete_corrected_record():
    example = load()
    engine, los = example.run()
    assert los.available and los.corrected and los.missing == {}
    assert los.signals == (ssrcore.Signal.GPS_L1C, ssrcore.Signal.GPS_L2W)
    assert engine.status().active.end_of_set
    lines = example.describe(los)
    assert len(lines) == 15 and not any("missing" in line for line in lines)


def test_the_example_runs_as_a_script(capsys):
    runpy.run_path(str(EXAMPLE), run_name="__main__")
    out = capsys.readouterr().out
    assert "corrected=True" in out and "troposphere" in out


def test_public_names():
    """DESIGN 2.1: ``from ssrcore import Engine, GnssTime, Sat, Signal`` and the rest of ``__all__``."""
    for name in ssrcore.__all__:
        assert hasattr(ssrcore, name), name
    from ssrcore import Engine, GnssTime, Sat, Signal  # noqa: F401  the form DESIGN 2.1 shows
