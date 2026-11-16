"""The tables of kernels/src/astro.c are their sources' terms, selected as stated.

tools/astro_tables.py regenerates them from ELP 2000-82B (ELP1 to ELP36) and
VSOP87D.ear when those files are staged in refs/Q11_sun_moon_tides/
(THIRD_PARTY_NOTICES.md, refs/README.md); without them the test skips with
the reason. The sources are never committed (DESIGN 5.1, Specifications are
not committed).
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "refs" / "Q11_sun_moon_tides"
NEEDED = [f"ELP{k}" for k in range(1, 37)] + ["VSOP87D.ear"]


def _tool():
    spec = importlib.util.spec_from_file_location("astro_tables", ROOT / "tools" / "astro_tables.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_astro_tables_regenerate_from_their_sources():
    """DESIGN 2.7 astro: every ELP term with |A| >= 0.2 and VSOP87D term with |A| >= 1e-6, verbatim."""
    absent = [name for name in NEEDED if not (SOURCES / name).is_file()]
    if absent:
        pytest.skip(f"refs/Q11_sun_moon_tides/ absent or incomplete (missing {absent[0]}); "
                    "the ELP 2000-82B and VSOP87D files are not distributed")
    tool = _tool()
    assert tool.generate(SOURCES) == tool.committed()


def test_astro_tables_markers_and_counts():
    """The committed tables sit between the generator's markers and declare their term counts."""
    tool = _tool()
    tables = tool.committed()
    for name in ("ELP_MAIN_LONGITUDE", "ELP_MAIN_LATITUDE", "ELP_MAIN_DISTANCE",
                 "ELP_PERTURBATIONS", "VSOP87D_EARTH"):
        assert f" {name}[" in tables
    assert "/* VSOP87D Earth: 69 of 2425 terms. */" in tables
