"""Every ```python block of README.md runs, in order, in one namespace, from the
repository's root, as the README says; and what it prints is what it claims."""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"


def python_blocks(text: str) -> list[str]:
    """[design: the fenced blocks marked python, in document order]"""
    return re.findall(r"^```python\n(.*?)^```$", text, flags=re.MULTILINE | re.DOTALL)


def test_the_readme_has_its_examples():
    blocks = python_blocks(README.read_text(encoding="utf-8"))
    assert len(blocks) == 4, "feed, query, virtual base station, atmosphere providers"


def test_every_python_block_of_the_readme_runs(monkeypatch, capsys):
    monkeypatch.chdir(ROOT)
    monkeypatch.setattr(sys, "path", list(sys.path))      # the README puts examples/ on the path
    namespace: dict = {"__name__": "__readme__"}
    for k, block in enumerate(python_blocks(README.read_text(encoding="utf-8"))):
        exec(compile(block, f"README.md python block {k + 1}", "exec"), namespace)
    out = capsys.readouterr().out.splitlines()
    # feed: the ephemeris and the three SPARTN messages are stored
    assert out[:4] == ["Msg1019 stored", "Message stored", "Message stored", "Message stored"]
    # query: a complete record, two signals
    assert out[4].startswith("elevation 70.") and out[6].startswith("1C ") and out[7].startswith("2W ")
    # virtual base station: 1005, 1033, then the GPS MSM7
    assert out[8].startswith("[1005, 1033, 1077]")
    # providers: the caller's troposphere, no ionosphere outside the area
    assert out[9].startswith("fixed ")
    assert out[10] == "{'ionosphere': 'outside_area'}"
