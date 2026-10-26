import os
import subprocess
import sys
from pathlib import Path

import pytest

GALLERY = Path(__file__).resolve().parents[1] / "besca" / "examples" / "gallery_examples"
EXAMPLES = sorted(GALLERY.glob("*/plot_*.py"))


def test_gallery_examples_found():
    assert EXAMPLES, f"no gallery examples found in {GALLERY}"


@pytest.mark.parametrize("example", EXAMPLES, ids=lambda path: f"{path.parent.name}/{path.name}")
def test_gallery_example(example: Path, tmp_path: Path):
    # run each example as a script in its own process, with a non-interactive backend
    env = dict(os.environ, MPLBACKEND="Agg")
    result = subprocess.run(
        [sys.executable, str(example)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"{example.name} failed:\n{result.stdout[-5000:]}\n{result.stderr[-5000:]}"
