import os

# tests must not open plot windows (plt.show() blocks with interactive backends, e.g. on macOS)
os.environ.setdefault("MPLBACKEND", "Agg")

# The gallery examples are scripts that load large datasets when imported, so
# `pytest --doctest-modules` must not collect them; tests/test_gallery_examples.py runs them.
collect_ignore = ["besca/examples/gallery_examples"]
