Scripts and tools used by BESCA developers
===
Klas Hatj, Roland Schmucki, Jitao David Zhang

## Set up the development environment

Dependencies are managed with [uv](https://docs.astral.sh/uv/). The exact versions of all
dependencies are pinned in `uv.lock` at the root of the repository. To create an environment
(in `.venv`) with besca installed in editable mode and all locked dependencies, go to the
root directory and run:

```bash
uv sync
```

`bash devtools/install_besca_editable.bash` additionally registers the environment as the
`besca_dev` Jupyter kernel.

## Test besca package

```bash
## run the unit tests and the gallery examples
uv run pytest
## run the doctests
cd besca && uv run pytest --doctest-modules && cd ..
## run notebooks (should run without errors before a major release)
bash devtools/run_workbooks.bash
```

`uv run pytest` runs every script in `besca/examples/gallery_examples` (with a non-interactive
matplotlib backend). The tests and examples download about 3.5 GB of datasets from Zenodo on the
first run and cache them in `~/.cache/besca/datasets` (see `BESCA_DATA_DIR` in the README);
`uv run pytest -k "not gallery"` skips the examples. The doctest run does not collect the gallery
examples (see `conftest.py`); the documentation build executes them with sphinx-gallery.

## Update dependencies

To update the locked dependencies, e.g. to pick up security fixes, and check them for known
vulnerabilities:

```bash
uv lock --upgrade
uv audit --preview-features audit-command
```

`uv audit` (still an experimental uv command, hence `--preview-features`) checks every locked
package, including the `dev` and `docs` groups, against the OSV vulnerability database and exits
with an error if it finds a known vulnerability.

If you need a `requirements.txt` to recreate the environment with pip, export it from the lock file:

```bash
uv export --no-emit-project --format requirements-txt -o requirements.txt
```

## Release a new version

The package version is derived from the latest git tag by [versioneer](https://github.com/python-versioneer/python-versioneer)
(tags have no `v` prefix, e.g. `2.5.5`). Untagged commits get a local version such as
`2.5.5+12.gabc1234`. To release, tag the release commit on `main` and push the tag; pushing a tag
publishes the package to PyPI (see `.github/workflows/pypi-auto-deploy.yml`):

```bash
git switch main && git pull
git tag -a 3.0.0 -m "besca 3.0.0"
git push origin 3.0.0
```

Afterwards, build the documentation with the manually triggered "Build Sphinx Documentation" workflow
on `main`; it publishes the docs to the `gh-pages` branch.
