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
## run the unit tests and doctests
uv run pytest
cd besca && uv run pytest --doctest-modules && cd ..
## run notebooks (should run without errors before a major release)
bash devtools/run_workbooks.bash
```

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
