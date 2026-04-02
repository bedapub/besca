# Besca v3.0.0 — Specification

## Problem Statement

Besca (BEyond Single Cell Analysis) is a single-cell RNA-seq analysis toolkit currently at v2.5, pinned to Python 3.9 with outdated dependencies, a legacy build system (setup.py + versioneer), 30 unresolved GitHub issues, stale documentation, and broken CI workflows. The project needs a major version upgrade to modernize the stack, fix all outstanding issues, and validate that the toolkit still produces correct results against the published benchmarks.

## Scope

This spec covers the full Besca v3.0.0 release on a new `besca3` branch:

1. Python 3.12 upgrade with all dependencies updated
2. Build system migration from setup.py/versioneer to pyproject.toml
3. Resolution of all 30 open GitHub issues
4. Removal of deprecated R dependencies
5. Rewrite of the C `reformat` binary in Python
6. CI/CD workflow modernization
7. Sphinx documentation rebuild and automated gh-pages deployment
8. Benchmark replication across all 10 Zenodo datasets (Supplementary Table 7)

---

## 1. Branch and Version Setup

### Requirements
- Create branch `besca3` from `main`
- Set version to `3.0.0`
- Remove versioneer; version managed via `pyproject.toml` metadata

### Acceptance Criteria
- `besca3` branch exists and diverges from `main` at HEAD
- `import besca; besca.__version__` returns `"3.0.0"`
- No versioneer files remain (`versioneer.py`, `besca/_version.py` generated content, `setup.cfg` versioneer section)

---

## 2. Build System Migration

### Requirements
- Replace `setup.py` + `setup.cfg` + `versioneer.py` with `pyproject.toml`
- Use `setuptools` as build backend (or `hatchling`/`flit` if simpler)
- Use `setuptools-scm` or static version in `pyproject.toml`
- Move all metadata (name, version, description, classifiers, URLs, license) to `pyproject.toml`
- Move `requirements.txt` contents into `[project.dependencies]`
- Update classifiers to `Programming Language :: Python :: 3.12`
- Keep `package_data` declarations for `.h5ad`, `.css`, `.tsv`, `.gmt` files
- Remove `MANIFEST.in` if no longer needed

### Acceptance Criteria
- `pip install .` succeeds from the repo root
- `pip install -e .` succeeds for development
- `python -m build` produces a valid sdist and wheel
- No `setup.py`, `setup.cfg`, `versioneer.py`, or `MANIFEST.in` remain

---

## 3. Python 3.12 and Dependency Upgrade

### Current State
- Python 3.9, numpy<1.24, scanpy>=1.7.2, anndata>=0.7.4, scvi-tools (old), jax 0.4.7, pytorch-lightning 1.6, torch (not in main deps but in devtools)

### Requirements
- Target Python >=3.12
- Update all dependencies to latest versions compatible with Python 3.12:
  - `scanpy` (latest, likely >=1.10)
  - `anndata` (latest, likely >=0.10)
  - `scvi-tools` (latest compatible, likely >=1.0 — note: API changed significantly)
  - `scipy`, `numpy`, `pandas`, `matplotlib`, `seaborn` — latest stable
  - `scikit-learn` — latest stable
  - `bbknn` — use official PyPI release, not git fork
  - `leidenalg`, `python-igraph` — latest
  - `scanorama`, `scvelo`, `umap-learn` — latest compatible
  - `gseapy` — latest
  - `plotly`, `flask-restful`, `dominate`, `requests` — latest
  - `session-info` — latest (or replace with `sinfo` if renamed)
  - `pydot` — latest
  - `nbclean` — check if still maintained, replace if not
  - `mygene` — latest
  - `deprecation` — latest
- Remove dependencies:
  - `rpy2`, `anndata2ri` — R integration removed
  - `ipython` (duplicate in requirements.txt) — keep one entry
- Update `environment.yml` to Python 3.12 with updated deps
- Update `environment.lock.yml` or remove if no longer used
- Update `devtools/requirements.txt` to match

### Key Migration Risks
- **scvi-tools**: v1.0+ has major API changes (module renames, training API). All `import scvi` usage in `_auto_annot.py` needs updating.
- **numpy**: v2.0+ drops some deprecated APIs. Check for `np.float`, `np.int`, `np.bool` usage.
- **anndata**: v0.10+ may change `.X` handling and views.
- **scanpy**: Check for deprecated function signatures.

### Acceptance Criteria
- `pip install .` succeeds on Python 3.12
- `import besca` succeeds without import errors
- All submodules (`besca.pl`, `besca.tl`, `besca.pp`, `besca.st`, `besca.export`, `besca.Import`, `besca.datasets`) import cleanly
- `pytest` passes (existing tests + doctests)

---

## 4. Remove R Dependencies

### Requirements
- Remove `Rlibs.R` from repo root and `besca/Rlibs.R`
- Remove or stub R-dependent functions:
  - `besca.pp.valOutlier` (uses `scater::isOutlier`)
  - `besca.pp.scTransform` (uses `Seurat::SCTransform`)
  - `besca.st.maxLikGlobalDimEst` (uses `intrinsicDimension`)
  - `besca.st.deviance` (uses `scry`)
  - `besca.st.dsb_normalize` (uses `dsb`)
- Remove R installation instructions from README.md
- Remove conda R dependency instructions
- If any of these functions have pure-Python alternatives, implement them

### Acceptance Criteria
- No R imports anywhere in the codebase
- No `rpy2` or `anndata2ri` in dependencies
- Functions that relied on R either removed or replaced with Python equivalents
- README has no R installation section

---

## 5. Rewrite `reformat` Binary in Python

### Requirements
- Understand what the C `reformat` binary does (located in `besca/export/reformat`)
- Rewrite equivalent functionality in Python
- Remove the compiled binary from the repo
- Update `besca/export/` to use the Python implementation
- Remove the "set executable flag" instructions from README

### Acceptance Criteria
- No compiled C binary in the repository
- Python `reformat` function produces identical output to the C version
- Export functionality works end-to-end
- README updated to remove binary instructions

---

## 6. Resolve All 30 Open Issues

Issues grouped by category:

### Bugs (fix in code)
| # | Title | Approach |
|---|-------|----------|
| 386 | Missing reformat Binary File | Covered by §5 (Python rewrite) |
| 384 | recluster() calculates HVG from previous ones | Fix HVG recalculation logic in `besca/tl/rc/` |
| 381 | Dependent versions of anndata and scipy too old | Covered by §3 (dependency upgrade) |
| 369 | seaborn jointplot fails with multiple positional args | Fix argument passing in plotting code |
| 367 | bc.st.clr_normalize() fails | Debug and fix CLR normalization |
| 363 | sig.read_annotconfig cell_type issue | Fix annotation config reader |
| 347 | besca installation failure | Covered by §2/§3 (build system + deps) |
| 346 | riverplot_2categories() figsize param broken | Fix figsize parameter handling |
| 345 | percent_mito reflects fraction not percentage | Fix to return actual percentage or document behavior |
| 344 | tl.annotate_cells_cluster function issue | Debug and fix annotation logic |
| 342 | new 10x mtx version | Support new 10x matrix format in Import |
| 341 | high memory consumption and NoneType error | Profile and fix memory issue |
| 323 | LFS quota exceeded | Remove LFS objects, use Zenodo for large files |
| 271 | KeyError 'base' in get_de() / additional_labeling() | Fix key access pattern |
| 265 | Remove duplicated log messages | Fix logging configuration |
| 81 | Display bug in gene_expr_split_stacked | Fix plotting function |
| 312 | use_raw=False not working in combined_signature_score() | Fix raw data handling |

### Enhancements (implement)
| # | Title | Approach |
|---|-------|----------|
| 320 | Allow flavor selection in highly_variable_genes() | Add `flavor` parameter |
| 180 | Additional params min_fract_pos/min_cells_per_group in perform_dge() | Add parameters |
| 188 | Upgrading plot functions with axes | Return/accept matplotlib axes objects |
| 181 | figure size for celllabel_quant_stackedbar | Add figsize parameter |
| 171 | Pending improvements in Silhouette Computation | Implement improvements per issue discussion |
| 156 | Add missing gene Ly6g and revise MGItoHGNC.csv | Update gene mapping file |
| 160 | Expanded cell filtering for datasets with protein data | Extend filtering logic |
| 34 | Automatic Detection of cutoff levels for signatures | Implement auto-cutoff |
| 18 | Multimodal Data | Add basic multimodal support or document scope |

### Documentation/Meta (update docs)
| # | Title | Approach |
|---|-------|----------|
| 317 | Move to poetry as dependency manager | Addressed by pyproject.toml migration (§2) |
| 298 | env yaml file | Update environment.yml (§3) |
| 272 | Storing public datasets in install dir problematic | Use XDG cache or configurable path |
| 177 | Getting started guide and interactive notebooks | Create/update getting started docs |
| 161 | Documentation for citeseq workflow missing | Add CITE-seq documentation |
| 135 | Documentation update | General docs refresh (§7) |

### Acceptance Criteria
- All 30 issues addressed with code changes, documentation, or justified closure
- Each issue gets a commit message referencing the issue number (`Fixes #N`)
- Issues closed on GitHub after fix is merged

---

## 7. Documentation and gh-pages

### Requirements
- Update `docs/source/conf.py`:
  - Version to `3.0.0`
  - Fix deprecated `app.add_stylesheet` → `app.add_css_file`
  - Update `github_version` from `master` to `besca3`
  - Update intersphinx mapping
  - Check `nbclean` usage (may need replacement)
- Update all RST/notebook documentation to reflect v3 changes
- Remove R-related documentation sections
- Add migration guide from v2.x to v3.0
- Build docs locally with `make html` to verify

### gh-pages Automation
- Create `.github/workflows/deploy-docs.yml`:
  - Trigger on push to `besca3` branch
  - Use Python 3.12
  - Install besca + sphinx dependencies
  - Build docs with `make html`
  - Deploy to `gh-pages` branch using `peaceiris/actions-gh-pages@v4` or similar
- Remove or update old `build_docs.yml` workflow

### Acceptance Criteria
- `make html` in `docs/` completes without errors or warnings
- gh-pages branch updated with v3.0 documentation
- GitHub Actions workflow auto-deploys on push to `besca3`
- https://bedapub.github.io/besca/ shows v3.0 docs after deployment

---

## 8. CI/CD Workflow Modernization

### Requirements
- **doc-tests.yml**: Update to Python 3.12, use pip instead of conda/mamba, update actions to v4
- **pypi-auto-deploy.yml**: Update to Python 3.12, use `pyproject.toml` build, update actions to v4
- **deploy-docs.yml**: New workflow (see §7)
- Add a basic CI workflow that runs on PR/push:
  - Install besca
  - Run `pytest`
  - Run `pytest --doctest-modules` on `besca/`

### Acceptance Criteria
- All workflows use `actions/checkout@v4`, `actions/setup-python@v5`
- No conda/mamba in CI (use pip)
- Tests run automatically on push to `besca3`
- PyPI deployment triggers on tags

---

## 9. Benchmark Replication (Supplementary Table 7)

### Requirements
- Download all 10 Zenodo datasets:
  - PBMC3k (doi:10.5281/zenodo.3948150)
  - Granja2019 (doi:10.5281/zenodo.3944753)
  - Kotliarov2020 (doi:10.5281/zenodo.3938290)
  - Smillie2019 (doi:10.5281/zenodo.3960617)
  - Martin2019 (doi:10.5281/zenodo.3862132)
  - Haber2017 (doi:10.5281/zenodo.3935782)
  - Lee2020 (doi:10.5281/zenodo.3967538)
  - Segerstolpe2016 (doi:10.5281/zenodo.3928276)
  - Peng2019 (doi:10.5281/zenodo.3969339)
  - Baron2016 (doi:10.5281/zenodo.3968315)

- Replicate analyses from `besca_publication_results` repository:
  - **Sig-annot comparison**: Signature-based annotation at levels 1-3 across Lee2020, Peng2019 datasets
  - **Auto-annot comparison**: Automated annotation using logistic regression (log reg, log reg th), SVM, scANVI across all datasets with various training data combinations
  - **CellAssign**: PBMC3k with Besca/BescaF/BescaFv2/BescaFvvs signatures
  - **SingleR**: With Monaco, Kotliarov2020, Granja2019 references
  - **Silhouette analysis**: Clustering quality metrics

- Compute all metrics from the benchmark table:
  - Accuracy (Acc)
  - F1 score
  - Adjusted Mutual Information (AMI)
  - Adjusted Rand Index (ARI)
  - Silhouette score (reference, predicted, random)

- Create a benchmark results notebook/script that:
  - Runs each analysis
  - Computes metrics
  - Produces a comparison table (v3 results vs. original Supp Table 7 values)
  - Documents any drift with explanations

### Acceptance Criteria
- All 10 datasets download and load successfully with besca v3
- All sig-annot and auto-annot analyses run to completion
- Metrics computed for all rows in Supplementary Table 7
- Results table produced comparing v3 vs. published values
- Drift documented (expected due to dependency version changes in sklearn, scanpy, etc.)
- Benchmark notebook committed to repo (e.g., `benchmarks/supp_table7_replication.py`)

---

## Implementation Plan (Ordered)

### Phase 1: Foundation
1. Create `besca3` branch from `main`
2. Migrate build system to `pyproject.toml`, remove versioneer
3. Upgrade Python to 3.12 and update all dependencies
4. Fix import errors and API breakages from dependency upgrades
5. Remove R dependencies and related code
6. Rewrite `reformat` C binary in Python

### Phase 2: Issue Resolution
7. Fix all 17 bug issues
8. Implement all 9 enhancement issues
9. Address all 6 documentation/meta issues

### Phase 3: CI/CD and Documentation
10. Modernize all GitHub Actions workflows
11. Update Sphinx documentation for v3
12. Create automated gh-pages deployment workflow
13. Build and deploy documentation

### Phase 4: Validation
14. Run existing test suite, fix failures
15. Download all 10 Zenodo datasets
16. Replicate besca_publication_results analyses
17. Compute Supplementary Table 7 benchmark metrics
18. Document results and drift

### Phase 5: Release Preparation
19. Update README.md for v3
20. Final test pass
21. Tag v3.0.0
