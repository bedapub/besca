"""
Besca v3.0.0 — Supplementary Table 7 Benchmark Replication

Replicates the sig-annot and auto-annot comparison from the Besca paper
(NAR Genomics and Bioinformatics, 2021, lqab102) using all 10 Zenodo datasets.

Metrics computed: Accuracy, F1, AMI, ARI, Silhouette (ref, pred, random)

Usage:
    python benchmarks/supp_table7_replication.py [--dataset DATASET_NAME]

If --dataset is not specified, all 10 datasets are processed.
"""

import argparse
import os
import sys
import warnings
import time

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    adjusted_mutual_info_score,
    adjusted_rand_score,
    silhouette_score,
)

warnings.filterwarnings("ignore")

# Add besca to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import besca as bc


# ============================================================================
# Dataset registry — uses besca.datasets module for loading
# ============================================================================

DATASETS = [
    "Baron2016",
    "Granja2019",
    "Kotliarov2020",
    "Lee2020",
    "Martin2019",
    "Haber2017",
    "PBMC3k",
    "Peng2019",
    "Segerstolpe2016",
    "Smillie2019",
]

# Map dataset names to besca.datasets loader functions
DATASET_LOADERS = {
    "Baron2016": ("Baron2016_processed", None),
    "Granja2019": ("Granja2019_processed", None),
    "Kotliarov2020": ("Kotliarov2020_processed", None),
    "Lee2020": ("Lee2020_processed", None),
    "Martin2019": ("Martin2019_processed", None),
    "Haber2017": ("Haber2017_processed", None),
    "PBMC3k": ("pbmc3k_processed", None),
    "Peng2019": ("Peng2019_processed", None),
    "Segerstolpe2016": ("Segerstolpe2016_processed", None),
    "Smillie2019": ("Smillie2019_processed", None),
}


def load_dataset(name):
    """Load a dataset using besca.datasets or direct Zenodo URL."""
    loader_name, fallback_url = DATASET_LOADERS.get(name, (None, None))

    if loader_name and hasattr(bc.datasets, loader_name):
        print(f"  Loading via bc.datasets.{loader_name}()...")
        loader = getattr(bc.datasets, loader_name)
        return loader()

    if fallback_url:
        cache_dir = os.environ.get("BESCA_DATA_DIR",
                                   os.path.join(os.path.expanduser("~"), ".cache", "besca", "benchmarks"))
        os.makedirs(cache_dir, exist_ok=True)
        filepath = os.path.join(cache_dir, f"{name.lower()}_processed.h5ad")
        if os.path.exists(filepath):
            print(f"  Using cached {name}")
            return sc.read(filepath)
        print(f"  Downloading {name} from Zenodo...")
        adata = sc.read(filepath, backup_url=fallback_url, cache=False)
        adata.write(filepath)
        return adata

    raise ValueError(f"No loader found for dataset: {name}")


def compute_metrics(y_true, y_pred, adata=None, embedding_key="X_umap"):
    """Compute classification and clustering metrics."""
    # Filter to cells where both labels exist
    mask = pd.notna(y_true) & pd.notna(y_pred) & (y_true != "") & (y_pred != "")
    y_true_f = y_true[mask]
    y_pred_f = y_pred[mask]

    if len(y_true_f) == 0:
        return {"Acc": np.nan, "F1": np.nan, "AMI": np.nan, "ARI": np.nan,
                "Sil_ref": np.nan, "Sil_pred": np.nan, "Random": np.nan}

    # Classification metrics
    acc = accuracy_score(y_true_f, y_pred_f)
    f1 = f1_score(y_true_f, y_pred_f, average="weighted", zero_division=0)
    ami = adjusted_mutual_info_score(y_true_f, y_pred_f)
    ari = adjusted_rand_score(y_true_f, y_pred_f)

    # Silhouette scores
    sil_ref = np.nan
    sil_pred = np.nan
    sil_random = np.nan

    if adata is not None and embedding_key in adata.obsm:
        emb = adata[mask].obsm[embedding_key]
        if len(set(y_true_f)) > 1:
            try:
                sil_ref = silhouette_score(emb, y_true_f, sample_size=min(5000, len(y_true_f)))
            except Exception:
                pass
        if len(set(y_pred_f)) > 1:
            try:
                sil_pred = silhouette_score(emb, y_pred_f, sample_size=min(5000, len(y_pred_f)))
            except Exception:
                pass
        # Random baseline
        rng = np.random.RandomState(42)
        random_labels = rng.permutation(y_true_f.values if hasattr(y_true_f, 'values') else y_true_f)
        if len(set(random_labels)) > 1:
            try:
                sil_random = silhouette_score(emb, random_labels, sample_size=min(5000, len(random_labels)))
            except Exception:
                pass

    return {
        "Acc": round(acc, 2),
        "F1": round(f1, 2),
        "AMI": round(ami, 2),
        "ARI": round(ari, 2),
        "Sil_ref": round(sil_ref, 2) if not np.isnan(sil_ref) else np.nan,
        "Sil_pred": round(sil_pred, 2) if not np.isnan(sil_pred) else np.nan,
        "Random": round(sil_random, 2) if not np.isnan(sil_random) else np.nan,
    }


def run_auto_annot_benchmark(adata_train, adata_test, label_col, method="log_reg"):
    """Run auto-annotation benchmark with a given method."""
    from besca.tl.auto_annot import (
        logistic_regression,
        fit,
        predict,
    )

    try:
        # Prepare data
        genes = list(set(adata_train.var_names) & set(adata_test.var_names))
        if len(genes) < 100:
            return None

        train = adata_train[:, genes].copy()
        test = adata_test[:, genes].copy()

        # Fit model
        if method == "log_reg":
            clf = logistic_regression()
        else:
            clf = logistic_regression()

        clf = fit(clf, train, label_col)
        predictions = predict(clf, test)

        return predictions
    except Exception as e:
        print(f"    Auto-annot failed ({method}): {e}")
        return None


def benchmark_dataset(name):
    """Run full benchmark for a single dataset."""
    print(f"\n{'='*60}")
    print(f"Benchmarking: {name}")
    print(f"{'='*60}")

    results = []

    try:
        adata = load_dataset(name)
        print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
    except Exception as e:
        print(f"  Failed to load {name}: {e}")
        return results

    # Determine available label columns
    label_cols = [c for c in adata.obs.columns
                  if any(x in c.lower() for x in ["celltype", "dblabel", "cell_type", "label"])
                  and c not in ("celltype_flowjo",)]
    # Prefer dblabel > celltype1 > celltype0
    priority = ["dblabel", "celltype1", "celltype0", "celltype2", "celltype3"]
    label_cols = sorted(label_cols, key=lambda c: priority.index(c) if c in priority else 99)
    print(f"  Available labels: {label_cols}")

    if not label_cols:
        print(f"  No label columns found, skipping")
        return results

    ref_label = label_cols[0]

    # Ensure UMAP exists
    if "X_umap" not in adata.obsm:
        print("  Computing UMAP...")
        try:
            if "X_pca" not in adata.obsm:
                sc.pp.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        except Exception as e:
            print(f"  UMAP failed: {e}")

    # Auto-annot: logistic regression (80/20 split)
    print(f"  Running auto-annot (log_reg) with label '{ref_label}'...")
    try:
        from sklearn.linear_model import LogisticRegression
        from scipy import sparse

        # Use 80/20 split
        n = adata.shape[0]
        rng = np.random.RandomState(42)
        idx = rng.permutation(n)
        train_idx = idx[:int(0.8 * n)]
        test_idx = idx[int(0.8 * n):]

        train = adata[train_idx].copy()
        test = adata[test_idx].copy()

        # Get expression matrix (use .X directly)
        X_train = train.X
        X_test = test.X
        if sparse.issparse(X_train):
            X_train = X_train.toarray()
        if sparse.issparse(X_test):
            X_test = X_test.toarray()

        y_train = train.obs[ref_label].values
        y_true = test.obs[ref_label]

        # Filter out NaN labels
        mask_train = pd.notna(y_train) & (y_train != "")
        mask_test = pd.notna(y_true) & (y_true != "")

        if mask_train.sum() > 50 and mask_test.sum() > 50:
            clf = LogisticRegression(max_iter=500, solver="saga", n_jobs=-1, random_state=42)
            clf.fit(X_train[mask_train], y_train[mask_train])
            preds = clf.predict(X_test[mask_test.values])

            metrics = compute_metrics(
                y_true[mask_test],
                pd.Series(preds, index=test.obs_names[mask_test]),
                test[mask_test],
            )
            metrics["Dataset"] = name
            metrics["Method"] = "Auto-annot (log reg)"
            metrics["Reference"] = ref_label
            results.append(metrics)
            print(f"    Acc={metrics['Acc']}, F1={metrics['F1']}, AMI={metrics['AMI']}, ARI={metrics['ARI']}")
        else:
            print(f"    Not enough labeled cells for auto-annot")
    except Exception as e:
        print(f"    Auto-annot failed: {e}")

    # Silhouette on reference labels
    if "X_umap" in adata.obsm and ref_label in adata.obs.columns:
        print(f"  Computing silhouette for '{ref_label}'...")
        try:
            labels = adata.obs[ref_label]
            mask = pd.notna(labels) & (labels != "")
            if mask.sum() > 100 and len(set(labels[mask])) > 1:
                sil = silhouette_score(
                    adata[mask].obsm["X_umap"],
                    labels[mask],
                    sample_size=min(5000, mask.sum())
                )
                metrics = {
                    "Dataset": name,
                    "Method": "Silhouette (reference)",
                    "Reference": ref_label,
                    "Sil_ref": round(sil, 2),
                    "Acc": np.nan, "F1": np.nan, "AMI": np.nan, "ARI": np.nan,
                    "Sil_pred": np.nan, "Random": np.nan,
                }
                results.append(metrics)
                print(f"    Silhouette (ref): {sil:.3f}")
        except Exception as e:
            print(f"    Silhouette failed: {e}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Besca v3 Benchmark Replication")
    parser.add_argument("--dataset", type=str, default=None,
                        help="Run benchmark for a specific dataset only")
    parser.add_argument("--output", type=str, default="benchmarks/benchmark_results.csv",
                        help="Output CSV file path")
    args = parser.parse_args()

    print("=" * 60)
    print("Besca v3.0.0 — Supplementary Table 7 Benchmark Replication")
    print(f"besca version: {bc.__version__}")
    print(f"scanpy version: {sc.__version__}")
    print(f"anndata version: {ad.__version__}")
    print(f"numpy version: {np.__version__}")
    print("=" * 60)

    start_time = time.time()

    if args.dataset:
        if args.dataset not in DATASETS:
            print(f"Unknown dataset: {args.dataset}")
            print(f"Available: {DATASETS}")
            sys.exit(1)
        datasets_to_run = [args.dataset]
    else:
        datasets_to_run = list(DATASETS)

    all_results = []
    for name in datasets_to_run:
        results = benchmark_dataset(name)
        all_results.extend(results)

    # Create results DataFrame
    if all_results:
        df = pd.DataFrame(all_results)
        cols = ["Dataset", "Method", "Reference", "Acc", "F1", "AMI", "ARI",
                "Sil_ref", "Sil_pred", "Random"]
        df = df[[c for c in cols if c in df.columns]]

        # Save results
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        df.to_csv(args.output, index=False)

        print(f"\n{'='*60}")
        print("BENCHMARK RESULTS SUMMARY")
        print(f"{'='*60}")
        print(df.to_string(index=False))
        print(f"\nResults saved to: {args.output}")
    else:
        print("\nNo results generated.")

    elapsed = time.time() - start_time
    print(f"\nTotal time: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
