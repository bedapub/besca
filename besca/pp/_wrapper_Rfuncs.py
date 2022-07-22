import importlib
import scanpy as sc


def valOutlier(adata, nmads=3, rlib_loc=""):
    """
    Estimates and returns the thresholds to use for gene/cell filtering based on outliers calculated from the deviation to the median QCs. Wrapper function based on 'isOutlier' function of the 'scater' R package.

    Parameters
    ----------
    adata: `AnnData`
        Unfiltered AnnData object of RNA counts.
    nmads: `int`
        Number of median absolute deviation to use as threshold for outlier detection. Lenient NMADS (3 to 5) generally yield the best results.
    rlib_loc: `str`
        R library location that will be added to the default .libPaths() to locate the required packages.

    Returns
    -------
    The estimated parameters to set in the besca workflow considering the QC distribution.
    """

    rpy2_import = importlib.util.find_spec("rpy2")
    if rpy2_import is None:
        raise ImportError("deviance requires rpy2. Install with pip install rpy2")
    import rpy2.robjects as ro
    import anndata2ri
    from scipy.sparse import issparse

    anndata2ri.activate()

    ro.globalenv["rlib_loc"] = rlib_loc
    ro.r(".libPaths(c(rlib_loc, .libPaths()))")
    ro.r("suppressPackageStartupMessages(library(scater))")
    ro.r("suppressPackageStartupMessages(library(Matrix))")
    ro.r("suppressPackageStartupMessages(library(Seurat))")

    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv["dat"] = adata
    ro.globalenv["sym"] = adata.var["SYMBOL"]
    ro.r('seurat_obj = as.Seurat(dat, counts="X", data = NULL)')
    ro.r(
        """"
            if (packageVersion("Seurat") >= "4.0" ) {
            print("Detected Seurat Version >= 4")
            dat =  SingleCellExperiment(assays = list(counts=seurat_obj@assays$originalexp@counts))
            } else {
            print("Detected Seurat Version < 4")
            dat =  SingleCellExperiment(assays = list(counts=seurat_obj@assays$RNA@counts))
            }
            """
    )
    ro.r("rownames(dat) <- sym")
    ro.r(
        """
    valOutlier <- function(dat, nmads = 3){

      mito <- grep('MT-', rownames(dat))
      if (length(mito) == 0){
        mito <- NULL
      } else {
        mito <- list(Mito = mito)
      }

      stats_cells <- perCellQCMetrics(dat, subsets = mito )
      stats_genes <- perCellQCMetrics(t(counts(dat)))

      lower_detected <- as.numeric(attr(isOutlier(stats_cells$detected, nmads = nmads, type = 'lower'), 'thresholds')['lower'])
      if(lower_detected < 0) lower_detected <- 0
      rm_detected <- sum(isOutlier(stats_cells$detected, nmads = nmads, type = 'lower'))

      lower_expressed <- as.numeric(attr(isOutlier(stats_genes$detected, nmads = nmads, type = 'lower'), 'thresholds')['lower'])
      if(lower_expressed < 0) lower_expressed <- 0
      rm_expressed <- sum(isOutlier(stats_genes$detected, nmads = nmads, type = 'lower'))

      lower_sum <- as.numeric(attr(isOutlier(stats_cells$sum, nmads = nmads, type = 'lower'), 'thresholds')['lower'])
      if(lower_sum < 0) lower_sum <- 0
      rm_sum <- sum(isOutlier(stats_cells$sum, nmads = nmads, type = 'lower'))

      higher_detected <- as.numeric(attr(isOutlier(stats_cells$detected, nmads = nmads, type = 'higher'), 'thresholds')['higher'])
      rm_high_detected <- sum(isOutlier(stats_cells$detected, nmads = nmads, type = 'higher'))

      if(!length(mito) == 0){
        max_mito <- as.numeric(attr(isOutlier(stats_cells$subsets_Mito_percent , nmads = nmads, type = 'higher'), 'thresholds')['higher'])/100
        if(max_mito>1) max_mito <- 1
        rm_mito <- sum(isOutlier(stats_cells$subsets_Mito_percent, nmads = nmads, type = 'higher'))
      }

      higher_sum <- as.numeric(attr(isOutlier(stats_cells$sum, nmads = nmads, type = 'higher'), 'thresholds')['higher'])
      if(is.na(higher_sum)) higher_sum <- as.numeric(attr(isOutlier(stats_cells$sum, nmads = nmads, type = 'higher'), 'thresholds')['higher', 1])
      rm_high_sum <- sum(isOutlier(stats_cells$sum, nmads = nmads, type = 'higher'))

      message('Advised parameters based on outliers with ',nmads, ' NMADS:')
      message('standard_min_genes: ', round(lower_detected,2), ', removing ', rm_detected, ' cells')
      message('standard_min_cells: ', round(lower_expressed, 2), ', removing ', rm_expressed, ' genes')
      message('standard_min_counts: ', round(lower_sum,2), ', removing ', rm_sum, ' cells')
      message('standard_n_genes: ', round(higher_detected, 2), ', removing ', rm_high_detected, ' cells')
      if(!length(mito) == 0) {
          message('standard_percent_mito: ', round(max_mito, 2), ', removing ', rm_mito, ' cells')
      } else {
          message('No mitochondrial gene detected.')
          max_mito <- 1
      }
      message('standard_max_counts: ', round(higher_sum, 2), ', removing ', rm_high_sum, ' cells')

      return(c(round(lower_detected,2),
            round(lower_expressed, 2),
            round(lower_sum,2),
            round(higher_detected, 2),
            round(max_mito, 2),
            round(higher_sum, 2)))
    }

     """
    )
    ro.globalenv["nmads"] = nmads
    return ro.r("valOutlier(dat, nmads = nmads)")


def scTransform(adata, hvg=False, n_genes=4000, rlib_loc=""):
    """
    Function to call scTransform normalization or HVG selection from Python. Modified from https://github.com/normjam/benchmark/blob/master/normbench/methods/ad2seurat.py.

    Parameters
    ----------
    adata: `AnnData`
        AnnData object of RNA counts.
    hvg: `boolean`
        Should the hvg method be used (returning a reduced adata object) or the normalization method (returning a normalized adata).
    n_genes: `int`
        Number of hvgs to return if the hvg method is selected. A selection of 4000-5000 generally yields the best results.
    rlib_loc: `str`
        R library location that will be added to the default .libPaths() to locate the required packages.

    Returns
    -------
    returns an AnnData object reduced to the highly-variable genes.
    """

    rpy2_import = importlib.util.find_spec("rpy2")
    if rpy2_import is None:
        raise ImportError("deviance requires rpy2. Install with pip install rpy2")
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri
    import anndata2ri
    from scipy.sparse import issparse

    ro.globalenv["rlib_loc"] = rlib_loc
    ro.r(".libPaths(c(rlib_loc, .libPaths()))")
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(scater))")
    anndata2ri.activate()

    sc.pp.filter_genes(adata, min_cells=5)

    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv["adata"] = adata

    ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')
    if hvg:
        numpy2ri.activate()
        ro.globalenv["n_genes"] = n_genes
        print("Reducing the data to", n_genes, "variable genes.")
        ro.r(
            "res <- SCTransform(object=seurat_obj, return.only.var.genes = TRUE, do.correct.umi = FALSE, variable.features.n = n_genes)"
        )
        hvgs_r = ro.r("res@assays$SCT@var.features")
        adata = adata[:, list(hvgs_r)]
        adata.var["highly_variable"] = True
        return adata
    else:
        ro.r(
            "res <- SCTransform(object=seurat_obj, return.only.var.genes = FALSE, do.correct.umi = FALSE)"
        )

        norm_x = ro.r("res@assays$SCT@scale.data").T

        adata.layers["counts"] = norm_x
        adata.raw = adata
        return adata
