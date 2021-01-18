
import os
from scanpy import read_csv as sc_read_csv
import importlib
from scanpy.tools import pca as sc_pca

def maxLikGlobalDimEst(adata, k = 20, nrpcs=50, rlib_loc = ''):
    """
    Estimates the intrinsic dimensionality of the data, based on the 'maxLikGlobalDimEst' function of the 'intrinsicDimension' R package.
    
    Parameters
    ----------
    adata: `AnnData`
        AnnData object of RNA counts.
    k: `
        Number of neighbours to use in the 'maxLikGlobalDimEst'. Choosing k between 10 and 20 generally yields the best results. 
    nrpcs:
        Number of PCs to compute initially before estimating the dimensionality. Consider increasing it for very high dimensional data. 
    rlib_loc: `str`
        R library location that will be added to the default .libPaths() to locate the required packages. 
  
    Returns
    -------
    Returns the estimated intrinsic dimensionality of the data that can be used for graph clustering. 
    """
    
    rpy2_import = importlib.util.find_spec('rpy2')
    if rpy2_import is None:
        raise ImportError(
            "maxLikGlobalDimEst requires rpy2. Install with pip install rpy2")
    import rpy2.robjects as ro
    import anndata2ri
    from scipy.sparse import issparse

    ro.globalenv['rlib_loc'] = rlib_loc
    ro.r('.libPaths(c(.libPaths(), rlib_loc))')
    ro.r('suppressPackageStartupMessages(library(intrinsicDimension))')

    random_state = 0
    print('Using random_state = 0 for all the following calculations')
    sc_pca(adata, svd_solver='arpack', random_state=0, n_comps=nrpcs)
    adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    
    ro.globalenv['pcs'] = adata.obsm['X_pca']
    ro.globalenv['k'] = k
    ro.r('n <- maxLikGlobalDimEst(as.matrix(pcs), k=k, unbiased=TRUE)')

    ro.r('message("Estimated dimensionality: ", round(n$dim.est))')

    n_dimest = ro.r('round(n$dim.est)')

    return int(n_dimest[0])



def deviance(adata, n_genes = 4000, rlib_loc = ''):
    """
    Wrapper of the 'deviance' method of highly-variable gene selection, included in the 'scry' R package.  
    
    Parameters
    ----------
    adata: `AnnData`
        AnnData object of RNA counts.
    n_genes: `int`
        Number of highly-variable genes to return. A selection of 4000-5000 generally yields the best results. 
    rlib_loc: `str`
        R library location that will be added to the default .libPaths() to locate the required packages. 
  
    Returns
    -------
    returns an AnnData object reduced to the highly-variable genes. 
    """
    rpy2_import = importlib.util.find_spec('rpy2')
    if rpy2_import is None:
        raise ImportError(
            "deviance requires rpy2. Install with pip install rpy2")
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    import anndata2ri
    from scipy.sparse import issparse
    
    anndata2ri.activate()
    
    ro.globalenv['rlib_loc'] = rlib_loc
    ro.r('.libPaths(c(.libPaths(), rlib_loc))')
    ro.r('suppressPackageStartupMessages(library(scry))')
    ro.r('suppressPackageStartupMessages(library(Seurat))')
    
    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()
    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv['adata'] = adata
    ro.globalenv['n'] = n_genes
    print('Reducing the data to', n_genes, 'variable genes.')
    ro.globalenv['rownam'] = adata.var.index
    ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')
    ro.r('adata <- t(as.matrix(seurat_obj@assays$RNA@counts))')
    ro.r('out <- devianceFeatureSelection(adata)')
    ro.r('out <- sort(devianceFeatureSelection(adata),decreasing = TRUE)[1:n] ')
    hvgs_r =  ro.r('rownam[order(out, decreasing = TRUE)][1:n]')
    adata = adata[:,list(hvgs_r)]
    adata.var['highly_variable'] = True
    
    return adata


def dsb_normalize(adata_prot, raw_path, ana_path, rlib_loc = '',  example_dataset = False, hto = False,  numi_min = 2, numi_max = 3.5):
    """Performs DSB normalization. If isotypes are present among the Ab, please make sure that the relevant Ab have 'isotype' in their names. The function also generate a QC plot when negative cells are imputed from UMI threshold. Please have a look at it and eventually adapt the numi_min and numi_max. It is highly advised to use this function if hto/ isotypes are available as they lead to higher-confidence negative droplets. The function is a wrapper adapter from https://github.com/niaid/dsb. 
    
    Parameters
    ----------
    adata_prot: `AnnData`
      AnnData object of protein counts.
    raw_path: `str` 
        Path to the 'raw' folder. 
    ana_path: `str`
        Path to the 'citeseqDSB' analysis folder. Default should be of form 'analyzed/ANALYSIS_NAME/citeseqDSB'
    rlib_loc: `str`
        R library location that will be added to the default .libPaths() to locate the required packages. 
    example_dataset: `bool`  
        Logical, whether example_dataset is being used or not. 
    hto: `str`
        List of string, either 'Negative' or 'Positive' for each cell, corresponding to the result of the HTO demultiplexing. NaN if the information is not provided. 
    numi_min: `int`
        Minimum log10 RNA count per cell to use as a threshold to select the negative droplets if HTOs are not given. 
    numi_max: `int`
        Maximum log10 RNA count per cell to use as a threshold to select the negative droplets if HTOs are not given. 
    

    Returns
    -------
    returns an AnnData object with DSB-normalized counts. 
    """    
    
    rpy2_import = importlib.util.find_spec('rpy2')
    if rpy2_import is None:
        raise ImportError(
            "dsb_normalize requires rpy2. Install with pip install rpy2")
    from rpy2.robjects.packages import importr
    import rpy2.robjects as ro
    import anndata2ri
    from scipy.sparse import issparse
    
    ro.globalenv['rlib_loc'] = rlib_loc
    ro.r('.libPaths(c(.libPaths(), rlib_loc))')
    ro.r('suppressPackageStartupMessages(library(vctrs))')
    ro.r('suppressPackageStartupMessages(library(ggplot2))')
    ro.r('suppressPackageStartupMessages(library(patchwork))')
    ro.r('suppressPackageStartupMessages(library(dsb))')
    ro.r('suppressPackageStartupMessages(library(tidyverse))')
    ro.r('suppressPackageStartupMessages(library(magrittr))')
    ro.r('suppressPackageStartupMessages(library(data.table))')
    ro.r('suppressPackageStartupMessages(library(Matrix))')
    ro.r('suppressPackageStartupMessages(library(DropletUtils))')
    ro.r('suppressPackageStartupMessages(library(readr))')
        
    ro.r('''
    do_dsb <- function(raw_path, ana_path, example_dataset = FALSE, hto = NA, numi_min = 2, numi_max = 3.5 ){
    message("Reading input data...")
    if (example_dataset){
        mat <- t(read.csv(file.path(raw_path,'mtx', 'adata_raw_x.csv'), row.names = 1))
        mat_prot <- t(read.csv(file.path(raw_path,'mtx_prot', 'X.csv'), header = FALSE))
        rownames(mat_prot) <- read.csv(file.path(raw_path,'mtx_prot', 'var.csv'))[,1]
        colnames(mat_prot) <- read.csv(file.path(raw_path,'mtx_prot', 'obs.csv'), row.names = NULL)[,1]
        cells2keep <- read.csv(file.path(raw_path,'mtx', 'filt_cells.csv'))$col1
    } else {
        mat <- read10xCounts(paste0(raw_path))
        colnames(mat) <- mat$Barcode
        gene_info <- read.csv(paste0(raw_path,"/genes.tsv"), sep = "\t", header = FALSE)
        mat_prot <- counts(mat[gene_info[,3] == "Antibody Capture"])
        mat <- mat[gene_info[,3] == "Gene Expression"]
        cells2keep <- read.csv(paste0(raw_path, "/filt_cells.csv" ))$col1
        mat <- as.matrix(counts(mat))
        mat_prot <- as.matrix(mat_prot)
    }
  
    if(length(hto) > 1) message('Using HTO to define negative droplets.') else hto <- NA
  
    # ------------------------
  
    # calculate metadata 
    rna_size = log10(Matrix::colSums(mat))
    prot_size = log10(Matrix::colSums(mat_prot))
    ngene = Matrix::colSums(mat > 0)
    mtgene = grep(pattern = "^MT-", rownames(mat), value = TRUE)
    propmt = Matrix::colSums(mat[mtgene, ]) / Matrix::colSums(mat)
    md = as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
    md$bc = rownames(md)
  
    # define a vector of cell-containing droplet barcodes 
    cells_mtx_rawprot = mat_prot[ , cells2keep] %>% as.matrix()
  
    if (is.na(hto)) {
      message("Defining negative droplets based on nUMI thresholds.")
      mat <- mat[ ,apply(mat, 2, function(x) sum(x != 0) > 0) | (colnames(mat) %in% cells2keep)]
      mat <- mat[apply(mat, 1, function(x) sum(x != 0) > 0), ]
      # Apply basic filtering recommended by DSB prior to selecting negative/positive drops
      # but by keeping the cells that were not filtered out by besca filtering
      md = md[colnames(mat), ]
      md$nCount_RNA <- colSums(mat)
      md$log10umi = log10(md$nCount_RNA  + 1) 


      # define negative based on an mRNA threshold to define neg and positive cells (see details below)
      neg_drops <- md$bc[(md$log10umi > numi_min & md$log10umi < numi_max) | ((md$bc %in% cells2keep))]
      if(length(neg_drops) == 0) stop('Not enough negative drops selected. Please increase numi_min/ decrease numi_max')
      negative_mtx_rawprot = mat_prot[ , neg_drops] %>%  as.matrix()

      # Plot
      ndrop = dim(md)[1]

      hist_attr = list( theme_bw() , theme(text = element_text(size = 8)) , geom_density(fill = "#3e8ede") )
      p1 = ggplot(md[c(colnames(cells_mtx_rawprot),neg_drops), ], aes(x = log10(nCount_RNA + 1 ) )) +
        hist_attr +
        ggtitle(paste0( " raw_feature_bc_matrix: ", ndrop, " droplets")) + 
        xlim(0, NA) +
        geom_vline(xintercept = c(numi_max, numi_min ),   linetype = "dashed") +
        annotate("text", x = 1, y=1.5, label = " region 1: \n void of data ") +
        annotate("text", x = numi_min, y=2, label = paste0(" region 2: \n",length(neg_drops), " background drops \n define 'empty_drop_matrix' \n with these drops ")) +
        annotate("text", x = (numi_min +2), y=2, label = " region 3: \n cell containing droplets \n zomed in on next plot")

      p2 = ggplot(md[log10(md$nCount_RNA + 1) > numi_max, ], aes(x = log10(nCount_RNA + 1 ) )) +
        hist_attr +
        ggtitle(paste0(ncol(cells_mtx_rawprot)," drops containing cells "))
      p3 = cowplot::plot_grid( p1 , p2 )
      message('Saving DSB QC plot in \n', file.path(ana_path, 'DSB_qc.png'))
      ggsave(file.path(ana_path, 'DSB_qc.png'), p3, width = 12, height = 6)


    } else {
      negative_mtx_rawprot = mat_prot[ , hto == 'Negative'] %>%  as.matrix()
    }
    
    message(paste("Using", ncol(negative_mtx_rawprot), "cells as negative droplets"))
    message(paste("Using", ncol(cells_mtx_rawprot), "cells as positive droplets"))
    
    # run DSB normalization
    message("DSB normalization")
    iso <- grep("[iI]sotype", rownames(mat_prot))
    if(length(iso) == 0) {
      message("No isotype found: removing background without isotype control")
      mtx = DSBNormalizeProtein(cell_protein_matrix = cells_mtx_rawprot,
                                empty_drop_matrix = negative_mtx_rawprot,
                                denoise.counts = TRUE, use.isotype.control = FALSE, 
                                isotype.control.name.vec = FALSE)
    } else {
      message(paste(length(iso), "isotype(s) found: using them as control for the normalization"))
      isotypes = rownames(mat_prot)[iso]
      mtx = DSBNormalizeProtein(cell_protein_matrix = cells_mtx_rawprot,
                                empty_drop_matrix = negative_mtx_rawprot,
                                denoise.counts = TRUE, use.isotype.control = TRUE, 
                                isotype.control.name.vec = isotypes)
    }
    
    message("Saving file with normalized proteins in")
    message(file.path(ana_path, 'citeseq', 'normalized_counts', "dsb_norm_matrix.csv"))
    # Done! you can also save the resulting normalized matrix for integration with scanpy etc
    write_delim(as.data.frame(t(mtx)), file = file.path(ana_path, "dsb_norm_matrix.csv"),delim = "," )
    
    
    }   
        ''')
    r_do_dsb = ro.r['do_dsb']
    r_do_dsb(raw_path = raw_path,
             ana_path = ana_path,
             example_dataset = example_dataset,
             hto = hto,  
             numi_min = numi_min,
             numi_max = numi_max)
    a = sc_read_csv(os.path.join(ana_path, 'dsb_norm_matrix.csv'))
    a.obs = adata_prot.obs
    a.var = adata_prot.var
    a.raw = adata_prot.copy()
    a.layers['counts'] = adata_prot.layers["counts"]
    adata_prot = a

    return(adata_prot)
