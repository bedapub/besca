# this contains wrapper functions to export data into the standard format for the standard pipeline

import logging
import os
# import other modules
import sys
from os import makedirs
from os.path import join
from time import time
import bbknn
import scipy
from matplotlib.pyplot import subplots
from numpy import cumsum
from pandas import DataFrame
# import scanpy functions
from scanpy.api.pl import highly_variable_genes as pl_highly_variable_genes
from scanpy.api.pl import pca as sc_pl_pca
from scanpy.api.pl import umap as sc_pl_umap
from scanpy.api.pp import highly_variable_genes as sc_highly_variable_genes
from scanpy.api.pp import log1p, neighbors, normalize_per_cell
from scanpy.api.pp import regress_out as sc_regress_out
from scanpy.api.pp import scale as sc_scale
from scanpy.api.tl import leiden as sc_leiden
from scanpy.api.tl import louvain as sc_louvain
from scanpy.api.tl import pca as sc_pca
from scanpy.api.tl import rank_genes_groups
from scanpy.api.tl import umap as sc_umap


# import other besca functions
from .. import _logging as logs
from ..export._export import labeling, labeling_info
from ..Import._read import read_mtx
from ..pp._filtering import filter
from ..pp._normalization import normalize_geometric
from ..tl.bcor import batch_correct, postprocess_mnnpy
from ._FAIR_export import (export_clr, export_clustering, export_cp10k,
                           export_metadata, export_rank, export_regressedOut)


def setup(results_folder,
          analysis_name,
          labeling_name,
          labeling_to_use,
          log_file,
          version,
          root_path,
          species,
          batch_to_correct,
          standard_min_genes,
          standard_min_cells,
          standard_min_counts,
          standard_n_genes,
          standard_percent_mito,
          standard_max_counts):
    '''
    This function is only intended for use in the standardpipeline! As such it has no input variables
    since all of these variables have been already defined in the standard pipeline. It does not appear
    in bescas documentation since it should not be modified.
    '''
    start = time()
    # generate the necessary directories
    makedirs(results_folder, exist_ok=True)
    makedirs(join(results_folder, 'figures'), exist_ok=True)
    makedirs(join(results_folder, 'labelings'), exist_ok=True)
    makedirs(join(results_folder, 'labelings', 'leiden'), exist_ok=True)
    makedirs(join(results_folder, 'labelings', 'louvain'), exist_ok=True)
    if labeling_to_use != 'None':
        makedirs(join(results_folder, 'labelings',
                      labeling_name), exist_ok=True)
    makedirs(join(results_folder, 'normalized_counts'), exist_ok=True)
    makedirs(join(results_folder, 'normalized_counts', 'cp10k'), exist_ok=True)
    makedirs(join(results_folder, 'normalized_counts',
                  'regressedOut'), exist_ok=True)
    print('all output directories created successfully')

    if os.path.exists(log_file):
        os.remove(log_file)
        open(log_file, "x")
    else:
        open(log_file, "x")

    # setup logging
    logs.initialize_logger(log_file)

    # initialize log file
    logs.initialize_log_file(analysis_name,
                             root_path,
                             species,
                             batch_to_correct,
                             standard_min_genes,
                             standard_min_cells,
                             standard_min_counts,
                             standard_n_genes,
                             standard_percent_mito,
                             standard_max_counts,
                             version)

    # output feedback to logfile
    logging.info('\tTime for creating all output directories and setting up logging: ' +
                 str(round(time()-start, 3))+'s')


def setup_citeseq(results_folder,
                  labeling_name,
                  labeling_to_use, ):
    '''This function generates the required file structure to save the citeseq data
    '''

    # time file creation
    logging.info(
        'CiteSeq Values present in  Dataset. Setting up analysis to utilize CiteSeq data.')
    start = time()

    # define file paths
    results_folder_citeseq = os.path.join(results_folder, 'citeseq')
    results_folder_merged = os.path.join(results_folder, 'citeseq_merged')

    # generate folder structure
    makedirs(results_folder_citeseq, exist_ok=True)
    makedirs(results_folder_merged, exist_ok=True)
    makedirs(join(results_folder_citeseq, 'figures'), exist_ok=True)
    makedirs(join(results_folder_merged, 'figures'), exist_ok=True)
    makedirs(join(results_folder_citeseq, 'labelings'), exist_ok=True)
    makedirs(join(results_folder_citeseq, 'labelings', 'leiden'), exist_ok=True)
    makedirs(join(results_folder_citeseq, 'labelings', 'louvain'), exist_ok=True)
    if labeling_to_use != 'None':
        makedirs(join(results_folder_citeseq, 'labelings',
                      labeling_name), exist_ok=True)
    makedirs(join(results_folder_citeseq, 'normalized_counts'), exist_ok=True)

    # generate log message
    print('all output directories for citeseq data created successfully')
    logging.info('\tTime for creating all citeseq output directories: ' +
                 str(round(time()-start, 3))+'s')


def read_matrix(root_path,
                citeseq=None,
                annotation=True,
                use_genes='SYMBOL',
                species='human'):
    """Read matrix file as expected for the standard workflow.
    ----------
    root_path: `str`
       root path of the analysis. Expected in this folder a raw folder containing
       containg the matrix.mtx, genes.tsv,
        barcodes.tsv and if applicable metadata.tsv
    annotation: `bool` (default = True)
        boolian identifier if an annotation file is also located in the folder
        and should be added to the AnnData object
    use_genes: `str`
        either SYMBOL or ENSEMBL. Other genenames are not yet supported.
    species: `str` | default = 'human'
        string specifying the species, only needs to be used when no Gene Symbols
        are supplied and you only have the ENSEMBLE gene ids to perform a lookup.
    citeseq: 'gex_only' or 'citeseq_only' or None | default = None
        string indicating if only gene expression values (gex_only) or only protein
        expression values ('citeseq_only') or everything is read if None is specified 
    Returns
    -------
    returns an AnnData object
  """
    start = time()
    input_path = join(root_path, 'raw')

    adata = read_mtx(input_path, citeseq=citeseq,
                     annotation=annotation,
                     use_genes=use_genes,
                     species=species)
    logging.info('After input: ' +
                 str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' genes')
    logging.info("\tTime for reading data: "+str(round(time()-start, 3))+'s')

    return(adata)


def filtering_cells_genes_min(adata, standard_min_cells, standard_min_genes, standard_min_counts):
    # record start time
    start = time()
    # perform first filtering
    adata = filter(adata, annotation_type='SYMBOL', min_cells=standard_min_cells,
                   min_genes=standard_min_genes, min_counts=standard_min_counts)
    # generate logs
    logging.info('After filtering for minimum number of cells and minimum number of expressed genes: ' +
                 str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' genes')
    logging.info("\tTime for filtering: "+str(round(time()-start, 3))+'s')
    return(adata)


def filtering_mito_genes_max(adata, standard_percent_mito, standard_n_genes, standard_max_counts):
    # record start time
    start = time()
    # perform first filtering
    adata = filter(adata, annotation_type='SYMBOL', max_mito=standard_percent_mito,
                   max_genes=standard_n_genes, max_counts=standard_max_counts)
    # generate logs
    logging.info('After filtering for maximum number of expressed genes and max percent mito: ' +
                 str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' genes')
    logging.info("\tTime for filtering: "+str(round(time()-start, 3))+'s')

    return adata


def per_cell_normalize(adata, results_folder):
    # get start time
    start = time()
    # normalize per cell
    # already normalize BEFORE saving "raw" - as recommended in the scanpy tutorial
    normalize_per_cell(adata, counts_per_cell_after=1e4)
    print('adata normalized per cell')

    # keep raw copy
    adata.raw = log1p(adata, copy=True)
    print('log1p values saved into adata.raw')

    # make log entries
    logging.info('Per cell normalization completed successfully.')
    logging.info("\tTime for per-cell normalization: " +
                 str(round(time()-start, 3))+'s')

    # export to file
    start = time()
    export_cp10k(adata, basepath=results_folder)

    logging.info('cp10k values exported to file.')
    logging.info("\tTime for cp10k export: "+str(round(time()-start, 3))+'s')

    return(adata)


def clr_normalize(adata, results_folder):
    """Perform clr normalization.
    """

    # get start time
    start = time()

    # normalize per cell
    # this also already applies log! Is not taken seperately
    # already normalize BEFORE saving "raw" - as recommended in the scanpy tutorial
    normalize_geometric(adata)
    print('clr normalization applied to adata')

    # keep raw copy
    adata.raw = adata.copy()
    print('normalized values saved into adata.raw')

    # make log entries
    logging.info('CLR normalization completed successfully.')
    logging.info("\tTime for CLR normalization: " +
                 str(round(time()-start, 3))+'s')

    # export to file
    start = time()
    export_clr(adata, basepath=results_folder)

    logging.info('CLR values exported to file.')
    logging.info("\tTime for CLR export: "+str(round(time()-start, 3))+'s')

    return(adata)


def highly_variable_genes(adata, batch_key=None):
    start = time()

    # take log1p
    log1p(adata)
    print('log1p taken of adata')

    filter_result = sc_highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, inplace=False, batch_key=batch_key)
    pl_highly_variable_genes(filter_result, save='.hvg.png', show=True)

    adata = adata[:, filter_result.highly_variable == True]

    # logging
    logging.info('After feature selection of highly variable genes: ' +
                 str(adata.shape[0])+' cells, ' + str(adata.shape[1])+' genes')
    logging.info('\tTime for feature selection: ' +
                 str(round(time()-start, 3))+'s')

    return(adata)


def regress_out(adata, results_folder):
    start = time()
    # regressout
    sc_regress_out(adata, ['n_counts', 'percent_mito'])
    print("'n_counts' and 'percent_mito' regressed out")
    sc_scale(adata, max_value=10)
    print('adata scaled with max_value set to 10')
    # logging
    logging.info(
        "Regression steps completed. 'n_counts' and 'percent_mito' regressed out. adata was log-normalized and scaled.")
    logging.info('\tTime for regression steps: ' +
                 str(round(time()-start, 3))+'s')
    return(adata)


def batch_correction(adata, batch_to_correct):
    start = time()
    # get number of batches
    batches = list(set(adata.obs[batch_to_correct]))
    batchNum = len(batches)

    # perform batch correction
    bdata = batch_correct(adata, batch_to_correct)
    # perform post processing to regenerate original raw
    adata = postprocess_mnnpy(adata, bdata)
    print('postprocessing performed. adata contains original .raw')

    logging.info("Batch correction on '"+batch_to_correct +
                 "' with '"+str(batchNum)+"' categories completed")
    logging.info('\tTime for batch correction: ' +
                 str(round(time()-start, 3))+'s')
    # return new AnnData object
    return(adata)


def pca_neighbors_umap(adata, results_folder, nrpcs=50, nrpcs_neigh=None, nrneigh=10, method='NULL'):
    '''
    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    results_folder: `str`
        path to the results folder 
    nrpcs: int | nrpcs = 50
        number of principle components to calculate
    nrpcs_neigh: int | nrpcs_neigh = 50
        number of principle components to use for nearest neighbor calculation. 
        When set to None the number is chosen automatically. For .n_vars < 50, .X is used, otherwise ‘X_pca’ is used with 50 components.
    nrneigh: int | nrpcs = None
        number of principle components to calculate
    method: `str`
        Method for nearest neighbor calculation.  Can be set to 'NULL' or bbknn
    '''
    start = time()
    random_state = 0
    print('Using random_state = 0 for all the following calculations')

    # PCA
    sc_pca(adata, svd_solver='arpack', random_state=random_state, n_comps=nrpcs)
    adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    print("PCA calculated using svd_solver = 'arpack'. PCA multiplied by -1 to match Seurat output.")

    # generate plot of PCA
    fig, (ax1, ax2) = subplots(ncols=2, nrows=1)
    fig.set_figwidth(12)
    fig.set_figheight(6)
    fig.tight_layout(pad=4.5)

    cumulative_variance = cumsum(adata.uns['pca']['variance_ratio'])
    x = list(range(nrpcs))
    data = DataFrame({'x': x, 'y': cumulative_variance})

    ax1.scatter(x=x, y=cumulative_variance)
    ax1.set_ylabel('cumulative explained variance')
    ax1.set_xlabel('PCA components')
    ax1.set_title('cumulative explained variance (as ratio)')

    sc_pl_pca(adata, ax=ax2)
    fig.savefig(join(results_folder, 'figures', 'PCA.png'))

    # display(fig)

    # neighbors
    if(method == 'bbknn'):
        if('batch' in adata.obs.columns):
            bbknn.bbknn(adata)
        else:
            sys.exit(
                'bbknn correction requires a column "batch" in the observations.')
    else:
        neighbors(adata, n_neighbors=nrneigh,
                  random_state=random_state, n_pcs=nrpcs_neigh)
        print('Nearest neighbors calculated with n_neighbors = '+str(nrneigh))
        if nrpcs_neigh == 0:
            print('Using .X to calculate nearest neighbors instead of PCs.')
            logging.info(
                'Neighborhood analysis performed with .X instead of PCs.')
    # umap
    sc_umap(adata, random_state=random_state)
    print('UMAP coordinates calculated.')

    logging.info('Neighborhood analysis completed, and UMAP generated.')
    logging.info('\t Time for PCA, nearest neighbor calculation and UMAP generation: ' +
                 str(round(time()-start, 3))+'s')

    # export metadata
    start = time()
    export_metadata(adata, basepath=results_folder,
                    n_pcs=3, umap=True, tsne=False)
    logging.info(
        'Metadata containing 3 PCAs and UMAP coordinates exported successfully to file.')
    logging.info('Time for export: '+str(round(time()-start, 3))+'s')

    return(adata)


def clustering(adata, results_folder, myres=1, method='leiden'):
    """ Perform adata clustering and write the corresponding results

    parameters
    ----------
    adata: `ÀnnData`
        AnnData object that is to be exported
    results_folder: `str`
        path to the results folder 
    myres: int
        resolution for the algorithm
    method: `str`
        clustering algorithm. Implemented: louvain/leiden

    returns
    -------
    None
        writes to file

    """
    if(not method in ['leiden', 'louvain']):
        raise ValueError("method argument should be leiden or louvain")
    random_state = 0
    start = time()
    if method == 'louvain':
        sc_louvain(adata, resolution=myres, random_state=random_state)
    if method == 'leiden':
        sc_leiden(adata, resolution=myres, random_state=random_state)
    print(method + ' clustering performed with a resolution of '+str(myres))
    clusNum = len(set(adata.obs[method]))

    sc_pl_umap(adata, color=[method],
               legend_loc='on data', save='.' + method + '.png')

    logging.info(method + 'clustering done. Found ' +
                 str(clusNum) + ' clusters.')
    logging.info('\tTime for ' + method + ' clustering: ' +
                 str(round(time()-start, 3))+'s')

    # detect marker genes
    start = time()
    rank_genes_groups(adata, method, method='wilcoxon',
                      use_raw=True, n_genes=adata.raw.X.shape[1])
    print('rank genes per cluster calculated using method wilcoxon.')

    logging.info(
        'Marker gene detection performed on a per-cluster basis using the method wilcoxon.')
    logging.info('\tTime for marker gene detection: ' +
                 str(round(time()-start, 3))+'s')

    # export  clustering to file
    start = time()
    export_clustering(adata, basepath=results_folder, method=method)
    export_rank(adata, basepath=results_folder,
                type='wilcox', labeling_name=method)
    logging.info('Cluster level analysis and marker genes exported to file.')
    logging.info('\tTime for export of cluster level analysis: ' +
                 str(round(time()-start, 3))+'s')

    return(adata)


def additional_labeling(adata, labeling_to_use, labeling_name, labeling_description, labeling_author, results_folder):
    """ Standard Workflow function to export an additional labeling besides louvain to FAIR format.

    This function calculated marker genes per label (using rank_genes_groups and the method 'wilcoxon'), exports the labeling,
    generates a labeling_info file, and exports the rank file.

    parameters
    ----------
    adata: `AnnData`
      AnnData object from which the labeling is to be exported
    labeling_to_use: `str`
      string identifying the column in adata.obs containing the labeling that is to be exported (also used
      to calculate the ranked_genes)
    labeling_name: `str`
      string identifiying under which name the labeling should be exported
    labeling_description: `str`
      string defining the description which should be saved in the labeling_info file for the exported labeling
    labeling_author: `str`
      string defining the author of the labeling which should be saved in the labeling_info file for the exported labeling
    results_folder: `str`
      string indicating the basepath to the results folder which is automatically generated when using the standard workflow (pass results_folder)

    returns
    -------
    None
      writes out several files to folder results_folder/labelings/<labeling_name>
    """
    start = time()

    # calculate marker genes for labeling
    rank_genes_groups(adata, labeling_to_use, method='wilcoxon',
                      use_raw=True, n_genes=adata.raw.X.shape[1])
    print('rank genes per label calculated using method wilcoxon.')

    logging.info(
        'Marker gene detection performed on a per-label basis using the method wilcoxon.')
    logging.info('\tTime for marker gene detection: ' +
                 str(round(time()-start, 3))+'s')

    outpath_ = os.path.join(results_folder, 'labelings', labeling_name)
    # export labeling
    start = time()
    labeling(adata, column=labeling_to_use, outpath=outpath_)
    # generate labelinfo.tsv file
    labeling_info(outpath=outpath_,
                  description=labeling_description,
                  public=True,
                  default=False,
                  expert=False,
                  reference=True,
                  method=labeling_author,
                  annotated_version_of=' -')
    export_rank(adata, basepath=results_folder,
                type='wilcox', labeling_name=labeling_name)

    logging.info('Label level analysis and marker genes exported to file.')
    logging.info('\tTime for export of cluster level analysis: ' +
                 str(round(time()-start, 3))+'s')
    return(adata)


def celltype_labeling(adata, labeling_author, results_folder, labeling_to_use='celltype', labeling_name='celltype', labeling_description='manual celltype annotation', cluster_method='louvain'):
    """ Standard Workflow function to export an additional labeling besides louvain to FAIR format.

    This function calculated marker genes per label (using rank_genes_groups and the method 'wilcoxon'), exports the labeling,
    generates a labeling_info file, and exports the rank file.

    parameters
    ----------
    adata: `AnnData`
      AnnData object from which the labeling is to be exported
    labeling_to_use: `str` | default = 'celltype'
      string identifying the column in adata.obs containing the labeling that is to be exported (also used
      to calculate the ranked_genes)
    labeling_name: `str` | default = 'celltype'
      string identifiying under which name the labeling should be exported
    labeling_description: `str` | default = 'manual celltype annotation'
      string defining the description which should be saved in the labeling_info file for the exported labeling
    labeling_author: `str`
      string defining the author of the labeling which should be saved in the labeling_info file for the exported labeling
    results_folder: `str`
      string indicating the basepath to the results folder which is automatically generated when using the standard workflow (pass results_folder)

    returns
    -------
    None
      writes out several files to folder results_folder/labelings/<labeling_name>
    """
    start = time()
    # calculate marker genes for labeling
    rank_genes_groups(adata, labeling_to_use, method='wilcoxon',
                      use_raw=True, n_genes=adata.raw.X.shape[1])
    print('rank genes per label calculated using method wilcoxon.')
    logging.info(
        'Marker gene detection performed on a per-label basis using the method wilcoxon.')
    logging.info('\tTime for marker gene detection: ' +
                 str(round(time()-start, 3))+'s')
    # export labeling
    outpath = os.path.join(results_folder, 'labelings', labeling_name)
    start = time()
    labeling(adata, column=labeling_to_use, outpath=outpath)
    # generate labelinfo.tsv file
    labeling_info(outpath=outpath,
                  description=labeling_description,
                  public=False,
                  default=False,
                  expert=True,
                  reference=False,
                  method=labeling_author,
                  annotated_version_of=cluster_method)

    export_rank(adata, basepath=results_folder,
                type='wilcox', labeling_name=labeling_name)

    logging.info('Label level analysis and marker genes exported to file.')
    logging.info('\tTime for export of cluster level analysis: ' +
                 str(round(time()-start, 3))+'s')

    return(None)
