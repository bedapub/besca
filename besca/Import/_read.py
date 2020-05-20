import pandas as pd
from scanpy.api import read
import os
from .._helper import convert_ensembl_to_symbol

import sys

def assert_filepath(filepath):
    """Asserts that the filepath contains the files required by read_mtx.

    Parameters
    ----------
    filepath: `str`
        filepath as string to the directory that is expected to contain at least
        following files: `matrix.mtx`, `genes.tsv`, and `barcodes.tsv`.

    Returns
    -------
    None if the filepath is valid. If not, an Exception will be raised
    """
    req_files = ['matrix.mtx', 'genes.tsv', 'barcodes.tsv']
    req_files_exist = [os.path.isfile(os.path.join(filepath, x)) for x in
            req_files]
    valid = all(req_files_exist)
    if valid:
        return(None)
    else:
        for ind, exist  in enumerate(req_files_exist):
            if not exist:
                raise FileNotFoundError('{} is not found in '
                        'the given path `{}`'.format(req_files[ind], filepath))

def read_mtx(
        filepath,
        annotation = True,
        use_genes = 'SYMBOL',
        species = 'human',
        citeseq = None):

    """Read matrix.mtx, genes.tsv, barcodes.tsv to AnnData object.
    By specifiying an input folder this function reads the contained matrix.mtx,
    genes.tsv and barcodes.tsv files to an AnnData object. In case annotation = True
    it also adds the annotation contained in metadata.tsv to the object.
    Parameters
    ----------
    filepath: `str`
        filepath as string to the directory containg the matrix.mtx, genes.tsv,
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

    assert_filepath(filepath)

    print('reading matrix.mtx')
    adata = read(os.path.join(filepath, 'matrix.mtx'), cache= True).T  #transpose the data

    print('adding cell barcodes')
    adata.obs_names = pd.read_csv(os.path.join(filepath, 'barcodes.tsv'), header=None)[0]

    print('adding genes')
    var_anno = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')
    symbols = var_anno[1]
    ensembl_id = var_anno[0]

    adata.var['ENSEMBL'] = ensembl_id.tolist()
    adata.var.index.names = ['index']

    if use_genes == 'SYMBOL':
        adata.var_names = symbols

        #make unique
        print('making var_names unique')
        adata.var_names_make_unique()

        if symbols.tolist() != ensembl_id.tolist():
            print('adding ENSEMBL gene ids to adata.var')
            adata.var['SYMBOL'] = symbols.tolist()
    elif use_genes == 'ENSEMBL':
        adata.var_names = ensembl_id

        if symbols[0] != ensembl_id[0]:
            #add symbols to object
            print('adding symbols to adata.var')
            adata.var['SYMBOL'] = symbols.tolist()
        else:
            #lookup the corresponding Symbols
            adata.var['SYMBOL'] = convert_ensembl_to_symbol(ensembl_id.tolist(), species = species)
    else:
        sys.exit('supplied unknown use_genes parameter')

    if citeseq is not None:
        features = var_anno[2]
        adata.var['feature_type'] = features.tolist()

        #only extract the information we want
        if citeseq == 'gex_only':
                adata = adata[:, adata.var.feature_type == 'Gene Expression'].copy()
        if citeseq == 'citeseq_only':
            adata = adata[:, adata.var.feature_type == 'Antibody Capture'].copy()
            
    if annotation == True:
        print('adding annotation')
        adata.obs = pd.read_csv(os.path.join(filepath, 'metadata.tsv'), sep='\\t',engine='python')
        if adata.obs.get('CELL') is not None:
            adata.obs.index = adata.obs.get('CELL').tolist()

    return(adata)
