import pandas as pd
from scanpy.api import read
import os
from .._helper import convert_ensembl_to_symbol

import sys

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
    if use_genes == 'SYMBOL':
        if os.path.isfile(os.path.join(filepath, 'matrix.mtx')):
            print('reading matrix.mtx')
            
            adata = read(os.path.join(filepath, 'matrix.mtx'), cache= True).T  #transpose the data
           
            print('adding genes')
            adata.var_names = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[1] # Use gene symbols
            print('adding cell barcodes')
            adata.obs_names = pd.read_csv(os.path.join(filepath, 'barcodes.tsv'), header=None)[0]

            #make unique
            print('making var_names unique')
            adata.var_names_make_unique()
            
            #add ENSEMBL ids if present in genes.tsv file
            symbols = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[1] #get gene symbols
            ENSEMBL_id = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[0] #get ENSEMBL Ids
            if symbols.tolist() != ENSEMBL_id.tolist():
                print('adding ENSEMBL gene ids to adata.var')
                adata.var['ENSEMBL'] = ENSEMBL_id.tolist()
                adata.var.index.names = ['index']
                adata.var['SYMBOL'] = symbols.tolist()
            
            if citeseq is not None:
                features = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[2]
                adata.var['feature_type'] = features.tolist()

            if annotation == True:
                print('adding annotation')
                adata.obs = pd.read_csv(os.path.join(filepath, 'metadata.tsv'), sep='\\t',engine='python')
                if adata.obs.get('CELL') is not None:
                    adata.obs.index = adata.obs.get('CELL').tolist()
            
            if citeseq == 'gex_only':
                adata = adata[:, adata.var.feature_type == 'Gene Expression'].copy()
            if citeseq == 'citeseq_only':
                adata = adata[:, adata.var.feature_type == 'Antibody Capture'].copy()
            
            return(adata)
        
        if not os.path.exists(os.path.join(filepath, 'matrix.mtx')):
            raise ValueError('Reading file failed, file does not exist.')

    elif use_genes == 'ENSEMBL':
        if os.path.isfile(os.path.join(filepath, 'matrix.mtx')):
            print('reading matrix.mtx')
            
            adata = read(os.path.join(filepath, 'matrix.mtx'), cache= True).T  #transpose the data
            print('adding genes')
            adata.var_names = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[0] # Use gene symbols
            print('adding cell barcodes')
            adata.obs_names = pd.read_csv(os.path.join(filepath, 'barcodes.tsv'), header=None)[0]

            #add symbols to object
            symbols = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[1] #get gene symbols
            ENSEMBL_id = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[0] #get ENSEMBL Ids
            if symbols[0] != ENSEMBL_id[0]:
                print('adding symbols to adata.var')
                adata.var['SYMBOL'] = symbols.tolist()
                adata.var.index.names = ['index']
                adata.var['ENSEMBL'] = ENSEMBL_id.tolist()

            #lookup the corresponding Symbols
            else:
                adata.var['ENSEMBL'] = ENSEMBL_id.tolist()
                adata.var.index.names = ['index']
                adata.var['SYMBOL'] = convert_ensembl_to_symbol(ENSEMBL_id.tolist(), species = species)

            if citeseq is not None:
                features = pd.read_csv(os.path.join(filepath, 'genes.tsv'), header=None, sep='\\t',engine='python')[2]
                adata.var['feature_type'] = features.tolist()

            if annotation == True:
                print('adding annotation')
                adata.obs = pd.read_csv(os.path.join(filepath, 'metadata.tsv'), sep='\\t',engine='python')
                if adata.obs.get('CELL') is not None:
                    adata.obs.index = adata.obs.get('CELL').tolist()

            if citeseq == 'gex_only':
                adata = adata[:, adata.var.feature_type == 'Gene Expression'].copy()
            if citeseq == 'citeseq_only':
                adata = adata[:, adata.var.feature_type == 'Antibody Capture'].copy()
            
            return(adata)
        
        if not os.path.exists(os.path.join(filepath, 'matrix.mtx')):
            raise ValueError('Reading file failed, file does not exist.')
    else:
        sys.exit('supplied unknown use_genes parameter')
