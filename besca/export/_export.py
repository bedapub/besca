import os
from pandas import DataFrame, read_csv, concat
from numpy import round, ndarray, where, arange, expm1, zeros, ix_
from scipy import io, sparse
import sys
from io import BytesIO
import pkg_resources

#Functions to export AnnData objects to FAIR dataformat

def X_to_mtx(adata,
             outpath = None,
             write_metadata = False, 
             geneannotation = 'SYMBOL',
             additional_geneannotation = 'ENSEMBL', 
             C_path = None):
    """export adata object to mtx format (matrix.mtx, genes.tsv, barcodes.tsv)

    exports the counts contained in adata.X to a matrix.mtx file (in sparse format),
    genes taken from adata.var_names (and if applicable adata.var) to genes.tsv and the 
    cellbarcodes from adata.obs_names to barcodes.tsv. If annotation = True, then the entire 
    pd.Dataframe contained in adata.obs is in addition exported to metadata.tsv.
    Through the parameter geneannotation you can specify which type of geneidentifer is saved in 
    adata.var_names. You can pass an additional string to the parameter additional_geneannotation
    which specifies under which column name in adata.var an additional geneannotation can be 
    found. Currently this function is only capable of dealing with geneannotation of the type 
    ENSEMBL or SYMBOL. This feature is intended to conserve the correct order of ENSEMBL IDs and 
    SYMBOLS in the genes.tsv file.
    If the outpath directory does not exist, this function automatically generates it.
    
    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    write_metadata: `bool` | default = False
        boolian indicator if the annotation contained in adata.obs should be exported as well 
    geneannotation: `'ENSEMBL'` or `'SYMBOL'`
        string indicator of the type of gene annotation saved in adata.var_names
    additional_geneannotation: `str` | default = None
        string identifying the coloumn name in which either the SYMBOL or the ENSEMBL geneids
        are contained as additional gene annotation in adata.var
    C_path: `string` | default = 'pkg_resources.resource_filename('besca', 'export/reformat')'
        File path to the location of the C-programm that is called to reformat the matrix.mtx file
    
    returns
    -------
    None
        writes out files to the specified output directory
    
    """
    if outpath is None:
        outpath = os.getcwd()
    if C_path is None:
        C_path = pkg_resources.resource_filename('besca', 'export/reformat')
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    ### write out matrix.mtx as float with 3 significant digits
    print('writing out matrix.mtx ...')

    if type(adata.X) == ndarray:
        E = sparse.csr_matrix(adata.X).T
    else:
        E = adata.X.tocsr().T #transpose back into the same format as we imported it

    
    #check if executable is actually executable on the system
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if is_exe(C_path):
        #create a temporary file to pipe the output from io.mmwrite to
        f_in = BytesIO()
        
        #write out temporary matrix to created pipe
        io.mmwrite(target=f_in, a=E, precision = 2)

        from subprocess import run
        f_out = open(os.path.join(outpath, 'matrix.mtx'), "w")
        run([C_path, "-"], input=f_in.getvalue(), stdout=f_out)

        print('adata.X successfully written to matrix.mtx')

    else:
        io.mmwrite(target=os.path.join(outpath, 'matrix.mtx'), a=E, precision = 2)

        print('adata.X successfully written to matrix.mtx. \n Warning: could not use reformat script.')
        
    ### export genes
    
    #get genes to export
    if geneannotation == 'SYMBOL':
        genes_SYMBOL = adata.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_ENSEMBL = adata.var.get(additional_geneannotation)
        else:
            print('No ensembel gene ids provided, will fill the respective columns in genes.tsv with NA')
            genes_ENSEMBL = ['NA']*len(genes_SYMBOL)

    elif geneannotation == 'ENSEMBL':
        genes_ENSEMBL = adata.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_SYMBOL = adata.var.get(additional_geneannotation)
        else:
            print('No SYMBOLSprovided, will fill the respective columns in genes.tsv with NA')
            genes_SYMBOL = ['NA']*len(genes_ENSEMBL)

    else:
        sys.exit('need to provide either \'ENSEMBL\' or \'SYMBOL\' gene annotation.')

    feature = None
    #add catch to write out type of annotation 
    if 'feature_type' in adata.var.columns:
        print('feature annotation is present and will be written out')
        feature = True
        gene_feature = adata.var.get('feature_type')

    #write the genes out in the correct format (first ENSEMBL THEN SYMBOL)
    with open (os.path.join(outpath, 'genes.tsv'), "w") as fp:
        if feature is not None:
            for ENSEMBL, symbol, feature in zip(genes_ENSEMBL, genes_SYMBOL, gene_feature):
                fp.write(ENSEMBL+"\t"+ symbol+"\t"+ feature +"\n")
        else:
            for ENSEMBL, symbol in zip(genes_ENSEMBL, genes_SYMBOL):
                fp.write(ENSEMBL+"\t"+ symbol+"\n")
        fp.close()
        print('genes successfully written out to genes.tsv')

    ### write out the cellbarcodes
    cellbarcodes = adata.obs_names.tolist()
    with open(os.path.join(outpath, 'barcodes.tsv'), "w") as fp:
        for barcode in cellbarcodes:
            fp.write(barcode+"\n")
        fp.close()
        print('cellbarcodes successfully written out to barcodes.tsv')

    ### export annotation

    if write_metadata == True:
        annotation = adata.obs
        annotation.to_csv(os.path.join(outpath, 'metadata.tsv'), sep = '\t', header = True, index = None)
        print('annotation successfully written out to metadata.tsv')

    return(None)
    sys.exit(0)

def raw_to_mtx(adata,
               outpath = os.getcwd(),
               write_metadata = False,
               geneannotation = 'ENSEMBL',
               additional_geneannotation = None,
               C_path = None ):

    """ export adata.raw to .mtx (matrix.mtx, genes.tsv, barcodes, tsv)

    exports the counts contained in adata.raw.X to a matrix.mtx file (in sparse format),
    genes taken from adata.raw.var_names to genes.tsv and the 
    cellbarcodes from adata.obs_names to barcodes.tsv. If annotation = True, then the entire 
    pd.Dataframe contained in adata.obs is in addition exported to metadata.tsv.

    Through the parameter genetannotation you can specify which type of geneidentifer is saved in 
    adata.raw.var_names. You can pass an additional string to the parameter additional_geneannotation
    which specifies under which column name in adata.raw.var an additional geneannotation can be 
    found. Currently this function is only capable of dealing with geneannotation of the type 
    ENSEMBL or SYMBOL. This feature is intended to conserve the correct order of ENSEMBL IDs and 
    SYMBOL symbols in the genes.tsv file.

    If the outpath directory does not exist, this function automatically generates it.

    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    write_metadata: `bool` | default = False
        boolian indicator if the annotation contained in adata.obs should be exported as well 
    geneannotation: `'ENSEMBL` or `'SYMBOL`
        string indicator of the type of gene annotation saved in adata.var_names
    additional_geneannotation: `str` | default = None
        string identifying the coloumn name in which either the SYMBOL symbols or the ENSEMBL geneids
        are contained as additional gene annotation in adata.var
    C_path: `string` | default = pkg_resources.resource_filename('besca', 'export/reformat')
        File path to the location of the C-programm that is called to reformat the matrix.mtx file

    returns
    -------
    None
        writes out files to the specified output directory

    """
    if C_path is None:
        C_path= pkg_resources.resource_filename('besca', 'export/reformat')
    ### write out matrix.mtx as float with 3 significant digits
    print('writing out matrix.mtx ...')
    if adata.raw.X is None:
        sys.exit(1, 'adata does not have .raw')
    elif type(adata.raw.X) == ndarray:
        E = sparse.csr_matrix(adata.raw.X).T
    else:
        E = adata.raw.X.tocsr().T #transpose back into the same format as we imported it
    
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        
    #create a temporary file to pipe the output from io.mmwrite to
    f_in = BytesIO()
    #write out temporary matrix to created pipe
    io.mmwrite(target=f_in, a=E, precision = 2)
    
    from subprocess import run
    f_out = open(os.path.join(outpath, 'matrix.mtx'), "w")
    run([C_path, "-"], input=f_in.getvalue(), stdout=f_out)
    
    print('adata.raw.X successfully written to matrix.mtx')

    ### export genes
    
    #get genes to export
    if geneannotation== 'SYMBOL':
        genes_SYMBOL = adata.raw.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_ENSEMBL = adata.raw.var.get(additional_geneannotation)
        else:
            print('No ensembel gene ids provided, will fill the respective columns in genes.tsv with NA')
            genes_ENSEMBL = ['NA']*len(genes_SYMBOL)

    elif geneannotation == 'ENSEMBL':
        genes_ENSEMBL = adata.raw.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_SYMBOL = adata.raw.var.get(additional_geneannotation)
        else:
            print('No SYMBOL symbols provided, will fill the respective columns in genes.tsv with NA')
            genes_SYMBOL = ['NA']*len(genes_ENSEMBL)
    else:
        sys.exit('need to provide either \'ENSEMBL\' or \'SYMBOL\' gene annotation.')

    feature = None
    #add catch to write out type of annotation 
    if 'feature_type' in adata.var.columns:
        print('feature annotation is present and will be written out')
        feature = True
        gene_feature = adata.var.get('feature_type')

    #write the genes out in the correct format (first ENSEMBL THEN SYMBOL)
    with open (os.path.join(outpath, 'genes.tsv'), "w") as fp:
        if feature is not None:
            for ENSEMBL, symbol, feature in zip(genes_ENSEMBL, genes_SYMBOL, gene_feature):
                fp.write(ENSEMBL+"\t"+ symbol+"\t"+ feature +"\n")
        else:
            for ENSEMBL, symbol in zip(genes_ENSEMBL, genes_SYMBOL):
                fp.write(ENSEMBL+"\t"+ symbol+"\n")
        fp.close()
        print('genes successfully written out to genes.tsv from adata.raw')

    ### write out the cellbarcodes
    cellbarcodes = adata.obs_names.tolist()
    with open(os.path.join(outpath, 'barcodes.tsv'), "w") as fp:
        for barcode in cellbarcodes:
            fp.write(barcode+"\n")
        fp.close()
        print('cellbarcodes successfully written out to barcodes.tsv')

    ### export annotation

    if write_metadata == True:
        annotation = adata.obs
        annotation.to_csv(os.path.join(outpath, 'metadata.tsv'), sep = '\t', header = True, index = None)
        print('annotation successfully written out to metadata.tsv')

    return(None)
    sys.exit(0)

def X_to_mtx_python(adata,
                    outpath = os.getcwd(),
                    write_metadata = False, 
                    geneannotation = 'SYMBOL',
                    additional_geneannotation = 'ENSEMBL'):
    """export adata object to mtx format (matrix.mtx, genes.tsv, barcodes.tsv) (ONLY USING PYTHON)

    exports the counts contained in adata.X to a matrix.mtx file (in sparse format),
    genes taken from adata.var_names (and if applicable adata.var) to genes.tsv and the 
    cellbarcodes from adata.obs_names to barcodes.tsv. If annotation = True, then the entire 
    pd.Dataframe contained in adata.obs is in addition exported to metadata.tsv.

    Through the parameter geneannotation you can specify which type of geneidentifer is saved in 
    adata.var_names. You can pass an additional string to the parameter additional_geneannotation
    which specifies under which column name in adata.var an additional geneannotation can be 
    found. Currently this function is only capable of dealing with geneannotation of the type 
    ENSEMBL or SYMBOL. This feature is intended to conserve the correct order of ENSEMBL IDs and 
    SYMBOLS in the genes.tsv file.

    If the outpath directory does not exist, this function automatically generates it.

    This is a slower version of this function that does not rely on external C-Programms. The default
    version of this function that is used in most cases 'X_to_mtx' does rely on this program.

    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    write_metadata: `bool` | default = False
        boolian indicator if the annotation contained in adata.obs should be exported as well 
    geneannotation: `'ENSEMBL'` or `'SYMBOL'`
        string indicator of the type of gene annotation saved in adata.var_names
    additional_geneannotation: `str` | default = None
        string identifying the coloumn name in which either the SYMBOL or the ENSEMBL geneids
        are contained as additional gene annotation in adata.var

    returns
    -------
    None
        writes out files to the specified output directory

    """

    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    ### write out matrix.mtx as float with 3 significant digits
    print('writing out matrix.mtx ...')

    if type(adata.X) == ndarray:
        E = sparse.csr_matrix(adata.X).T
    else:
        E = adata.X.tocsr().T #transpose back into the same format as we imported it
    #write out temporary matrix
    io.mmwrite(os.path.join(outpath, 'tmp.mtx'), E, precision = 2)
    
    #transform the temporary matrix into the correct version that shows floats with 3 significant digits
    tmp = read_csv(os.path.join(outpath, 'tmp.mtx'), skiprows=3, sep=' ', index_col=None, header=None)
    with open (os.path.join(outpath, 'tmp.mtx'),"r") as fp_in:
        with open (os.path.join(outpath, 'matrix.mtx'),"w") as fp_out:
    # first, copy the first three lines
            for _ in range(3):
                line = fp_in.readline()
                fp_out.write(line)
    with open (os.path.join(outpath, 'matrix.mtx'),"a") as fp_out:
        tmp.to_csv(fp_out, sep=' ', header=None, float_format="%.3f", index=False)

    
    #remove temporary file
    os.remove(os.path.join(outpath, 'tmp.mtx'))
    print('adata.X successfully written to matrix.mtx')

    ### export genes
    
    #get genes to export
    if geneannotation == 'SYMBOL':
        genes_SYMBOL = adata.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_ENSEMBL = adata.var.get(additional_geneannotation)
        else:
            print('No ensembel gene ids provided, will fill the respective columns in genes.tsv with NA')
            genes_ENSEMBL = ['NA']*len(genes_SYMBOL)

    elif geneannotation == 'ENSEMBL':
        genes_ENSEMBL = adata.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_SYMBOL = adata.var.get(additional_geneannotation)
        else:
            print('No SYMBOLSprovided, will fill the respective columns in genes.tsv with NA')
            genes_SYMBOL = ['NA']*len(genes_ENSEMBL)

    else:
        sys.exit('need to provide either \'ENSEMBL\' or \'SYMBOL\' gene annotation.')

    feature = None
    #add catch to write out type of annotation 
    if 'feature_type' in adata.var.columns:
        print('feature annotation is present and will be written out')
        feature = True
        gene_feature = adata.var.get('feature_type')

    #write the genes out in the correct format (first ENSEMBL THEN SYMBOL)
    with open (os.path.join(outpath, 'genes.tsv'), "w") as fp:
        if feature is not None:
            for ENSEMBL, symbol, feature in zip(genes_ENSEMBL, genes_SYMBOL, gene_feature):
                fp.write(ENSEMBL+"\t"+ symbol+"\t"+ feature +"\n")
        else:
            for ENSEMBL, symbol in zip(genes_ENSEMBL, genes_SYMBOL):
                fp.write(ENSEMBL+"\t"+ symbol+"\n")
        fp.close()
        print('genes successfully written out to genes.tsv')


    ### write out the cellbarcodes
    cellbarcodes = adata.obs_names.tolist()
    with open(os.path.join(outpath, 'barcodes.tsv'), "w") as fp:
        for barcode in cellbarcodes:
            fp.write(barcode+"\n")
        fp.close()
        print('cellbarcodes successfully written out to barcodes.tsv')

    ### export annotation

    if write_metadata == True:
        annotation = adata.obs
        annotation.to_csv(os.path.join(outpath, 'metadata.tsv'), sep = '\t', header = True, index = None)
        print('annotation successfully written out to metadata.tsv')

    return(None)
    sys.exit(0)

def raw_to_mtx_python(adata,
                      outpath = os.getcwd(),
                      write_metadata = False,
                      geneannotation = 'ENSEMBL',
                      additional_geneannotation = None):
    """ export adata.raw to .mtx (matrix.mtx, genes.tsv, barcodes, tsv) (ONLY USING PYTHON)

    exports the counts contained in adata.raw.X to a matrix.mtx file (in sparse format),
    genes taken from adata.raw.var_names to genes.tsv and the 
    cellbarcodes from adata.obs_names to barcodes.tsv. If annotation = True, then the entire 
    pd.Dataframe contained in adata.obs is in addition exported to metadata.tsv.

    Through the parameter genetannotation you can specify which type of geneidentifer is saved in 
    adata.raw.var_names. You can pass an additional string to the parameter additional_geneannotation
    which specifies under which column name in adata.raw.var an additional geneannotation can be 
    found. Currently this function is only capable of dealing with geneannotation of the type 
    ENSEMBL or SYMBOL. This feature is intended to conserve the correct order of ENSEMBL IDs and 
    SYMBOL symbols in the genes.tsv file.

    If the outpath directory does not exist, this function automatically generates it.

    This is a slower version of this function that does not rely on external C-Programms. The default
    verison of this function that is used in most cases 'X_to_mtx' does rely on this program.

    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    write_metadata: `bool` | default = False
        boolian indicator if the annotation contained in adata.obs should be exported as well 
    geneannotation: `'ENSEMBL` or `'SYMBOL`
        string indicator of the type of gene annotation saved in adata.var_names
    additional_geneannotation: `str` | default = None
        string identifying the coloumn name in which either the SYMBOL symbols or the ENSEMBL geneids
        are contained as additional gene annotation in adata.var

    returns
    -------
    None
        writes out files to the specified output directory

    """

    ### write out matrix.mtx as float with 3 significant digits
    print('writing out matrix.mtx ...')
    if adata.raw.X is None:
        sys.exit(1, 'adata does not have .raw')
    elif type(adata.raw.X) == ndarray:
        E = sparse.csr_matrix(adata.raw.X).T
    else:
        E = adata.raw.X.tocsr().T #transpose back into the same format as we imported it
    
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    #write out temporary matrix
    io.mmwrite(os.path.join(outpath, 'tmp.mtx'), E, precision = 2)
    
    #transform the temporary matrix into the correct version that shows floats with 3 significant digits
    with open(os.path.join(outpath, 'tmp.mtx'), "r") as fp_in:
        with open(os.path.join(outpath, 'matrix.mtx'), "w") as fp_out:
            #copy the first three lines
            for i in range(3):
                line = fp_in.readline()
                fp_out.write(line)
            #now read the rest of the lines
            while True:
                line = fp_in.readline()
                if not line:
                    break
                n1, n2, val = line.split(' ')
                fp_out.write('{} {} {:.3f}'.format(n1, n2, float(val)))
                fp_out.write('\n')
        fp_out.close()
    fp_in.close()
    
    #remove temporary file
    os.remove(os.path.join(outpath, 'tmp.mtx'))
    print('adata.raw.X successfully written to matrix.mtx')

    ### export genes
    
    #get genes to export
    if geneannotation== 'SYMBOL':
        genes_SYMBOL = adata.raw.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_ENSEMBL = adata.raw.var.get(additional_geneannotation)
        else:
            print('No ensembel gene ids provided, will fill the respective columns in genes.tsv with NA')
            genes_ENSEMBL = ['NA']*len(genes_SYMBOL)

    elif geneannotation == 'ENSEMBL':
        genes_ENSEMBL = adata.raw.var_names.tolist()
        #get additional annotation saved in adata.var
        if additional_geneannotation is not None:
            genes_SYMBOL = adata.raw.var.get(additional_geneannotation)
        else:
            print('No SYMBOL symbols provided, will fill the respective columns in genes.tsv with NA')
            genes_SYMBOL = ['NA']*len(genes_ENSEMBL)
    else:
        sys.exit('need to provide either \'ENSEMBL\' or \'SYMBOL\' gene annotation.')

    feature = None
    #add catch to write out type of annotation 
    if 'feature_type' in adata.var.columns:
        print('feature annotation is present and will be written out')
        feature = True
        gene_feature = adata.var.get('feature_type')

    #write the genes out in the correct format (first ENSEMBL THEN SYMBOL)
    with open (os.path.join(outpath, 'genes.tsv'), "w") as fp:
        if feature is not None:
            for ENSEMBL, symbol, feature in zip(genes_ENSEMBL, genes_SYMBOL, gene_feature):
                fp.write(ENSEMBL+"\t"+ symbol+"\t"+ feature +"\n")
        else:
            for ENSEMBL, symbol in zip(genes_ENSEMBL, genes_SYMBOL):
                fp.write(ENSEMBL+"\t"+ symbol+"\n")
        fp.close()
        print('genes successfully written out to genes.tsv')

    ### write out the cellbarcodes
    cellbarcodes = adata.obs_names.tolist()
    with open(os.path.join(outpath, 'barcodes.tsv'), "w") as fp:
        for barcode in cellbarcodes:
            fp.write(barcode+"\n")
        fp.close()
        print('cellbarcodes successfully written out to barcodes.tsv')

    ### export annotation

    if write_metadata == True:
        annotation = adata.obs
        annotation.to_csv(os.path.join(outpath, 'metadata.tsv'), sep = '\t', header = True, index = None)
        print('annotation successfully written out to metadata.tsv')

    return(None)
    sys.exit(0)


def clustering(adata,  
            outpath = None,
            export_average = True,
            export_fractpos = True,
            method = 'leiden'):               
    """ export mapping of cells to clusters to .tsv file

    This function exports the labels saved in adata.obs.method and the corresponding cell barcode to the file cell2labels.tsv.

    parameters
    ----------
    adata:
    outpath: `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    export_average: `bool` | default = True
        boolian indicator if the average gene expression of each cluster should be exported to file
    export_fractpos: `bool` | default = True
        boolian indicator if the fraction of positive cells (i.e. cells that express the gene) should
        be exported to file
    method: `str1 | default = 'leiden'
        string indicating the clustering method used and where to store the results. Shuold be either louvain or leiden
    returns
    -------
    None
        files are written out.

    """
    if outpath is None:
        outpath = os.getcwd()
    if( not method in ['leiden', 'louvain']):
        raise ValueError("method argument should be leiden or louvain")
    cluster_data = adata.obs.get(method).to_frame(name='LABEL')
    if cluster_data is None:
        sys.exit('need to perform ' + method +' clustering before exporting')
    #perform export calling the general export function
    labeling(adata,
             outpath = outpath, 
             column = method, 
             label = 'LABEL', 
             filename = 'cell2labels.tsv', 
             export_average = export_average, 
             export_fractpos= export_fractpos)
    return(None)

    

def labeling_info(outpath = None,
                  description = 'leiden clustering',
                  public = False,
                  default = True,
                  expert = False,
                  reference = False,
                  method = 'leiden',
                  annotated_version_of = '-',
                  filename = 'labelinfo.tsv'):
    """ write out labeling info for uploading to database

    This functions outputs the file labelinfo.tsv which is needed to annotate a written out 
    labeling in the scseq database.
    
    parameters
    ----------
    outpath: `str` | default = current working directory
        The filepath as a string indicating the location where the file should be written out to.
    description: `str` | default = 'leiden clustering'
        string describing what type of information is saved in the corresponding labeling.
    public: `bool` | default = False
        boolian indicator if the contained labeling information is available in the public domain.
    default_ `bool` | default = True
        boolian indicator if the labeling was created using a standardized process e.g. the leiden 
        clusters outputed by the standard pipeline (this should be false if expert is true)
    expert: `bool` | default = False
        boolian indicator if the labeling was created by an 'expert' i.e. manually done (this should 
        be false if default is true)
    reference: `bool` | default = True
        boolian indicator if this is the labeling (e.g. celltype annotation) that should be used for further analysis
        (there should only be one reference labeling per study)
    method: `str` | default = 'leiden'
        string indicating the type of method that was applied, e.g. if the labeling is of a clustering
        which clustering algorithm was used.
    annotated_version_of: `str` | default = '-'
        string identifying of what othe labeling/dataset this is an annotated version of (so for 
        example if the labeling is celltype annotation of a leiden clustering then this would 
        reference the leiden clsutering that was used to obtain the clusters that were then 
        labeled here)
    filename: `str` | default = 'labelinfo.tsv'
        string indicating the filename that should be used. This is per default set to the correct 
        file name for uploading to the scseq database.

    returns
    -------
    None
        results are written out to a file instead

    """
    if outpath is None:
        outpath = os.getcwd()
    if public:
        Public = 'TRUE'
    else:
        Public = 'FALSE'

    if default:
        Default = 'TRUE'
    else:
        Default = 'FALSE'

    if expert:
        Expert = 'TRUE'
    else:
        Expert = 'FALSE'

    if reference:
        Reference = 'TRUE'
    else:
        Reference = 'FALSE'

    ciFile = os.path.join(outpath, filename)
    with open (ciFile,"w") as fp:
        fp.write("description\tisPublic\tisDefault\tisExpert\tisReference\tmethod\tannotated_version_of\n")
        fp.write(description +'\t'+Public+'\t'+Default+'\t'+Expert+'\t'+Reference+'\t'+method+'\t'+annotated_version_of+"\n")
    fp.close()
    print('labelinfo.tsv successfully written out')

    return(None)
    sys.exit(0)

def labeling(adata,
             outpath = None,
             column = 'leiden',
             label  = 'LABEL',
             filename = 'cell2labels.tsv',
             export_average = True,
             export_fractpos = True,
             use_raw = True):
    """export mapping of cells to specified label to .tsv file

    This is a function with which any type of labeling (i.e. celltype annotation, leiden 
    clustering, etc.) can be written out to a .tsv file. The generated file can then also be easily 
    uploaded to the database since it fullfilles the FAIR document standards. 

    To ensure FAIR compatbility label, and file name should not be changed.

    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    column: `str` | default = 'leiden'
        Name of the column in adata.obs that is to be mapped to cell barcodes and written out to file.
    label: `str` | default = 'LABEL'
        label above the column when it is written out to file
    filename: `str` | default = 'cell2labels.tsv'
        Filename that is written out.

    returns
    -------
    None
        files are written out.

    """
    if outpath is None:
        outpath = os.getcwd()
    data = adata.obs.get(column).to_frame(name=label)
    if data is None:
        sys.exit('please specify column name that is present in adata.obs')
    
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    data.to_csv(os.path.join(outpath, filename), sep = '\t', header=True, index_label='CELL')
    print('mapping of cells to ', column, 'exported successfully to', filename)
    
    if export_average:
        if (column == 'louvain' or column == 'leiden'):
            labeling = sorted(set(adata.obs[column].astype(int)))
        else:
            labeling = sorted(set(adata.obs[column]))

        label_names = []
        for i in range(len(labeling)):
            label_names.append(str(labeling[i]))

        num_labels = len(label_names)

        if use_raw:
            gct = DataFrame(index = adata.raw.var_names, columns = label_names)
            mydat = adata.raw.var.copy()
            if type(adata.raw.X) == ndarray:
                E = sparse.csr_matrix(adata.raw.X).T
            else:
                E = adata.raw.X.tocsr().T

        else:
            gct = DataFrame(index = adata.var_names, columns = label_names)
            mydat = adata.var.copy()
            if type(adata.X) == ndarray:
                E = sparse.csr_matrix(adata.X).T
            else:
                E = adata.X.tocsr().T

        #revert to linear scale
        E = E.expm1()

        for i in range(num_labels):
            cells = where(adata.obs.get(column) == label_names[i])[0]
            gct.iloc[:,i] = E[:,cells].mean(axis= 1) #get mean expression per gene
            labeling_size = len(cells)

        mydat = mydat.loc[:,['SYMBOL', 'ENSEMBL']]
        mydat.rename(columns={'SYMBOL':'Description'}, inplace=True)
        gct = mydat.merge(gct, how='right', left_index=True, right_index=True)
        
        gct.set_index('ENSEMBL', inplace=True)
        gct.index.names = ['NAME']

        #write out average expression
        gctFile_average = os.path.join(outpath, 'average.gct')
        with open (gctFile_average,"w") as fp:
            fp.write("#1.2"+"\n")
            fp.write(str(gct.shape[0])+'\t'+str(gct.shape[1] - 1)+'\n') # "description" already merged in as a column
        fp.close()
        #...and then the matrix
        gct.to_csv(gctFile_average, sep = '\t', index=True, index_label='NAME', header=True, mode = 'a',  float_format='%.3f')
        print('average.gct exported successfully to file')

    if export_fractpos:
        if (column == 'louvain' or column == 'leiden'):
            labeling = sorted(set(adata.obs[column].astype(int)))
        else:
            labeling = sorted(set(adata.obs[column]))

        label_names = []
        for i in range(len(labeling)):
            label_names.append(str(labeling[i]))

        num_labels = len(label_names)

        if use_raw:
            f = DataFrame(index=adata.raw.var_names, columns=label_names)
            mydat = adata.raw.var.copy()
            if type(adata.raw.X) == ndarray:
                E = sparse.csr_matrix(adata.raw.X).T
            else:
                E = adata.raw.X.tocsr().T
        else:
            gct = DataFrame(index = adata.var_names, columns = label_names)
            f = DataFrame(index=adata.var_names, columns=label_names)
            mydat = adata.var.copy()
            if type(adata.X) == ndarray:
                E = sparse.csr_matrix(adata.X).T
            else:
                E = adata.X.tocsr().T
        #revert to linear scale
        E = E.expm1()

        for i in range(num_labels):
            cells = where(adata.obs.get(column) == label_names[i])[0]
            a = E[:,cells].getnnz(axis = 1) #get number of values not 0
            f.iloc[:,i] = a.copy()
        f = f.astype(float)
        for i in range(num_labels):
            cells = where(adata.obs.get(column) == label_names[i])[0]
            labeling_size   = len(cells)
            f[label_names[i]] = f[label_names[i]]/labeling_size

        mydat = mydat.loc[:,['SYMBOL', 'ENSEMBL']]
        mydat.rename(columns={'SYMBOL':'Description'}, inplace=True)

        f = mydat.merge(f, how='right', left_index=True, right_index=True)
        f.set_index('ENSEMBL', inplace=True)
        f.index.names = ['NAME']

        #write out frac_pos.gct
        gctFile_fracpos = os.path.join(outpath, 'fract_pos.gct')
        with open (gctFile_fracpos,"w") as fp:
            fp.write("#1.2"+"\n")
            fp.write(str(f.shape[0])+'\t'+str(f.shape[1] - 1)+'\n') # "description" already merged in as a column
        fp.close()
        #...and then the matrix
        f.to_csv(gctFile_fracpos, sep = '\t', index=True, index_label='NAME', header=True, mode = 'a',  float_format='%.3f')
        print('fract_pos.gct exported successfully to file')

    return(None)
    sys.exit(0)


def analysis_metadata(adata,
                      outpath= None,
                      filename='analysis_metadata.tsv',
                      total_counts=True,
                      n_pcs=3,
                      umap=True,
                      tsne=False):
    """export plotting coordinates to analysis_metadata.tsv

    This function exports the indicated plotting coordinates or calculated PCAs to a .tsv file. This
    can be used to either transfer the data between different analysis platforms but can also be 
    uploaded to the scseq database sicne the file follows the FAIR document format.

    To ensure FAIR compatibility the filename should not be changed.
    
    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    filename: `str` | default = 'analysis_metadata.tsv'
        filename of the file that is to be written out
    total_counts: `bool` | default = True
        boolian indicator if total counts should be written out or not
    n_pcs: `int` | default = 3
        indicates number of PCA components that should be written out. If 0 no PCA components will be written out
    umap: `bool` | default = True
        boolian indicator if UMAP coordinates should be written out.
    tsne: `bool` | default = False
        boolian indicator if tSNE coordinates should be written out.

    returns
    -------
    None
        file is written out.

    """
    if outpath is None:
        outpath = os.getcwd()
    #add cellbarcodes to index
    data = DataFrame(data=None, index=adata.obs_names)

    if total_counts:
        if adata.obs.n_counts is None:
            sys.exit('need to have calculated \'n_counts\' and stored it in adata.obs, consider running \"\"')
        else:
            data['totalCounts'] = adata.obs.n_counts.copy()    

    obsm = adata.obsm.to_df()        
    
    if n_pcs > 0:
        for i in range(1, n_pcs+1):
            if obsm.get('X_pca'+str(i)) is None:
                sys.exit('number of PCA components requested not saved in adata.obsm, please ensure that PCA components have been calculated')
            else: 
                data['PCA.PC'+str(i)] = obsm.get('X_pca'+str(i)).tolist()

    if umap:
        if obsm.get('X_umap1') is None:
            sys.exit('no UMAP coordinates found in adata.obsm, please ensure that UMAP coordinates have been calculated')
        else:
            data['UMAP.c1'] = obsm.get('X_umap1').tolist()
            data['UMAP.c2'] = obsm.get('X_umap2').tolist()

    if tsne:
        if obsm.get('X_tsne1') is None:
            sys.exit('no tSNE coordinates found in adata.obsm, please ensure that tSNE coordinates have been calculated')
        else:
            data['tSNE.c1'] = obsm.get('X_tsne1').tolist()
            data['tSNE.c2'] = obsm.get('X_tsne2').tolist()

    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    data.to_csv(os.path.join(outpath, 'analysis_metadata.tsv'), sep='\t', index_label='CELL', float_format='%.3f')
    print('results successfully written out to \'analysis_metadata.tsv\'')
    return(None)
    sys.exit(0)

def ranked_genes(adata, 
                 type = 'wilcox',
                 outpath = None,
                 geneannotation = 'SYMBOL',
                 additional_geneannotation = 'ENSEMBL'):
    """export marker genes for each cluster to .gct file

    This function exports the results of scanpy.api.tl.rank_genes_groups() on your AnnData object to a .gct
    file. This file can easily be uploaded into the scsqe database since it follows the FAIR data
    formats. 

    A prerequisit for executing this function is that sc.tl.rank_genes_groups() has already been run.
    Through the variables geneannotation and additional_geneannotation you can specify the type of
    gene annotationi that is saved in adata.var_names and any additional geneannotation columns saved
    in adata.vars. 

    parameters
    ----------
    adata:
        AnnData object on which scanpy.api.tl.rank_genes_groups has been executed
    type: `str` | 'wilcox' or 't-test overest var'  or 't-test'
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    geneannotation: `'ENSEMBL'` or `'SYMBOL'`
        type of gene annotation that is located in adata.var_names
    additional_geneannotation:

    returns
    -------
    None
        writes results to file in output directory

    """
    if outpath is None:
        outpath = os.getcwd()
    if adata.uns.get('rank_genes_groups') is None:
        sys.exit('need to rank genes before export, please run: scanpy.api.tl.rank_genes() before proceeding with export')
    else:
        #extract relevant data from adata object
        rank_genes = adata.uns['rank_genes_groups']

    #get group names
    groups = rank_genes['names'].dtype.names

    #get number of groups
    group_number = len(groups)

    #get gene information
    mydat = adata.raw.var.loc[:, ['SYMBOL', 'ENSEMBL']]
    mydat.rename(columns={'SYMBOL':'Description'}, inplace=True)
    
    #initialize dataframe for storing the scores
    scores = DataFrame({'n_' + key[:1]: rank_genes[key][groups [0]] for key in ['names', 'scores']})
    scores.rename(columns={'n_s':groups[0], 'n_n':'NAME'}, inplace=True)
    scores.set_index('NAME', inplace=True)
    scores.index = scores.index.astype(str)

    for i in range(1,group_number):
        n_scores = DataFrame({'n_' + key[:1]: rank_genes[key][groups [i]] for key in ['names', 'scores']})
        n_scores.set_index('n_n', inplace=True)
        n_scores.index = n_scores.index.astype(str)
        scores = scores.merge(n_scores, how='left', left_index=True, right_index=True)
        newname = groups[i]
        scores.rename(columns={'n_s':newname}, inplace=True)
    scores = scores.astype(float)

    #merge in gene annotation
    scores = mydat.merge(scores, how = 'right', left_index=True, right_index=True)

    #make index into ENSEMBL instead of symbol
    scores.set_index('ENSEMBL', inplace=True)
    scores.index.names = ['NAME']

    # get pvalues
    #initialize dataframe for storing the pvalues
    pvalues = DataFrame({'n_' + key[:1]: rank_genes[key][groups [0]] for key in ['names', 'pvals']})
    pvalues.rename(columns={'n_p':groups[0], 'n_n':'NAME'}, inplace=True)
    pvalues.set_index('NAME', inplace=True)
    pvalues.index = pvalues.index.astype(str)

    for i in range(1,group_number):
        n_pvalues = DataFrame({'n_' + key[:1]: rank_genes[key][groups [i]] for key in ['names', 'pvals']})
        n_pvalues.set_index('n_n', inplace=True)
        n_pvalues.index = n_pvalues.index.astype(str)
        pvalues = pvalues.merge(n_pvalues, how='left', left_index=True, right_index=True)
        newname = groups[i]
        pvalues.rename(columns={'n_p':newname}, inplace=True)
    pvalues = pvalues.astype(float)

    #merge in pvalues
    pvalues = mydat.merge(pvalues, how = 'right', left_index=True, right_index=True)

    #make index into ENSEMBL instead of symbol
    pvalues.set_index('ENSEMBL', inplace=True)
    pvalues.index.names = ['NAME']

    logFC = DataFrame({'n_' + key[:1]: rank_genes[key][groups [0]] for key in ['names', 'logfoldchanges']})
    logFC.rename(columns={'n_l':groups[0], 'n_n':'NAME'}, inplace=True)
    logFC.set_index('NAME', inplace=True)
    logFC.index = logFC.index.astype(str)

    for i in range(1,group_number):
        n_logFC = DataFrame({'n_' + key[:1]: rank_genes[key][groups [i]] for key in ['names', 'logfoldchanges']})
        n_logFC.set_index('n_n', inplace=True)
        n_logFC.index = n_logFC.index.astype(str)
        logFC = logFC.merge(n_logFC, how='left', left_index=True, right_index=True)
        newname = groups[i]
        logFC.rename(columns={'n_l':newname}, inplace=True)
    logFC = logFC.astype(float)

    #merge in pvalues
    logFC = mydat.merge(logFC, how = 'right', left_index=True, right_index=True)

        #make index into ENSEMBL instead of symbol
    logFC.set_index('ENSEMBL', inplace=True)
    logFC.index.names = ['NAME']

    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if type == 'wilcox':
        gct_rank_File = os.path.join(outpath, "WilxRank.gct")
        gct_pvalue_File = os.path.join(outpath, "WilxRank.pvalues.gct")
        gct_logFC_File = os.path.join(outpath, "WilxRank.logFC.gct")
    elif type == 't-test overest var':
        gct_rank_File= os.path.join(outpath, "tTestRank.gct")
        gct_pvalue_File = os.path.join(outpath, "tTestRank.pvalues.gct")
        gct_logFC_File = os.path.join(outpath, "tTestRank.logFC.gct")
    elif type == 't-test':
        gct_rank_File= os.path.join(outpath, "OverestVarRank.gct")
        gct_pvalue_File = os.path.join(outpath, "OverestVar.pvalues.gct")
        gct_logFC_File = os.path.join(outpath, "OverestVar.logFC.gct")
    else:
        sys.exit('need to specify type as one of \'wilcox\' or \'t-test overest var\'  or \'t-test\'')
    
    #write out rankfile
    with open (gct_rank_File,"w") as fp:
        fp.write("#1.2"+"\n")
        fp.write(str(scores.shape[0])+'\t'+str(scores.shape[1] - 1)+'\n') # "description" already merged in as a column
    fp.close()
    scores.to_csv(gct_rank_File, sep = '\t', index=True, header=True, mode = 'a',  float_format='%.3f')
    print(gct_rank_File, 'written out')

    #write out pvalues
    with open (gct_pvalue_File,"w") as fp:
        fp.write("#1.2"+"\n")
        fp.write(str(pvalues.shape[0])+'\t'+str(pvalues.shape[1] - 1)+'\n') # "description" already merged in as a column
    fp.close()
    pvalues.to_csv(gct_pvalue_File, sep = '\t', index=True, header=True, mode = 'a',  float_format='%.3e')
    print(gct_pvalue_File, 'written out')

    #write out logFC
    with open (gct_logFC_File,"w") as fp:
        fp.write("#1.2"+"\n")
        fp.write(str(logFC.shape[0])+'\t'+str(logFC.shape[1] - 1)+'\n') # "description" already merged in as a column
    fp.close()
    logFC.to_csv(gct_logFC_File, sep = '\t', index=True, header=True, mode = 'a',  float_format='%.3f')
    print(gct_logFC_File, 'written out')

    return(None)
    sys.exit(0)

def pseudobulk(adata,
             outpath = None,
             column = 'celltype0',
             label  = 'celltype0',
             split_condition  = 'donor',
             todrop =['CELL','input.path','percent_mito','n_counts','n_genes','leiden','celltype0','celltype1','celltype2','celltype3','dblabel'], main_condition='CONDITION'):
    """export pseudobulk profiles of cells to .gct files

    This is a function with which any type of labeling (i.e. celltype annotation, louvain 
    clustering, etc.) can be written out to several .gct files as well as a single metadata file. 

    To ensure FAIR compatbility label, and file name should not be changed.

    parameters
    ----------
    adata:
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.
    column: `str` | default = 'celltype0'
        Name of the column in adata.obs that is to be mapped to cell barcodes and written out to file.
    label: `str` | default = 'celltype0'
        label above the column when it is written out to several files
    split_condition: `str` | default = 'experiment'
        the experimental unit, e.g. sample ID
    todrop: `list` 
        Several column headers to be excluded from metadata
    main_condition: `str` | default = 'CONDITION'
        main condition to be outputed in the metadata file
    returns
    -------
    dfmerge: `pd.DataFrame`
        merged dataframe

    """
    if outpath is None:
        outpath = os.getcwd()
    data = adata.obs.get(column)
    if data is None:
        sys.exit('please specify a column name that is present in adata.obs')
        
    data = adata.obs.get(main_condition)
    if data is None:
        sys.exit('please specify a condition name that is present in adata.obs')
    
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    ### create adata subsets for each column value
    adata.obs[split_condition]=adata.obs[split_condition].astype('str')
    adata.obs[split_condition]=adata.obs[split_condition].astype('category')
    adata.obs[column]=adata.obs[column].astype('category')
    
    bulks={}
    myset=list(set(adata.obs[column]))
    for i in myset:
        ii=i.replace(" ", "_") ## to avoid spaces in cell names
        bulks[ii]=adata[adata.obs[column].isin([i])].copy()
    bulks['all']=adata.copy()
    
    ### go through each adata subset and export pseudobulk
    dfbulks={}
    for x in bulks.keys():   
        # sum expression
        auxdata=bulks[x].copy()
        myexp=list(auxdata.obs[split_condition].cat.categories) ### these are all different levels for experiments
        mysums=zeros((len(auxdata.raw.var.index),len(myexp)))
        for i in range(len(myexp)):
            mysums[:,i]=expm1(auxdata[auxdata.obs[split_condition]==myexp[i]].raw.X).sum(axis=0)
        mysums=DataFrame(mysums)
        mysums.index=adata.raw.var.index
        mysums.columns=[x+'.'+y for y in myexp]
        dfbulks[x]=mysums
    
        mydat = auxdata.raw.var.loc[:,['SYMBOL', 'ENSEMBL']]
        mydat.rename(columns={'SYMBOL':'Description'}, inplace=True)
        gct = mydat.merge(dfbulks[x], how='right', left_index=True, right_index=True)
        gct.set_index('ENSEMBL', inplace=True)
        gct.index.names = ['NAME']
        gct.columns=['Description']+myexp
    
        #write out average expression
        gctFile_pseudo = outpath+ 'Pseudobulk-'+label+'-'+x+'.gct'
        with open (gctFile_pseudo,"w") as fp:
            fp.write("#1.2"+"\n")
            fp.write(str(gct.shape[0])+'\t'+str(gct.shape[1] - 1)+'\n') # "description" already merged in as a column
        fp.close()
        #...and then the matrix
        gct.to_csv(gctFile_pseudo, sep = '\t', index=True, index_label='NAME', header=True, mode = 'a',  float_format='%.3f')
        print('Pseudobulk-'+label+'-'+x+'.gct exported successfully to file')

    #### Output into single .tsv file
    dfmerge=concat(dfbulks,axis=1)
    dfmerge.columns = dfmerge.columns.droplevel()
    dfmerge.to_csv(outpath+ 'Pseudobulk-'+label+'.tsv',sep='\t',index_label=False)

    ### Export one metadata file
    myexp=list(adata.obs[split_condition].cat.categories) 
    colindex=range(0,len(adata.obs.columns)) ### replace if only a subset of metadata should be used
    mysums=[]
    for i in range(len(myexp)):
        mysums.append(list(adata[adata.obs[split_condition]==myexp[i]].obs.iloc[:,colindex].iloc[0,:]))
    mysums=DataFrame(mysums).transpose()
    mysums.index=adata[adata.obs[split_condition]==myexp[i]].obs.iloc[:,colindex].columns
    mysums.columns=myexp
    mysums=mysums.transpose().drop(labels=todrop,axis=1,errors='ignore')
    mysums['ID']=list(mysums.index)
    colorder = ['ID',main_condition] + (mysums.columns.drop(['ID',main_condition]).tolist())
    mysums.loc[:,colorder].to_csv(outpath+ 'Pseudobulk.meta',sep='\t',index=False)

    return(dfmerge)
    sys.exit(0)
    
def generate_gep(adata,
                 filename='gep_basis_vector.csv',
                 column='<last_column>',
                 annot='ENSEMBL',
                 outpath = None):
    """Generate Gene Expression Profile (GEP) from scRNA-seq annotations

    Reads in the AnnData object, taking only the pre-filtered highly variable genes 
    (determined in the BESCA annotation workflow), to index the adata.raw expression matrix. 
    The adata.raw matrix is log1p from the workflow, thus we linearize it before summing the values 
    across each cell type.

    We generate the GEP from adata.raw and not from adata because we need CP10k normalised values, 
    which adata doesnt contain as it will have gone through several normalisation steps further downstream.
    At the same time we are subsetting adata.raw by the highly variable genes present only present in adata.

    For each cell type, its gene expression is calculted by summing up all values for given gene 
    and given cell type. A mean value is not taken as there are many cell with zero expresion for given gene.

    parameters
    ----------
    adata:
    filename 'str' | default = 'gep_basis_vector.csv'
        name of output file
    column: `str` | default = '<last_column>'
        Name of the column in adata.obs that contains cell-type annotations based on which 
        the GEP is supposed to be generated. The default value chooses the last column in the adata.obs
    annot: 'str' | default = 'ENSEMBL'
        Choose which gene annotation to use ['ENSEMBL', 'SYMBOL'] for the exported GEP
    outpath `str` | default = current working directory
        filepath to the directory in which the results should be outputed, if no directory is 
        specified it outputs the results to the current working directory.        

    returns
    -------
    None
        files are written out

    """
    if outpath is None:
        outpath = os.getcwd()
    ### check if the outdir exists if not create
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if (column == 'louvain' or column == 'leiden'):
        labeling = sorted(set(adata.obs[column].astype(int)))
    elif (column == '<last_column>'):
        column = adata.obs.columns[-1]
        labeling = sorted(set(adata.obs[column]))
    else:
        labeling = sorted(set(adata.obs[column]))
    print('Choosing column: [', column, ']for cell annotations')
    
    label_names = []
    for i in range(len(labeling)):
        label_names.append(str(labeling[i]))

    num_labels = len(label_names)

    gct = DataFrame(index=adata.var_names, columns=label_names)
    mydat = adata.raw.var.copy()

    print('Loading gene expression from adata.raw.X')
    if type(adata.raw.X) == ndarray:
        E = sparse.csr_matrix(adata.raw.X).T
    else:
        E = adata.raw.X.tocsr().T
    # mask to subset only highly variable genes
    mask = adata.raw.var_names.isin(adata.var_names)
    mask = arange(len(adata.raw.var_names))[mask]
    #revert to linear scale
    E = E.expm1()

    print('Calculating average expression per cell type per gene')
    for i in range(num_labels):
        cells = where(adata.obs.get(column) == label_names[i])[0]
        #gct.iloc[:,i] = E[:,cells].mean(axis= 1) #get mean expression per gene
        #gct.iloc[:, i] = E[ix_(mask, cells)].mean(axis=1) # get mean expression per *highly variable* gene per cell type
        gct.iloc[:, i] = E[ix_(mask, cells)].sum(axis=1) # get the sum of expression per *highly variable* gene per cell type
        labeling_size = len(cells)

    mydat = mydat.loc[:, [annot]]
    gct = mydat.merge(gct, how='right', left_index=True, right_index=True)

    gct.set_index(annot, inplace=True)
    gct.index.names = ['NAME']


    gctFile_average = os.path.join(outpath, filename)
    gct.to_csv(gctFile_average,
               index=True,
               index_label='NAME',
               header=True,
               float_format='%.3f')
    print('{} exported successfully to file'.format(filename))

    return (None)
    sys.exit(0)

