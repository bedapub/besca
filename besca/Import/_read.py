import os
import sys

import pandas as pd
from anndata import AnnData
from scanpy import read
from scipy.sparse import csr_matrix, issparse

from .._helper import convert_ensembl_to_symbol


def assert_filepath(filepath):
    """Asserts that the filepath contains the files required by read_mtx.

    Parameters
    ----------
    filepath: `str`
        filepath as string to the directory that is expected to contain at least
        following files: `matrix.mtx`, `genes.tsv`, and `barcodes.tsv` OR compressed
        versions of those files, i.e. `matrix.mtx.gz`, `genes.tsv.gz`, and `barcodes.tsv.gz`

    Returns
    -------
    None if the filepath is valid, returns 'gz' if valid and compressed files,
    otherwise an Exception will be raised
    """
    req_files = ["matrix.mtx", "genes.tsv", "barcodes.tsv"]
    req_files_exist = [os.path.isfile(os.path.join(filepath, x)) for x in req_files]
    valid = all(req_files_exist)
    if valid:
        return None

    # Check whether gzipped versions exist:
    req_files = ["matrix.mtx.gz", "genes.tsv.gz", "barcodes.tsv.gz"]
    req_files_exist = [os.path.isfile(os.path.join(filepath, x)) for x in req_files]
    valid = all(req_files_exist)
    if valid:
        return "gz"
    else:
        for ind, exist in enumerate(req_files_exist):
            if not exist:
                raise FileNotFoundError(
                    "{} is not found in "
                    "the given path `{}`".format(req_files[ind], filepath)
                )


def add_var_column(adata, colname="SYMBOL", attempFix=True):
    if not colname in adata.var.columns:
        if attempFix:
            print(
                f"Creating empty {colname} column;  please check if {colname} is not in the index of adata.var"
            )
            adata.var[colname] = "N.A."
        else:
            raise Exception("adata.var needs column named {colname}")
    return adata


def assert_adata(adata: AnnData, attempFix=True):
    """Asserts that an adata object is containing information needed for the besca pipeline to run and export information.
    This is particularly usefull when loading public data
    The parameter attempFix will try to fix the issue by itself.
    However, we advise the user to check by himself what is the leading problem.

    Parameters
    ----------
    adata: AnnData
    attempFix: `bool` if True will transform adata object to match requirements.
    Returns
    -------
    returns an AnnData object
    """
    if not "CELL" in adata.obs.columns:
        if attempFix:
            adata.obs["CELL"] = adata.obs.index
            print("Creating columns CELL in adata.obs using adata.obs.index.")
        else:
            raise Exception("Required CELL columns in observations")

    if not all(adata.obs_names == adata.obs["CELL"]):
        raise Exception("Required indexing of adata.obs by CELL column")
    if not issparse(adata.X):
        if attempFix:
            print("Required count matrix to be sparse, X transformed to sparse")
            try:
                adata.X = sparse.csr_matrix(adata.X.copy())
            except:
                raise Exception("X transformation to sparse failed.")
        else:
            raise Exception("adata.X needs to be sparse.")
    # checking adata.var concordance
    for x in ["SYMBOL", "ENSEMBL"]:
        adata = add_var_column(adata, x, attempFix)
        if not all(isinstance(el, str) for el in adata.var.get(x)):
            raise Exception("In {x} non string values will create an issue for export")
    return adata


def read_mtx(
    filepath, annotation=True, use_genes="SYMBOL", species="human", citeseq=None
):
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
    gzfiles = assert_filepath(filepath)
    if gzfiles == "gz":
        print("reading matrix.mtx.gz")
        adata = read(
            os.path.join(filepath, "matrix.mtx.gz"), cache=True
        ).T  # transpose the data
        print("adding cell barcodes")
        adata.obs_names = pd.read_csv(
            os.path.join(filepath, "barcodes.tsv.gz"), compression="gzip", header=None
        )[0]
        print("adding genes")
        var_anno = pd.read_csv(
            os.path.join(filepath, "genes.tsv.gz"),
            compression="gzip",
            header=None,
            sep="\\t",
            engine="python",
        )
    else:
        print("reading matrix.mtx")
        adata = read(
            os.path.join(filepath, "matrix.mtx"), cache=True
        ).T  # transpose the data
        print("adding cell barcodes")
        adata.obs_names = pd.read_csv(
            os.path.join(filepath, "barcodes.tsv"), header=None
        )[0]
        print("adding genes")
        var_anno = pd.read_csv(
            os.path.join(filepath, "genes.tsv"), header=None, sep="\\t", engine="python"
        )

    symbols = var_anno[1]
    ensembl_id = var_anno[0]
    adata.var["ENSEMBL"] = ensembl_id.tolist()
    adata.var.index.names = ["index"]

    if use_genes == "SYMBOL":
        adata.var_names = symbols

        # make unique
        print("making var_names unique")
        adata.var_names_make_unique()

        if symbols.tolist() != ensembl_id.tolist():
            print("adding ENSEMBL gene ids to adata.var")
            adata.var["SYMBOL"] = symbols.tolist()
    elif use_genes == "ENSEMBL":
        adata.var_names = ensembl_id

        if symbols[0] != ensembl_id[0]:
            # add symbols to object
            print("adding symbols to adata.var")
            adata.var["SYMBOL"] = symbols.tolist()
        else:
            # lookup the corresponding Symbols
            adata.var["SYMBOL"] = convert_ensembl_to_symbol(
                ensembl_id.tolist(), species=species
            )
    else:
        sys.exit("supplied unknown use_genes parameter")

    if citeseq is not None:
        features = var_anno[2]
        adata.var["feature_type"] = features.tolist()

        # only extract the information we want
        if citeseq == "gex_only":
            adata = adata[:, adata.var.feature_type == "Gene Expression"].copy()
        if citeseq == "citeseq_only":
            adata = adata[:, adata.var.feature_type == "Antibody Capture"].copy()

    if annotation == True:
        print("adding annotation")
        adata.obs = pd.read_csv(
            os.path.join(filepath, "metadata.tsv"), sep="\\t", engine="python"
        )
        if adata.obs.get("CELL") is not None:
            adata.obs.index = adata.obs.get("CELL").tolist()

    return adata
