import pytest
import os
from pandas import read_csv


def get_mito_genes(species: str = "human", annotation_type: str = "ENSEMBL"):
    """Returns the array of genes annotated as mitochondrial in species.
    Parameters
    ----------
    species:`str`| default = human ; species of the datasets.
    Accepted: cyno, cynomolgus, human, mouse, rat, pig
    annotation_type:`str`| default = ENSEMBL.  ENSEMBL or SYMBOL accepted

    Returns
    -------
    mito_genes : array of str

    Example
    -------
    >>> pytest.skip('Test does not work')
    >>> import besca as bc
    >>> mito_genes = bc.datasets.get_mito_genes('human')
    >>> mito_genes

    """
    valid = {"cyno", "cynomolgus", "human", "mouse", "rat", "pig"}
    if species not in valid:
        raise ValueError("species must be one of %s." % valid)
    ref_mito_file = os.path.dirname(__file__) + "/mito_files/" + species + ".mito.tsv"
    # ENS_GENE_ID  GENE_SYMBOL (2 cols)
    if annotation_type == "SYMBOL":
        mito_list = list(read_csv(ref_mito_file, header=None, sep="\t")[1])
    elif annotation_type == "ENSEMBL":
        mito_list = list(read_csv(ref_mito_file, header=None, sep="\t")[0])
    else:
        raise ValueError("annotation_type must be either SYMBOL or ENSEMBL")
    return mito_list
