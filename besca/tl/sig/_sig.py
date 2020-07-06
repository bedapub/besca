# this file contains the main functions for signature scoring analysis in python using scanpy
# import using the python version 1.3.2 at least is a prerequisit! needs to be checked in all functions

# for conversion
from itertools import repeat

from ._helper import _to_geneid
from ._io_sig import read_GMT_sign
from ._metrics import _handle_signature


def filter_siggenes(adata, signature_dict):
    """Filter all signatures in signature_dict to remove genes not present in adata 

    Parameters
    ----------
    adata: class:`~anndata.AnnData`
        An AnnData object (from scanpy). Following besca convention, var names (gene) are HGNC symbol and should match the signatures values.
    signature_dict: `dict`
        a dictionary of signatures. Nested dictionnaries, key: signature names
        Values are a dict  with keys as the directions (UP/DN) and genes names in values.

    Returns
    -------
    signature_dict_filtered: `dict`
        a dictionary of signatures. Nested dictionnaries, key: signature names
        Values are a dict  with keys as the directions (UP/DN) and genes names in values.

    """

    # to add: which signatures were discarded and print (if verbose) genes not found
    # to check: what happens if no signature is left?

    signature_dict_filtered = {}
    for key, value in signature_dict.items():
        mym = []
        for item in value:
            if (sum(adata.raw.var.index.isin([item.strip()])*1) > 0):
                mym.append(item.strip())
        if (len(mym) > 1):
            signature_dict_filtered[key] = mym
    return signature_dict_filtered


def combined_signature_score(adata, GMT_file,
                             UP_suffix='_UP', DN_suffix='_DN', method='scanpy',
                             overwrite=False, verbose=False,
                             use_raw=True, conversion=None):
    """Super Wrapper function to compute combined signature score for UP and DN scores.
    This function combines genesets (signatures) scores compose of UP and DN.
    Results are stored in adata.obs with the key: "score_"+ signature_name+"_" + method  .
    The scanpy method is the score_gene method from the scanpy python package.
    Combination of the scores is done substracting UP and DN (scanpy = UP - DN).

    Parameters
    ----------
    adata: class:`~anndata.AnnData`
        An AnnData object (from scanpy). Following besca convention, 
        var names (gene) are HGNC symbol and should match the signatures values.
    GMT_file: `str` | default = None
        gmt file location containing the geneset
    UP_suffix : str` | default = "_UP"
        str suffix indicating that the suffix indicating the signature in a UP direction (end of the signature).
        Can be replaced by "None" (quoted) or any kind of unexpected string to avoid combination.
    DN_suffix : str` | default = "_DN"
        str suffix indicating that the suffix indicating the signature in a DN direction (end of the signature). 
        Can be replaced by "None" (quoted) or any kind of unexpected string  to avoid combination.
    method: `str` | default = "scanpy"
        a string indicating which method to use ('scanpy' available)
    overwrite: `boolean` | default = False
        If False, will parse the data.obs to only recompute scores that are not present.
    verbose: `boolean` | default = False
        If True, will print the signature reads. This does not overwrite the scanpy verbosity parameter that should be set separately
    use_raw: `boolean` | default = False
        If True, computation will be done on adata.raw.X (on adata.X otherwise).
    conversion: `a panda serie' | default = None.
        If not none, this should contain a serie indexed by x with column containing values y.
        This indicate the how to transpose the signatures (from x to y).
        Classical example would be indexed by Ensembl id and column would contains HGNC symbol.

    Returns
    -------
    None (the adata obs is modified within the function) 
    Example
    -------

    >>> #insert example code here
    >>> adata = bc.datasets.pbmc3k_filtered()
    >>> gmt_file=   'genesets/Immune.txt' # Provided by BESCA datasets check the path
    >>> combined_signature_score( adata, GMT_file = gmt_file)
    >>> # this code is only displayed not executed

    """
    # RMK : Score could be computed while reading the gmt (one loop less).
    # However here we divided geneset provenance and computation.
    signature_dict = read_GMT_sign(
        GMT_file, UP_suffix, DN_suffix, True, verbose)
    if (verbose):
        print(str(len(signature_dict)) + " signatures obtained")
    # More readable than in signatures read. This forces a second loop.
    #  Should be optimized later on
    if (conversion is not None):
        for signature in signature_dict.keys():
            for direction in signature_dict[signature]:
                signature_dict[signature][direction] = [i for i in map(_to_geneid, repeat(
                    conversion), signature_dict[signature][direction]) if i is not None]
    # Filter out signature genes not present in adata
    # signature_dict=filter_siggenes(adata, signature_dict)
    compute_signed_score(adata=adata, signature_dict=signature_dict,
                         method=method, overwrite=overwrite, verbose=verbose, use_raw=use_raw)
    return None


def compute_signed_score(adata, signature_dict, method='scanpy',
                         overwrite=False, verbose=False, use_raw=False):
    """Compute signed score combining UP and DN for all signatures in signature_dict
    This function combines genesets (signatures) scores.
    Results are stored in adata.obs with the key: "score_" + method + signature_name.
    Multiples methods can be used to compute geneset scores.
    The scanpy method is the score_gene method. 
    Combination of the scores is done substracting UP and DN (scanpy = UP - DN).
    Method in development. Not all options implemented yet.

    Parameters
    ----------
    adata: class:`~anndata.AnnData`
        An AnnData object (from scanpy). Following besca convention, var names (gene) are HGNC symbol and should match the signatures values.
    signature_dict: `dict`
        a dictionary of signature. Nested dictionnaries, key: signature names
        Values are a dict  with keys as the directions (UP/DN) and genes names in values.
    method: `str` | default = scanpy
        a string indicating which method to use ('scanpy' available)
    overwrite: `boolean` | default = False
        If False, will parse the data.obs to only recompute scores that are not present.
    use_raw: `boolean` | default = False

    Returns
    -------
    None

    """
    # Filter out signature genes not present in adata
    # signature_dict=filter_siggenes(adata, signature_dict)

    [_handle_signature(signature_dict, method, adata, signature_name, overwrite,
                       verbose, use_raw) for signature_name in signature_dict.keys()]
    return None
