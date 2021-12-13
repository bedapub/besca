# this file contains the metric functions for signature scoring analysis using scanpy
# those functions compute scores for each signature (genesets)
# import using the python version 1.3.2 at least is a prerequisit!
# needs to be checked in all functions

from scanpy.tools  import score_genes


def _handle_signature(signature, method, adata, signature_name,
                      overwrite, verbose, use_raw):
    """Compute signed score combining UP and DN for a speficic signature
    This function combines geneset (signature) scores compose of UP and DN.
    The result is stored in adata.obs with the key:
    "score_" + method + signature_name.
    Multiples methods can be used to compute geneset scores.
    The scanpy method is the score_gene method.
    Combination of the scores is done substracting UP and DN (scanpy = UP - DN)
    Method in development. Not all options implemented yet.
    TODO : Score could be computed while reading the gmt (one loop less).
    However here we divided geneset provenance and computation (to discuss)

    Parameters
    ----------
    adata:class:`~anndata.AnnData`
        An AnnData object (from scanpy).
        Following besca convention, var names (gene) are HGNC symbol. It should match the signatures values.
    signature: `dict`
       dictionnary,; keys are the directions (UP/DN) and values are the genes.
    method: `str` | default = scanpy
        a string indicating which method to use ('scanpy' available)
    use_raw: `boolean`
    Returns
    -------
    None
    The adata object is modified.

    Example
    -------

    >>> #insert example code here
    >>> 1 + 1
    >>> # this code is only displayed not executed
    """
    scoreName = 'score_' + signature_name + "_" + method
    if (verbose):
        print('Computing ' + scoreName )
    if not overwrite:
        if scoreName in adata.obs:
            if (verbose):
                    print(signature_name + ' skipped was already  precomputed with this method (', scoreName, ')')
            return None
    score = [0] * adata.n_obs
    # Not copying columns but rather deleting them afterward is faster. For each direction the column is computed in adata.obs.
    # If the genes are not found, exception to manage as data.obs is not updated
    # The scores are combined and then deleted.
    for direction in signature[signature_name].keys():
            if method == 'scanpy':
                currentName = scoreName + direction
                score_genes(adata, signature[signature_name][direction],
                            score_name=currentName, use_raw=use_raw, copy=False)
                try:
                    scoreTMP = adata.obs[currentName]
                    adata.obs = adata.obs.drop(columns=currentName)
                    if direction == 'DN':
                        scoreTMP = [x * -1 for x in scoreTMP]
                except:
                    scoreTMP = [0] * adata.n_obs
                    if(verbose):
                        print("score " + currentName + " is 0. Exception")
                # SUMMING
                score = list(map(sum, zip(score, scoreTMP)))
    adata.obs[scoreName] = score.copy()
    return(None)
