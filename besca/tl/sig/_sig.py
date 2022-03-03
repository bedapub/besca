# this file contains the main functions for signature scoring analysis in python using scanpy

import sys
import logging

# for conversion
from itertools import repeat

from besca.tl.sig._helper import _to_geneid
from besca.tl.sig._io_sig import read_GMT_sign
from besca.tl.sig._metrics import _handle_signature

def filter_by_set(strs, universe_set):
    """Remove strings from the list that are not in the universe set
    
   Parameters
    ----------
    strs: a `list` or `tuple` or `set` of `str`
        A sequence (ordered list) or unordered set of strings, which often
        are HGNC symbol and should match the signatures values. Empty 
        characters around characters are ignored.
    universe_set: a `set` of strings
        a set of strings. The strings in `strs` that are not in this set 
        are filtered.

    Returns
    -------
    filter_by_set: `list` of `str`
        A list of gene names that are detected, in the same order as the input
        except for those that are filtered. In case the input is None or an 
        empty list, an empty list is returned.
        
    Example
    -------

    >>> import besca as bc
    >>> detected = list('ABCDE')
    >>> bc.tl.sig.filter_by_set(['D', 'E', 'A'], detected)
    ['D', 'E', 'A']
    >>> bc.tl.sig.filter_by_set(['D', 'E', 'A', 'F'], detected)
    ['D', 'E', 'A']
    >>> bc.tl.sig.filter_by_set(None, detected)
    []
    >>> bc.tl.sig.filter_by_set([], detected)
    []
    """
    if strs is None or len(strs)==0:
        return []
    genes = [gene.strip() for gene in strs] ## sometimes the gene names have empty spaces around
    int_gene_set = set(genes).intersection(universe_set)
    int_gene = sorted(int_gene_set, key=genes.index) ## keep the input order
    return(int_gene)
                
    
def filter_siggenes(adata, signature_dict):
    """Filter all signatures in signature_dict to remove genes not present in adata

    Parameters
    ----------
    adata: class:`~anndata.AnnData`
        An AnnData object (from scanpy). Following besca convention, var names
        (genes) are HGNC symbol and should match the signatures values.
    signature_dict: `dict`
        a dictionary of signatures. Keys signature names, and values can come 
        in two variants:
        * Variant 1: values are a list of gene names as strings. An example: 
        `{'gs1': ['A', 'B'], 'gs2':['B', 'C']}`;
        * Variant 2: Values are a dict with keys, for instance directions 
        (UP/DN), and genes names in values. An example: `{'gs1': {'UP': 'A', 
        'DN': 'B'}, 'gs2': {'UP': ['C', 'D'], 'DN': ['E', 'A']}}`.

    Returns
    -------
    signature_dict_filtered: `dict`
        a dictionary of signatures in the same format as the input.
        
    Example
    -------

    >>> import besca as bc
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> sigs = {'GeneSet1': ['JUNB', 'ALAD', 'ZNF559', 'NoSuchAGene'],
    ...         'GeneSet2': ['ITGA1', 'CTCF', 'CAPN10', 'UnknownGene']}
    >>> filtered_sigs = bc.tl.sig.filter_siggenes(adata, sigs)
    >>> filtered_sigs
    {'GeneSet1': ['JUNB', 'ALAD', 'ZNF559'], 'GeneSet2': ['ITGA1', 'CTCF', 'CAPN10']}
    >>> signed_sigs = {'GeneSet1': {'UP': ['JUNB', 'ALAD', 'ZNF559',
    ...                                              'NoSuchAGene'],
    ...                               'DN': ['REM2', 'AKT1', 'GPD2', 
    ...                                              'UnknownGene']},
    ...                  'GeneSet2': {'UP': ['TOP2A', 'TOP3A', 'MDM2', 
    ...                                               'NoSuchAGene2'],
    ...                               'DN': ['ITGA1', 'CTCF', 'CAPN10', 
    ...                                               'UnknownGene2']}}
    >>> filtered_signed_sigs = bc.tl.sig.filter_siggenes(adata, signed_sigs)
    >>> filtered_signed_sigs
    {'GeneSet1': {'UP': ['JUNB', 'ALAD', 'ZNF559'], 'DN': ['REM2', 'AKT1', 'GPD2']}, 'GeneSet2': {'UP': ['TOP2A', 'TOP3A', 'MDM2'], 'DN': ['ITGA1', 'CTCF', 'CAPN10']}}
    """

    signature_dict_filtered = {}
    raw_index_set = set(adata.raw.var.index)
    for geneset, dir_dict in signature_dict.items():
        if type(dir_dict) is dict:
            signature_dict_filtered[geneset] = {}
            for direction, genes in dir_dict.items():
                int_gene = filter_by_set(genes, raw_index_set)
                if len(int_gene) >= 1:
                    signature_dict_filtered[geneset][direction] = int_gene
                else:
                    logging.info('No genes are left after filtering in '
                          + geneset + ' direction ' + direction)
        elif type(dir_dict) is list or type(dir_dict) is tuple:
            int_gene = filter_by_set(dir_dict, raw_index_set)
            if len(int_gene) >= 1:
                signature_dict_filtered[geneset] = int_gene
            else:
                logging.info('No genes are left after filtering in '
                      + geneset)
        else:
            raise ValueError('The signature_dict is not in a valid format.')
    return signature_dict_filtered


def combined_signature_score(
    adata,
    GMT_file=None,
    signature_dict=None,
    UP_suffix="_UP",
    DN_suffix="_DN",
    method="scanpy",
    overwrite=False,
    verbose=False,
    use_raw=None,
    conversion=None,
):
    """Super Wrapper function to compute combined signature score for UP and DN scores.
    This function combines genesets (signatures) scores compose of UP and DN.
    Results are stored in adata.obs with the key: "score_"+ signature_name+"_" + method  .
    The scanpy method is the score_gene method from the scanpy python package.
    Combination of the scores is done substracting UP and DN (scanpy = UP - DN).

    If you do not have a signature dictionary composed with direction; please see bc.tl.sig.convert_to_directed

    Parameters
    ----------
    adata: class:`~anndata.AnnData`
        An AnnData object (from scanpy). Following besca convention,
        var names (gene) are HGNC symbol and should match the signatures values.
    GMT_file: `str` | default = None
        gmt file location containing the geneset
    signature_dict: `str` | default = None
        pre-loaded signature dictionnary using read_GMT_sign or get_GEMS_sign
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
    use_raw: `boolean` | default = None
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

    >>> import os
    >>> import besca as bc
    >>> bescapath = os.path.split(os.path.dirname(bc.__file__))[0]
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> gmt_file= bescapath + '/besca/datasets/genesets/Immune.gmt'
    >>> bc.tl.sig.combined_signature_score( adata, GMT_file = gmt_file)
    >>> # this code is only displayed not executed

    """
    # RMK : Score could be computed while reading the gmt (one loop less).
    # However here we divided geneset provenance and computation.
    if GMT_file is None and signature_dict is None:
        sys.exit("need to provide either GMT_file or signature_dict gene annotation.")
    if not GMT_file is None:
        if signature_dict is not None:
            signature_dict.update(
                read_GMT_sign(GMT_file, UP_suffix, DN_suffix, True, verbose)
            )
        else:
            signature_dict = read_GMT_sign(
                GMT_file, UP_suffix, DN_suffix, True, verbose
            )
    if verbose:
        print(str(len(signature_dict)) + " signatures obtained")
    # More readable than in signatures read. This forces a second loop.
    #  Should be optimized later on
    if conversion is not None:
        for signature in signature_dict.keys():
            for direction in signature_dict[signature]:
                signature_dict[signature][direction] = [
                    i
                    for i in map(
                        _to_geneid,
                        repeat(conversion),
                        signature_dict[signature][direction],
                    )
                    if i is not None
                ]
    # compute signed score
    compute_signed_score(
        adata=adata,
        signature_dict=signature_dict,
        method=method,
        overwrite=overwrite,
        verbose=verbose,
        use_raw=use_raw,
    )
    return None


def compute_signed_score(
    adata, signature_dict, method="scanpy", overwrite=False, verbose=False, use_raw=None
):
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
    use_raw: `boolean` | default = None

    Returns
    -------
    None

    """
    # Filter out signature genes not present in adata
    signature_dict_filtered=filter_siggenes(adata, signature_dict)

    [
        _handle_signature(
            signature_dict_filtered, method, adata, signature_name, overwrite, verbose, use_raw
        )
        for signature_name in signature_dict.keys()
    ]
    return None


def make_gmtx(
    setName,
    desc,
    User,
    Source,
    Subtype,
    domain,
    genesetname,
    genes,
    studyID=None,
    analysisID=None,
    application=None,
    celltype=None,
    coef_type="logFC",
):
    """Construct a gmtx file according to format conventions for import into Gems.
    Parameters
    ----------
    setName: `str`
        informative set name e.g. Pembro_induced_MC38CD8Tcell, Plasma_mdb, TGFB_Stromal_i
    desc : `str`
        informative and verbose signature description; for cell type signatures use nomenclature,
        if coef used explain what it represents; link to study if present; e.g. Genes higher expressed in
        Pembro vs. vehicle in non-naive CD8-positive T cells in MC38 in vivo exp. ID time T2; coefs are log2FC
    User : `str`
        related to signature origin e.g. Public (for literature-derived sets), own user ID for analysis-derived sets,
        rtsquad, scsquad, gred, other
    Source : `str`
        source of the signature,  one of Literature scseq, Literature, besca, scseqmongodb,
        internal scseq, pRED, Chugai, gRED, other
    Subtype : `str`
        specific subtype e.g. onc, all, healthy, disease
    domain : `str`
        one of pathway, biological process, cellular component,molecular function, phenotype,
        perturbation, disease, misc, microRNA targets, transcription factor targets, cell marker, tissue marker
    genesetname: `str`
        shared across different signatures of a specific type e.g. besca_marker, dblabel_marker,
        Pembro_induced_MC38CD8Tcell, FirstAuthorYearPublication
    genes: `str`
        tab-separated list of genes with/without a coefficient e.g. Vim | 2.4\tBin1 | 2.02 or Vim\tBin1
    studyID: `str` | default = None
        study name as in scMongoDB/bescaviz; only when source=internal scseq
    analysisID: `str` | default = None
        analysis name as in scMongoDB/bescaviz; only when source=internal scseq
    application: `str` | default = None
        specify which application will read the geneset e.g. rtbeda_CIT, bescaviz, celltypeviz
    celltype: `str` | default = None
        for cell markers, specify celltype according to dblabel_short convention to facilitate reuse
    coef_type: `str` | default = score
        specify what the coefficient corresponds too, e.g. logFC, gini, SAM, score, ...

    Returns
    -------
    geneset
        a dictionary with populated fields needed to later export the signature to gmtx

    Example
    -------
    >>> import besca as bc
    >>> User = 'nouser'
    >>> Source = 'pbmc3k_processed'
    >>> Subtype = 'public'
    >>> domain = 'perturbation'
    >>> studyID = 'pbmc3k_processed'
    >>> analysisID = 'default'
    >>> genesetname = 'pbmc3k_processed_cluster0'
    >>> setName = 'pbmc3k_processed_cluster0'
    >>> desc = 'Genes higher expressed in cluster 0; coefs are log2FC'
    >>> adata = bc.datasets.pbmc3k_processed()
    >>> myfc = 1
    >>> mypval = 0.05
    >>> DEgenes=bc.tl.dge.get_de(adata,'leiden',demethod='wilcoxon',topnr=5000, logfc=myfc,padj=mypval)
    >>> pdout=DEgenes['0'].sort_values('Log2FC', ascending=False)
    >>> genes="\t".join(list(pdout['Name'].astype(str) + " | " + pdout['Log2FC'].round(2).astype(str)))
    >>> signature_dict = bc.tl.sig.make_gmtx(setName,desc,User,Source,Subtype,domain,genesetname,genes,studyID,analysisID)
    Prefered source names: Literature scseq, Literature, besca, scseqmongodb, internal scseq, pRED, Chugai, gRED, other
    Metadata for signature pbmc3k_processed_cluster0 successfully captured
    """

    geneset = {}
    geneset["setName"] = setName
    geneset["desc"] = desc
    geneset["User"] = User
    geneset["Source"] = Source
    geneset["Subtype"] = Subtype
    geneset["geneset"] = genesetname

    if domain in [
        "pathway",
        "biological process",
        "cellular component",
        "molecular function",
        "phenotype",
        "perturbation",
        "disease",
        "misc",
        "microRNA targets",
        "transcription factor targets",
        "cell marker",
        "tissue marker",
    ]:
        geneset["domain"] = domain
    else:
        if domain == None:
            domain = "misc"
            print("You did not specify a domain name, it will be set to misc")
        print(
            "Recommended domain names are pathway, biological process, cellular component, molecular function, phenotype, perturbation, disease, misc, microRNA targets, transcription factor targets, cell marker, tissue marker"
        )
        geneset["domain"] = domain

    if Source == "internal scseq":
        if studyID == None:
            print(
                "Signatures of type "
                + Source
                + " require a studyID, please provide one."
            )
        else:
            geneset["studyID"] = studyID
        if analysisID == None:
            print(
                "Signatures of type "
                + Source
                + " require a analysisID, please provide one."
            )
        else:
            geneset["analysisID"] = analysisID
    elif Source == "besca":
        geneset["geneset"] = "besca_marker"
        geneset["domain"] = "cell marker"
        geneset["application"] = "rtbeda_CIT, bescaviz, celltypeviz"
    elif not Source in [
        "Literature scseq",
        "Literature",
        "besca",
        "scseqmongodb",
        "internal scseq",
        "pRED",
        "Chugai",
        "gRED",
        "other",
    ]:
        print(
            "Prefered source names: Literature scseq, Literature, besca, scseqmongodb, internal scseq, pRED, Chugai, gRED, other"
        )
    # pd1il2vsigsdetails['application']=''

    if domain == "cell marker":
        if celltype == None:
            print(
                setName
                + "is missing a celltype. Please specify celltype according to dblabel_short convention."
            )

    if "|" in genes:
        geneset["genes | " + coef_type] = genes
    else:
        geneset["genes"] = genes

    print("Metadata for signature " + setName + " successfully captured")
    return geneset

if __name__ == "__main__":
    import doctest
    doctest.testmod()
