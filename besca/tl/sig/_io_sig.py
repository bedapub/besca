# this file contains the functions to read / import signatures
# for signature score computations
import sys
import logging
from re import compile, match

from requests import post

def convert_to_directed(signature_dict, direction="UP"):
    """ Convert a simple dictionary into one with direction compatible with combined_signature_score

    Parameters
    ----------
    signature_dict: `str` | default = None
        gmt file location containing the geneset
    direction : `str` | default = "UP"
        str suffix indicating that the suffix indicating the signature direction. 
        Expected UP or DN

    Returns
    -------
    a dictionnary containing the signature names as key. \
    Value are a subdictionnary where key are direction(UP or DN).\

    Example
    -------
    >>> import os
    >>> import besca as bc
    >>> bescapath = os.path.split(os.path.dirname(bc.__file__))[0]
    >>> gmt_file= bescapath + '/besca/datasets/genesets/Immune.gmt'
    >>> mymarkers = bc.tl.sig.read_GMT_sign(gmt_file,directed=False)
    >>> bc.tl.sig.convert_to_directed( mymarkers, 'UP')
    {'lymphocyte': {'UP': ['PTPRC']}, 'myeloid': {'UP': ['S100A8', 'S100A9', 'CST3']}, 'Bcell': {'UP': ['CD19', 'CD79A', 'MS4A1']}, 'Tcells': {'UP': ['CD3E', 'CD3G', 'CD3D']}, 'CD4': {'UP': ['CD4']}, 'CD8': {'UP': ['CD8A', 'CD8B']}, 'NKcell': {'UP': ['NKG7', 'GNLY', 'NCAM1']}, 'monocyte': {'UP': ['CST3', 'CSF1R', 'ITGAM', 'CD14', 'FCGR3A', 'FCGR3B']}, 'macrophage': {'UP': ['CD14', 'IL1B', 'LYZ', 'CD163', 'ITGAX', 'CD68', 'CSF1R', 'FCGR3A']}}
    """
    if (direction != "UP") and (direction != "DN"):
        sys.exit("expecting direction UP or DN for directed signature dictionnary.")
    signed_sign = {}
    for signature, genes in signature_dict.items():
        signed_sign[signature] = {direction: genes}
    return signed_sign


def read_GMT_sign(
    GMT_file, UP_suffix="_UP", DN_suffix="_DN", directed=True, verbose=False
):
    """ Read gmt file to extract signed genesets.
    This function combines genesets scores composed of
    UP and DN regulated genes.
    Non directional geneset are by default considered as UP.

    Parameters
    ----------
    GMT_file: `str` | default = None
        gmt file location containing the geneset
    UP_suffix : `str` | default = "_UP"
        str suffix indicating that the suffix indicating the signature is UP.
        This should be the end of the signatures names ($).
        Indicate a dummy string to avoid combination.
    DN_suffix : `str` | default = "_DN"
        str suffix indicating that the suffix indicating the signature is DN.
        This should be the end of the signatures names ($).
        Indicate a dummy string to avoid combination.
    Returns
    -------
    a dictionnary containing the signature names as key. \
    Value are a subdictionnary where key are direction(UP or DN).\
    Values are then the gene names.
    Example
    -------
    >>> import besca as bc
    >>> import pkg_resources
    >>> gmt_file='datasets/genesets/Immune.gmt' # provided in besca
    >>> gmt_file_abs_path=pkg_resources.resource_filename('besca', gmt_file)
    >>> bc.tl.sig.read_GMT_sign(gmt_file_abs_path)
    {'lymphocyte': {'UP': ['PTPRC']}, 'myeloid': {'UP': ['S100A8', 'S100A9', 'CST3']}, 'Bcell': {'UP': ['CD19', 'CD79A', 'MS4A1']}, 'Tcells': {'UP': ['CD3E', 'CD3G', 'CD3D']}, 'CD4': {'UP': ['CD4']}, 'CD8': {'UP': ['CD8A', 'CD8B']}, 'NKcell': {'UP': ['NKG7', 'GNLY', 'NCAM1']}, 'monocyte': {'UP': ['CST3', 'CSF1R', 'ITGAM', 'CD14', 'FCGR3A', 'FCGR3B']}, 'macrophage': {'UP': ['CD14', 'IL1B', 'LYZ', 'CD163', 'ITGAX', 'CD68', 'CSF1R', 'FCGR3A']}}

    """
    signFile = open(GMT_file, "r")
    text_gmt = signFile.read().split("\n")
    signFile.close()
    signed_sign = {}
    # Here \S is used as signature might have '-' in their name
    #  (\w is not suficient if number in signature for EX.)
    pattern_DN = compile("(\S+)" + DN_suffix + "$")
    pattern_UP = compile("(\S+)" + UP_suffix + "$")
    # TODO: remove this for loop.
    for i in range(0, len(text_gmt)):
        temp_split = text_gmt[i].split("\t")
        signature_full_name = temp_split[0]
        if len(temp_split) < 3:
            if verbose:
                print(
                    "Skipping empty entry; less than 3 fields for "
                    + signature_full_name
                )
            continue
        # Skipping empty lines in gmt files
        if len(signature_full_name):
            z = match(pattern_DN, signature_full_name)
            if z:
                signatureName = z.groups()[0]
                direction = "DN"
            else:
                z = match(pattern_UP, signature_full_name)
                if z:
                    signatureName = z.groups()[0]
                    direction = "UP"
                else:
                    signatureName = signature_full_name
                    direction = "UP"
            # Get the gene names removing empty entry
            initialValue = temp_split[2 : len(temp_split)]
            geneArray = [x for x in initialValue if len(x)]
            if signatureName in signed_sign.keys():
                signed_sign[signatureName][direction] = geneArray
            else:
                signed_sign[signatureName] = {direction: geneArray}
            if verbose:
                print(i, ": ", signature_full_name)

    ### remove UP in case one just wants the signature as a dictionary
    if directed == False:
        mymarkersk = {}
        for key, value in signed_sign.items():
            mymarkersk[key] = value["UP"]
        signed_sign = mymarkersk.copy()
    return signed_sign


def write_gmtx_forgems(signature_dict, GMT_file):
    """Writes a gmtx file that can later be uploaded to GeMS.
    The input should be standardised, as facilitated by bc.tl.sig.make_gmtx.

    Parameters
    ----------
    signature_dict : `dictionary`
        the dictionary containing the signatures to be outputed to a gmtx file
        prepared with bc.tl.sig.make_gmtx, each signature is itself a dictionary

    GMT_file: `str`
        gmt file location containing the geneset


    Example
    -------
    >>> signature_dict = {'setName': 'Th17Tcell_mc38_user', 'desc': 'T-helper 17 cell markers; coefs are log2FC', 'User': 'user', 'Source': 'internal scseq', 'Subtype': 'onc', 'geneset': 'Pembro_MC38-tumor_dblabel', 'domain': 'cell marker', 'studyID': 'Pembro_MC38-tumor', 'analysisID': 'sw_besca24', 'genes|score': 'Cd163l1 | 10.67\tGm9961 | 10.49\tCdh10'}
    >>> outgmtfile='Celltypemarkers.gmtx'
    >>> write_gmtx_forgems(signature_dict, outgmtfile)

    """

    with open(GMT_file, "w") as f:
        f.writelines('\t'.join(list(signature_dict.keys()))+'\n')
        f.writelines('\t'.join(list(signature_dict.values()))+'\n')

    logging.info("Successfully written all signatures to " + GMT_file)
