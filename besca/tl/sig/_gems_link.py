import pytest
import sys
from re import compile, match

from requests import get, post


def insert_gems(BASE_URL, genesets, params, headers=["setName", "desc", "genes"]):
    """Insert genesets into the local gems server
    url_host will depend on GeMs deployement. Could be stored in crendential files.
    Parameters
    ----------
    BASE_URL: class:`str`
        an string  'http://' + hostname + ':' +  localport
    genesets: `list`
        a list of dict; each dict is a signature; key values should mapp the headers
    params: `list` of strings.
        The command-line arguments for GMTx file upload (see below) based on GeMs structure
    headers: `list` of string
        each element is a key of the GEMs setup in place. Minimal requirement for a geneset would be
        setName, desc and genes (minimal GMT)

    Returns
    -------
    results: `int`
        call results 200 indicates a good request

    """
    BASE_URL += "/api/insert"
    if isinstance(genesets, dict):  # only one geneset as dict
        genesets = [genesets]
    parsed = []
    for geneset in genesets:
        sParsed = []
        for header in headers:
            if header != "genes":
                sParsed.append(geneset[header])
            else:
                sParsed += geneset[header]
    parsed.append(sParsed)
    dataIn_1 = {"headers": headers, "parsed": parsed, "params": params}
    returnJSON = post(BASE_URL, json=dataIn_1).json()
    return returnJSON["response"]


def get_gems(setName, BASE_URL, UP_suffix="_UP", DN_suffix="_DN", verbose=False):
    """Connect to GEMS, dowload related geneset (specified by setName, can be a prefix/suffix)
    and return them
    This function combines genesets (signatures) scores (UP and DN) genes.
    Non directionaly geneset are by default considered as UP.
    Parameters
    ----------
    setName: `str`
        setName to find in GeMs (can be a subset)
    BASE_URL: `str`
        GeMS url for the api. Should look like: 'http://' + hostname + ':' +  localport
    UP_suffix : `str` | default = "_UP"
        str suffix indicating that the suffix indicating the signature  is in UP direction.
        This should be the end of the signatures names ($)
    DN_suffix : `str` | default = "_DN"
        str suffix indicating that the suffix indicating the signature is in DN direction.
        This should be the end of the signatures names ($)
    Returns
    -------
    a dictionnary containing the signature names as key, and subdictionnary as direction
    (Key could then be : UP or DN). values are then the gene names.
    Example
    -------
    >>> pytest.skip('Test is only for here as example and should not be executed')
    >>> import yaml
    >>> with open('mongocredentials/credential.yml') as f:
    ...     cred = yaml.safe_load(f)
    >>> get_GEMS_sign('Tcell', BASE_URL =  'http://' + cred['hostname'] + ':' +  cred['localport'])
    """
    # Calling GEMS
    BASE_URL += "/api/genesets"
    GEMS_command = {"setName": setName, "returnParams": ["setName", "desc", "genes"]}
    returnJSON = post(BASE_URL, json=GEMS_command).json()
    # Pattern compilation
    pattern_DN = compile("(\w+)" + DN_suffix + "$")
    pattern_UP = compile("(\w+)" + UP_suffix + "$")
    # Signature dict
    signed_sign = {}
    if "response" not in returnJSON.keys():
        sys.exit(print(returnJSON))
    for i in range(0, len(returnJSON["response"])):
        signature_full_name = returnJSON["response"][i]["setName"]
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
            if signatureName in signed_sign.keys():
                signed_sign[signatureName][direction] = [
                    x[0][0] for x in returnJSON["response"][i]["genes"]
                ]
            else:
                signed_sign[signatureName] = {
                    direction: [x[0][0] for x in returnJSON["response"][i]["genes"]]
                }
            if verbose:
                print(i, ": ", signature_full_name)
    return signed_sign


def get_similar_geneset(
    request, BASE_URL, similarity_coefficient=0.5, method="overlap", outputGeneset=False
):
    """Encapsulating small similary research. Will look for simalirity within GeMs and the mongoDB collections
    and returns the associated geneseets.
    Parameters
    ----------
    request: `string`
        request specificity, if the hosted collection is large, one might need to specify more into details the geneset.
    BASE_URL: `str`
        GeMS url for the api. Should look like: 'http://' + hostname + ':' +  localport
    UP_suffix : `str` | default = "_UP"
        str suffix indicating that the suffix indicating the signature  is in UP direction.
        This should be the end of the signatures names ($)
    DN_suffix : `str` | default = "_DN"
        str suffix indicating that the suffix indicating the signature is in DN direction.
        This should be the end of the signatures names ($)
    Returns
    -------
    a dictionnary containing the signature names as key, and subdictionnary as direction
    (Key could then be : UP or DN). values are then the gene names.
    Example
    -------
    >>> pytest.skip('Test is only for here as example and should not be executed')
    >>> import yaml
    >>> with open('mongocredentials/credential.yml') as f:
    ...     cred = yaml.safe_load(f)
    >>> get_similar_geneset(request='?setName=dz:770_UP&source=CREEDS&user=Public&subtype=disease',
            BASE_URL =  'http://' + cred['hostname'] + ':' +  cred['localport'], outputGeneset = True)
    >>> # this code is only displayed not executed
    """
    if not method in ["jaccard", "overlap"]:
        sys.exit(
            "Method should be  jaccard or overlap. " + method + " is not a valid choice"
        )
    BASE_URL_similar = BASE_URL + "/api/similar"
    request = (
        BASE_URL_similar
        + request
        + "&method="
        + method
        + "&threshold="
        + str(similarity_coefficient)
    )
    returnJSON = get(request).json()
    if "message" in returnJSON.keys():
        print(returnJSON["message"])
    if "response" in returnJSON.keys():
        getChecking = returnJSON["response"]
        if not outputGeneset:
            return getChecking
        else:
            print("Catching all genesets related to the request")
            signature_dict = {}
            for el in getChecking:
                signature_dict.update(get_GEMS_sign(el["setName"], BASE_URL=BASE_URL))
            return signature_dict
