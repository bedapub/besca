from re import compile, match

from requests import post


def insert_GeMs(url_host_insert, genesets, params, headers=['setName', 'desc', 'genes']):
    """Insert genesets into the local gems server
    url_host will depend on GeMs deployement. Could be stored in crendential files.
    Parameters
    ----------
    url_host_insert: class:`str`
        an string  'http://' + hostname + ':' +  localport + '/api/insert'
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
        call results

    """
    if isinstance(genesets, dict):  # only one geneset as dict
        genesets = [genesets]
    parsed = []
    for geneset in genesets:
        sParsed = []
        for header in headers:
            if header != 'genes':
                sParsed.append(geneset[header])
            else:
                sParsed += geneset[header]
    parsed.append(sParsed)
    returnJSON = post(url_host_insert, json=dataIn_1).json()
    return(returnJSON)


def get_GEMS_sign(GEMS_command, BASE_URL,  UP_suffix='_UP', DN_suffix='_DN', verbose=False):
    """ Connect to GEMS, dowload related geneset (specified by GEMS command)
    and return them
    This function combines genesets (signatures) scores (UP and DN) genes.
    Non directionaly geneset are by default considered as UP.
    Parameters
    ----------
    GEMS_command: `dict` | default = None
        dictionary specifying values possible for GEMS
    BASE_URL: `str` 
        GeMS url for the api. Should look like: 'http://' + hostname + ':' +  localport + '/api/genesets'
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
    >>> import yaml
    >>> with open('mongocredentials/credential.yml') as f:
            cred = yaml.safe_load(f)
    >>> dataIn = {
    'setName': 'HumanCD45p',
    'returnParams': ['setName','desc', 'genes']
    }
    >>> get_GEMS_sign(dataIn, BASE_URL =  'http://' + cred['hostname'] + ':' +  cred['localport'] + '/api/genesets')
    >>> # this code is only displayed not executed
    """
    # Calling GEMS
    returnJSON = post(BASE_URL, json=GEMS_command).json()
    # Pattern compilation
    pattern_DN = compile("(\w+)" + DN_suffix + "$")
    pattern_UP = compile("(\w+)" + UP_suffix + "$")
    # Signature dict
    signed_sign = {}
    for i in range(0, len(returnJSON['response'])):
        signature_full_name = returnJSON['response'][i]['setName']
        if len(signature_full_name):
            z = match(pattern_DN, signature_full_name)
            if(z):
                signatureName = z.groups()[0]
                direction = "DN"
            else:
                z = match(pattern_UP, signature_full_name)
                if(z):
                    signatureName = z.groups()[0]
                    direction = "UP"
                else:
                    signatureName = signature_full_name
                    direction = "UP"
            if(signatureName in signed_sign.keys()):
                signed_sign[signatureName][direction] = [x[0][0]
                                                         for x in returnJSON['response'][i]['genes']]
            else:
                signed_sign[signatureName] = {direction: [
                    x[0][0] for x in returnJSON['response'][i]['genes']]}
            if(verbose):
                print(i, ": ", signature_full_name)
    return(signed_sign)

