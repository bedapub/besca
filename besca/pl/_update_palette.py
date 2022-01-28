import sys
import re


def check_colors(aColor):
    """
    convert the color given in hex if needed. This avoid warning message a posteriori
    parameters
    ----------
    aColor: ``
      color to check; expected tupple. Hex would be returned as input

    returns
    -------
       the color in hex
    """
    if isinstance(aColor, tuple) and len(aColor) == 3:
        if aColor < (1, 1, 1):
            r = round(aColor[0] * 255)
            g = round(aColor[1] * 255)
            b = round(aColor[2] * 255)
        else:  # assuming rgb
            r = aColor[0]
            g = aColor[1]
            b = aColor[1]
        return "#{:02x}{:02x}{:02x}".format(r, g, b)
    else:
        matchingHex = re.search(r"^#(?:[0-9a-fA-F]{3}){1,2}$", aColor)
        if matchingHex:
            return aColor
        else:
            sys.exit("Color " + str(aColor) + "could not be converted")


def update_qualitative_palette(adata, palette, group="leiden", checkColors=True):
    """Update adata object such that the umap will adhere to the palette provided.
    

    parameters
    ----------
    adata: `AnnData`
      the AnnData object 
    palette: `dict`
        dict with keys as the values of the group observation. To avoid warning from matlib it is advised to have \
            hex color values
    group: `str`
        string identifying the column name of adata.obs where colors will be set.
    checkColors: `boolean`
        check the colors inputed to transform them if needed into a hex values. tupple of RBG of 0-1 values cna be converted. 
     returns
    -------
    None; update the adata object.
    """
    cell_pop = set(adata.obs.get(group))
    # Checking Validity
    if not all(elem in palette.keys() for elem in cell_pop):
        sys.exit(
            "Please provide a palette dict containing all element of the group " + group
        )
    if checkColors:
        palette = {k: check_colors(color) for k, color in palette.items()}
    sortedKey = sorted(k for k in palette.keys() if k in cell_pop)
    sortedKey.reverse()
    newColors = [palette[i] for i in sortedKey]
    newColors.reverse()
    adata.uns[group + "_colors"] = newColors.copy()
    return None
