import sys
import re
import pandas as pd
from anndata import AnnData


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


def update_qualitative_palette(
    adata: AnnData,
    palette: dict[str, str],
    group: str = "leiden",
    checkColors: bool = True,
) -> None:
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
        Used internally like this: `pd.Categorical(adata.obs[<group>]).categories.tolist()`
    checkColors: `boolean`
        check the colors inputed to transform them if needed into a hex values. tupple of RBG of 0-1 values cna be converted.
     returns
    -------
    None; update the AnnData object, that the color order matches the order of the AnnData object categories
    """

    # get the groups/categories in the same way as scanpy does it in scanpy/plotting/_tools/scatterplots.py: _get_palette
    category_list = pd.Categorical(adata.obs[group]).categories.tolist()

    # Checking Validity
    if not all(elem in palette.keys() for elem in category_list):
        sys.exit(
            "Please provide a palette dict containing all element of the group " + group
        )
    if checkColors:
        palette = {k: check_colors(color) for k, color in palette.items()}

    newColorList = []
    for category_name in category_list:
        newColorList.append(palette[category_name])

    adata.uns[group + "_colors"] = newColorList.copy()
    return None
