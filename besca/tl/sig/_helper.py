# this file contains the helper functions
# for signature scoring analysis in python using scanpy


def _to_geneid(conversionTable, symbol):
    """Convert the symbol into another using the conversionSymbol table,
    converting symbol from the first column of the table to the index of such.
    To use to convert from HGNC symobl to Ensembl and vice-versa e.g.

    Parameters
    ----------
    conversionTable:class:`~Serie`
        A serie containing
    symbol:class`~str`: symbol to convert.
        Should match in the first column of conversionTable

    Returns
    -------
    a str, the converted symbol.        

    Example
    -------

    >>> x = pd.Series(['aa', 'bb', 'cc', 'dd', 'ee'], index=['a', 'b', 'c', 'd', 'e'])
    >>> _to_geneid( x, 'bb')
    """
    try:
        res = conversionTable[conversionTable == symbol].index.format()[0]
        return res
    except IndexError:
        return None
