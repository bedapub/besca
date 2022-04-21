import numpy as np
from scipy.sparse.csr import csr_matrix
from anndata._core.views import SparseCSRView

import pytest


def closure(mat):
    """
    Performs closure to ensure that all elements add up to 1.

    Parameters
    ----------
    mat : array_like
       a matrix of proportions where
       rows = compositions
       columns = components
    Returns
    -------
    array_like, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1
    Raises
    ------
    ValueError
       Raises an error if any values are negative.
    ValueError
       Raises an error if the matrix has more than 2 dimension.
    ValueError
       Raises an error if there is a row that has all zeros.
    Examples
    --------
    >>> import numpy as np
    >>> import skbio
    >>> from skbio.stats.composition import closure
    >>> X = np.array([[2, 2, 6], [4, 4, 2]])
    >>> closure(X)
    array([[0.2, 0.2, 0.6],
           [0.4, 0.4, 0.2]])
    """
    mat = np.atleast_2d(mat)
    if np.any(mat < 0):
        raise ValueError("Cannot have negative proportions")
    if mat.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")
    if np.all(mat == 0, axis=1).sum() > 0:
        raise ValueError("Input matrix cannot have rows with all zeros")
    mat = mat / mat.sum(axis=1, keepdims=True)
    return mat.squeeze()


def clr(mat):
    r"""
    Performs centre log ratio transformation.
    This function transforms compositions from Aitchison geometry to
    the real space. The :math:`clr` transform is both an isometry and an
    isomorphism defined on the following spaces
    :math:`clr: S^D \rightarrow U`
    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`
    It is defined for a composition :math:`x` as follows:
    .. math::
        clr(x) = \ln\left[\frac{x_1}{g_m(x)}, \ldots, \frac{x_D}{g_m(x)}\right]
    where :math:`g_m(x) = (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.

    This function is based on the implementation of the clr function within the skbio package.

    Parameters
    ----------
    mat : array_like, float
       a matrix of proportions where
       rows = compositions and
       columns = components
    Returns
    -------
    numpy.ndarray
         clr transformed matrix
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])
    """
    mat = closure(mat)
    lmat = np.log(mat)
    gm = lmat.mean(axis=-1, keepdims=True)
    return (lmat - gm).squeeze()


def multiplicative_replacement(mat, delta=None):
    r"""Replace all zeros with small non-zero values
    It uses the multiplicative replacement strategy [1]_ ,
    replacing zeros with a small positive :math:`\delta`
    and ensuring that the compositions still add up to 1.

    This function is based on the implementation of the multiplicative_replacement function within the skbio package.

    Parameters
    ----------
    mat: array_like
       a matrix of proportions where
       rows = compositions and
       columns = components
    delta: float, optional
       a small number to be used to replace zeros
       If delta is not specified, then the default delta is
       :math:`\delta = \frac{1}{N^2}` where :math:`N`
       is the number of components
    Returns
    -------
    numpy.ndarray, np.float64
       A matrix of proportions where all of the values
       are nonzero and each composition (row) adds up to 1
    Raises
    ------
    ValueError
       Raises an error if negative proportions are created due to a large
       `delta`.
    Notes
    -----
    This method will result in negative proportions if a large delta is chosen.
    References
    ----------
    .. [1] J. A. Martin-Fernandez. "Dealing With Zeros and Missing Values in
           Compositional Data Sets Using Nonparametric Imputation"
    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import multiplicative_replacement
    >>> X = np.array([[.2,.4,.4, 0],[0,.5,.5,0]])
    >>> multiplicative_replacement(X)
    array([[0.1875, 0.375 , 0.375 , 0.0625],
           [0.0625, 0.4375, 0.4375, 0.0625]])
    """

    mat = closure(mat)
    z_mat = mat == 0

    num_feats = mat.shape[-1]
    tot = z_mat.sum(axis=-1, keepdims=True)

    if delta is None:
        delta = (1.0 / num_feats) ** 2

    zcnts = 1 - tot * delta
    if np.any(zcnts) < 0:
        raise ValueError(
            "The multiplicative replacment created negative "
            "proportions. Consider using a smaller `delta`."
        )
    mat = np.where(z_mat, delta, zcnts * mat)
    return mat.squeeze()


def normalize_geometric(adata):
    """Perform geometric normalization on CITEseq data.

    Add description of why geometric normalization

    parameters
    ----------
    adata: :class:`~anndata.AnnData`
        The annotated data matrix.

    returns
    -------
    AnnData
        updates adata.X to contain the geometrically normalized values.

    Example
    -------

    NEED TO IMPLEMENT AN EXAMPLE

    """

    X = adata.X

    X = np.nan_to_num(X)

    # if the matrix is sparse make to array
    if type(X) == csr_matrix:
        X = X.todense()
    # need to add a catch for newly encountered datatype
    elif type(X) == SparseCSRView:
        X = X.todense()

    # ensure that X is an array otherwise this will cause type issue with multiplicative replacement function
    X = np.array(X)

    # replacement of zero values with very small numbers without changing the overall sums
    X = multiplicative_replacement(X)

    # centre log ratio transformation
    X = clr(X)

    adata.X = X

    return None  # adata object is automatically updated
