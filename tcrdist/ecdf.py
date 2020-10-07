import numpy as np
import pandas as pd

__all__ = ['distance_ecdf',
           'make_ecdf_step']

def distance_ecdf(pwrect, thresholds=None, weights=None):
    """Computes the empirical cumulative distribution function (ECDF) for
    each TCR in a set of target TCRs [rows of pwrect] as the proportion
    of reference TCRs [columns of pwrect] within a distance radius less
    than or equal to a threhold d_i, over a range of
    D = [d_1, d_2, ..., d_i]. The distances between pairs of TCRs in the
    target and reference set are contained in the elements of pwrect.

    Optionally, relative weights can be supplied for each reference TCR.
    These can be TCR counts or other weights andthe ECDF will still
    be a probability on [0, 1].

    Parameters
    ----------
    clone_df : pd.DataFrame
    pwrect : np.ndarray or scipy.sparse.csr_matrix, (clone_df.shape[0], n_ref)
    thresholds : np.ndarray
        Vector of thresholds at which the ECDF should be evaluated.
        By default will use all unique values in pwrect.
    weights : np.ndarray or list, (clone_df.shape[0], )
        Relative weight of each TCR in the reference (column dimension of pwrect)

    Returns
    -------
    thresholds : vector, thresholds.shape[0]
    ecdf : vector, thresholds.shape[0]"""
    if weights is None:
        weights = np.ones(pwrect.shape[1])

    if thresholds is None:
        thresholds = np.unique(pwrect[:])

    """Vectorized and faster, using broadcasting for the <= expression"""
    """ecdf = np.sum((pwrect[:, :, None] <= thresholds[None, None, :]) * \
            weights[None, :, None], axis=1) / np.sum(weights)"""
    
    """Decided not to vectorize in 3D because it can create an array that's
    too big for memory"""
    ecdf = np.zeros((pwrect.shape[0], thresholds.shape[0]))
    for i in range(pwrect.shape[0]):
        ecdf[i, :] = np.sum((pwrect[i, :][:, None] <= thresholds[None, :]) * \
                        weights[:, None], axis=0) / np.sum(weights)
    
    return thresholds, ecdf

def make_ecdf_step(thresholds, ecdf, add00=False):
    """Create stepped vector for plotting an ECDF,
    since the ECDF should naturally have discrete steps
    but will not unless they are explictly added prior
    to plotting.

    Takes outputs from distance_ecdf function, but
    is general for creating stepped vectors for plotting.

    Parameters
    ----------
    thresholds : np.ndarray [n_thresholds, ]
        Vector of thresholds
        (typically the x-axis of the plot)
    ecdf : np.ndarray [n_thresholds, ]
        Vector of increasing probabilities
        (typically the y-axis of the plot)
    add00 : bool
        Will add a step to (x=0, y=0) to complete
        the ECDF.

    Returns
    -------
    x, y : np.ndarray
        Vectors for plotting"""
    y = np.asarray(ecdf)
    t = np.asarray(thresholds)
    if add00:
        t = np.concatenate(([0], t.ravel()))
        y = np.concatenate(([0], y.ravel()))

    t = np.concatenate(([t[0]], np.repeat(t[1:].ravel(), 2)))
    y = np.repeat(y.ravel(), 2)[:-1]
    return t, y

def gmean10(vec, axis=0):
    """Geometric mean which may be useful for
    summarizing many ECDF functions"""
    return 10 ** (np.mean(np.log10(vec), axis=axis))