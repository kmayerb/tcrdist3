import numpy as np
import pandas as pd
from scipy import sparse
from matplotlib import pyplot as plt
import matplotlib as mpl

__all__ = ['distance_ecdf',
           'make_ecdf_step',
           'plot_ecdf']

def distance_ecdf(pwrect, thresholds=None, weights=None, pseudo_count=0, skip_diag=False, absolute_weight = False):
    """Computes the empirical cumulative distribution function (ECDF) for
    each TCR in a set of target TCRs [rows of pwrect] as the proportion
    of reference TCRs [columns of pwrect] within a distance radius less
    than or equal to a threhold d_i, over a range of
    D = [d_1, d_2, ..., d_i]. The distances between pairs of TCRs in the
    target and reference set are contained in the elements of pwrect.

    Optionally, relative weights can be supplied for each reference TCR.
    These can be TCR counts or other weights and the ECDF will still
    be a probability on [0, 1].

    Parameters
    ----------
    pwrect : np.ndarray or scipy.sparse.csr_matrix, (clone_df.shape[0], n_ref)
        Matrix of pairwise distances among TCRs (can be rectangular). Sparse matrix may encode
        true zeros (i.e. identical TCRs) as -1 or 0, without having effect on the ECDF.
    thresholds : np.ndarray
        Vector of thresholds at which the ECDF should be evaluated.
        By default will use all unique values in pwrect. For sparse arrays the thresholds
        should not exceed the largest distance encoded.
    weights : np.ndarray or list, (clone_df.shape[0], )
        Relative weight of each TCR in the reference (column dimension of pwrect)
    pseudo_count : int
        Added to the numerator and denominator at each threshold
        to avoid zero. Useful if end goal is a log-scale plot.
    skip_diag : bool
        Skip counting the diagonal for computing ECDF of seqs against same seqs.
    absolute_weight : bool
        if True, denominator is number of total sequences in pwrect.shape[2] rather than sum of the weights.
        
    Returns
    -------
    thresholds : vector, thresholds.shape[0]
    ecdf : vector, thresholds.shape[0]"""
    if weights is None:
        weights = np.ones(pwrect.shape[1])
    else:
        weights = np.asarray(weights)

    if thresholds is None:
        if sparse.issparse(pwrect):
            thresholds = np.unique(pwrect.data)
        else:
            thresholds = np.unique(pwrect[:])
    else:
        thresholds = np.asarray(thresholds)

    """Vectorized and faster, using broadcasting for the <= expression"""
    """ecdf = np.sum((pwrect[:, :, None] <= thresholds[None, None, :]) * \
            weights[None, :, None], axis=1) / np.sum(weights)"""
    
    """Decided not to vectorize in 3D because it can create an array that's
    too big for memory"""
    ecdf = np.zeros((pwrect.shape[0], thresholds.shape[0]))
    sum_weights = np.sum(weights)
    if absolute_weight:
        denom = pwrect.shape[1]

    for i in range(pwrect.shape[0]):
        if sparse.issparse(pwrect):
            """By using the row.data attribute we only sum over distances
            that were stored: large distances that were discarded from the sparse
            representation (which are sometimes thought of as zeros) are simply
            not present and therefore never counted as below threshold[i]"""
            row = pwrect[i, :]
            numer = np.sum((row.data[:, None] <= thresholds[None, :]) * weights[row.indices, None], axis=0)
        else:
            row = np.reshape(pwrect[i, :], (pwrect.shape[1], 1))
            numer = np.sum((row <= thresholds[None, :]) * weights[:, None], axis=0)
        if absolute_weight:
            denom = pwrect.shape[1]
        else:
            denom = sum_weights
        if skip_diag:
            numer = numer - weights[i]
            denom = denom - weights[i]
        ecdf[i, :] = (numer + pseudo_count) / (denom + pseudo_count)
    return thresholds, ecdf

def make_ecdf_step(thresholds, ecdf, add_mnx=False, add_mny=False, add_mnmn=False,enforce_mn=False , mn=(0, 0), xjitter=0):
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
    add_mnmn : bool
        Will add a step to (x=mn[0], y=mn[1]) to complete
        the ECDF.

    Returns
    -------
    x, y : np.ndarray
        Vectors for plotting"""
    y = np.asarray(ecdf)
    t = np.asarray(thresholds)

    if enforce_mn:
        t[t<mn[0]] = mn[0]
        y[y<mn[1]] = mn[1]

    if add_mnx:
        t = np.concatenate(([mn[0]], t.ravel()))
        y = np.concatenate(([y[0]], y.ravel()))
    if add_mny:
        t = np.concatenate(([t[0]], t.ravel()))
        y = np.concatenate(([mn[1]], y.ravel()))

    if add_mnmn:
        t = np.concatenate(([mn[0]], t.ravel()))
        y = np.concatenate(([mn[1]], y.ravel()))

    t = np.concatenate(([t[0]], np.repeat(t[1:].ravel(), 2)))
    y = np.repeat(y.ravel(), 2)[:-1]

    jx = t + (np.random.rand(1) - 0.5) * xjitter
    jx[0] = t[0]
    jx[-1] = t[-1]
    return jx, y

def gmean10(vec, axis=0):
    """Geometric mean which may be useful for
    summarizing many ECDF functions"""
    return 10 ** (np.mean(np.log10(vec), axis=axis))


def plot_ecdf(
    thresholds,
    ecdf_mat, 
    ax = None,
    ylim = None,
    xlabel = f'Distance From Target TCR Clone',
    ylabel = f'Proportion of Reference TCRs', 
    plot_mean = True,
    min_freq = 1E-10):
    """
    A very basic ecdf plot

    Parameters
    ----------
    thresholds : np.array
        1D array tcrdistance threshold trdistance units (tdus)
    ecdf_mat : np.array
        2D array, each row repressents a TCR and the proportion of neighbors with tdus
    ax : matplotlib axes or None

    ylabel : str
        'Proportion of X TCRs'
    xlabel : str
        'Proportion of Reference TCRs'
    plot_mean : bool
        If true, plot the mean value
    """
    if ax is None:
        ax = plt.gca()
    if not ylim is None:
        ax.set_ylim(ylim)

    ax.set_yscale('log')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    
    for row_i in range(ecdf_mat.shape[0]):
        make_ecdf_step(thresholds, ecdf_mat[row_i, :], add_mnx=True, enforce_mn=True, mn=(0, min_freq), xjitter=1)
        x, y = make_ecdf_step(thresholds, ecdf_mat[row_i, :])
        ax.plot(x, y, color='k', alpha=0.2)
    if plot_mean:
        x, y = make_ecdf_step(thresholds, np.mean(ecdf_mat, axis=0))
        ax.plot(x, y, color='r', alpha=1)
    return True



def _plot_manuscript_ecdfs(
    thresholds, 
    ecdf_mat, 
    ylab = 'Proportion of X TCRs', 
    cdr3_len=None, 
    min_freq=1e-6,
    cdr3_len_min=10., 
    cdr3_len_max=16.,
    low_pass = -5,
    cmap=mpl.cm.viridis_r):
    """
    _plot_manuscript_ecdfs, recreate the manuscript type ecdf

    Parameters
    ----------
    thresholds : np.array
        1D array tcrdistance threshold trdistance units (tdus)
    ecdf_mat : np.array
        2D array, each row repressents a TCR and the proportion of neighbors with tdus
    ylab : str
        'Proportion of X TCRs'
    cdr3_len : list or pd.Series
        Lengths of cdr3, must match order and dimension of rows in ecdf_mat
    min_freq : float
        Minimum frequency or proportion to display
    cdr3_len_min : float 

    cdr3_len_max : float
    
    low_pass : int
           controls threshold frequencyh to show on lower gray bar (-5) corresponds with 10^-5

    """
    # Define colormap range for CDR3 length
    low_pass_num = 10**low_pass # -5 -> 10^-5
    norm = mpl.colors.Normalize(vmin=cdr3_len_min, vmax=cdr3_len_max)
    if cdr3_len is None:
        cdr3_len = 10 * np.ones(ecdf_mat.shape[0])

    # Declare Figure <figh>
    figh = plt.figure(figsize=(11, 8))
    # Declare Grid
    gs = plt.GridSpec(nrows=2, ncols=2, width_ratios=[10, 1], height_ratios=[7, 1], hspace=0.15)
        
    # [0,1] : Assign color bar axis
    cbar = mpl.colorbar.ColorbarBase(figh.add_subplot(gs[0, 1]), cmap=cmap,
                            norm=norm,
                            orientation='vertical')
    cbar.set_label('CDR3 length')

    # [0,0] : Asign ECDF axis
    axh = figh.add_subplot(gs[0, 0], yscale='log')
    axh.set_ylabel(ylab)
    # plt.xlabel(f'Neighborhood radius\n(tcrdist units)')
    for tari in range(ecdf_mat.shape[0]):
        x, y = make_ecdf_step(thresholds, ecdf_mat[tari, :], add_mnx=True, enforce_mn=True, mn=(0, min_freq), xjitter=1)
        axh.plot(x, y, color=mpl.cm.viridis_r(norm(cdr3_len[tari]))[:3], alpha=0.4)
    plt.annotate(text=f'n={ecdf_mat.shape[0]}',
                 xy=(min_freq, 0.5),
                 xytext=(3, -3),
                 textcoords='offset points',
                 ha='left',
                 va='top')
    axh.set_xticklabels([])
    axh.set_ylim((min_freq, 0.5))
    axh.set_xlim((0, thresholds[-1]))

    # [1, 0] : Assign underplot
    axh2 = figh.add_subplot(gs[1, 0])
    x, y = make_ecdf_step(thresholds, np.mean(ecdf_mat<low_pass_num, axis=0), add_mnx=False)
    axh2.fill_between(x, np.zeros(y.shape[0]), 100*y, color='gray')
    axh2.set_ylim((0, 100))
    axh2.set_ylabel(f'% TCRs\n$ECDF < 10^{low_pass}$')
    axh2.set_xlabel(f'Neighborhood radius\n(tcrdist units)')
    axh2.set_xlim((0, thresholds[-1]))
    return figh
