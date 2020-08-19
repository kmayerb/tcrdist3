import numpy as np
import pandas as pd
from scipy import stats

__all__ = ['generalized_simpsons_entropy',
           'simpsons_difference',
           'effective_number',
           'fuzzy_diversity']

def _calc_Z_base(counts, r_vec):
    """Python function for computing Z that is compatible with
    numba compilation."""
    n = np.sum(counts)
    
    K = len(counts)
    p = counts / n
    
    Z_out = np.zeros(len(r_vec))
    for r_i in prange(len(r_vec)):
        r = r_vec[r_i]
        Z = 0
        for i in range(K):
            if counts[i] == 0:
                continue
            else:
                prod = 1
                for k in range(1, r):
                    prod *= 1 - (counts[i] - 1) / (n - k)
                Z += prod * p[i]
        Z_out[r_i] = Z
    return Z_out

def _calc_stdev_base(counts, r_vec):
    """Python function for computing Z standard deviation that is
    compatible with numba compilation."""
    n = np.sum(counts) # samples1
    p = counts / n
    K = len(p)

    stdev = np.zeros(len(r_vec))
    for r_i in prange(len(r_vec)):
        r = r_vec[r_i]
        """Assumes that the last count is not 0"""
        # index = c1.shape[0] - 1
        h_hat = np.zeros(K - 1)
        tmp = (1 - p[:-1])**r + r*p[:-1]*(1 - p[:-1])**(r-1) - (1 - p[-1])**r - r*p[-1]*(1 - p[-1])**(r-1)
        h_hat[p[:-1] > 0] = tmp[p[:-1] > 0]

        sigma = np.zeros((K-1, K-1))
        for i in range(K-1):
            for j in range(K-1):
                if i == j:
                    sigma[i, j] = p[i] * (1 - p[i])
                else:
                    sigma[i, j] = -p[i] * p[j]
        
        v = np.dot(np.dot(h_hat.T, sigma), h_hat)
        """
        v = 0
        for i in range(K-1):
            for j in range(K-1):
                v += h_hat[i] * sigma[i, j] * h_hat[j]
        """
        stdev[r_i] = np.sqrt(v)
    return stdev

try:
    """Try to import and compile using jit, otherwise fall back on
    numpy and python loops (slow for large datasets)"""
    from numba import jit, prange

    _calc_Z = jit(_calc_Z_base, nopython=True, error_model='numpy', parallel=True)
    _calc_stdev = jit(_calc_stdev_base, nopython=True, error_model='numpy', parallel=True)
except ImportError:
    prange = range
    _calc_Z = _calc_Z_base
    _calc_stdev = _calc_stdev_base

def simpsons_difference(counts1, counts2, orders=[2], alpha=0.05):
    """Difference in diversity between two communities.

    Parameters
    ----------
    counts1, counts2 :  : np.ndarray or pd.Series
        Vector of counts for each species.
    orders : np.ndarray or list of integers
        Order for calculation. r = 2 is equivalent to common Simpson's entropy computations.
        Increasing r gives more relative importance to rare species.
    alpha : float
        Upper and lower confidence levels define the 1 - alpha/2 confidence interval

    Returns
    -------
    difference, lcl, ucl : np.ndarray, shape len(orders)"""
    orders = np.asarray(orders).astype(float)

    counts1 = np.asarray(counts1).astype(float)
    n1 = np.sum(counts1)

    counts2 = np.asarray(counts2).astype(float)
    n2 = np.sum(counts2)

    Z1 = _calc_Z(counts1, orders)
    sdev1 = _calc_stdev(counts1, orders)

    Z2 = _calc_Z(counts2, orders)
    sdev2 = _calc_stdev(counts2, orders)

    d = Z1 - Z2

    criticalz = -stats.norm.ppf(alpha / 2)
    
    lcl = d - criticalz * np.sqrt(sdev1**2 / n1 + sdev2**2 / n2)
    ucl = d + criticalz * np.sqrt(sdev1**2 / n1 + sdev2**2 / n2)
    return d, lcl, ucl

def effective_number(Z, orders=[2]):
    """Effective number of species is the number of equiprobable
    species that would yield the same diversity as a given distribution.

    As it is a monotonic transformation this function can be used to
    transform confidence intervals on Z as well.

    Parameters
    ----------
    Z : np.ndarray of floats, [0, 1]
        Generalized Simpson's entropies
    orders : np.ndarray of integers
        Order of the generalized Simpson's entropy.
        Must match the orders used in the calculation of Z's.

    Returns
    -------
    D : float
        Effective number"""
    orders = np.asarray(orders).astype(float)
    return 1 / (1 - Z**(1 / orders))

def generalized_simpsons_entropy(counts, orders=[2], alpha=0.05):
    """Generalized Simpson’s entropy of order r can be interpreted as
    the average information brought by the observation of an individual/species.
    
    Its information function I(p) = (1 − p)*r represents the probability of
    not observing a single individual of a species with proportion p in a sample
    of size r. Thus I is an intuitive measure of rarity.

    Above is quoted from:
    Grabchak M, Marcon E, Lang G, Zhang Z (2017) The generalized Simpson’s entropy
        is a measure of biodiversity. PLoS ONE 12(3): e0173305.
        https://doi.org/10.1371/journal.pone.0173305

    It is common to evaluate using a series of orders as long as r < len(counts) - 1

    Parameters
    ----------
    counts : np.ndarray or pd.Series
        Vector of counts for each species.
    orders : np.ndarray of integers
        Order for calculation. r = 2 is equivalent to common Simpson's entropy computations.
        Increasing r gives more relative importance to rare species.
    alpha : float
        Upper and lower confidence levels define the 1 - alpha/2 confidence interval

    Returns
    -------
    Z, LCL, UCL : floats or as pd.Series if counts is a pd.Series
        Generalized Simpson's entropy Z and lower (upper) confidence limits."""
    orders = np.asarray(orders).astype(float)

    if type(counts) is pd.Series:
        return_series = True
        name = counts.name
        counts = counts.values.astype(float)
    else:
        counts = np.asarray(counts).astype(float)
    n = np.sum(counts)

    Z = _calc_Z(counts, orders)
    sdev = _calc_stdev(counts, orders)
    criticalz = -stats.norm.ppf(alpha / 2)

    lcl = Z - criticalz * sdev / np.sqrt(n)
    ucl = Z + criticalz * sdev / np.sqrt(n)

    if return_series:
        return pd.DataFrame({'order':orders.astype(int), 'Z':Z, 'Z_LCL':lcl, 'Z_UCL':ucl,
                              'D':effective_number(Z, orders),
                              'D_LCL':effective_number(lcl, orders),
                              'D_UCL':effective_number(ucl, orders)}).set_index('order')
    else:
        return Z, lcl, ucl

def _non_generalized_simpsons_index(counts):
    p = counts / np.sum(counts)
    D = (p * p).sum()
    return 1 - D


def fuzzy_diversity(counts, pwmat, order=2, threshold=1, nsamples=1000, force_sampling=False):
    """Compute a "fuzzy" diversity index that takes into account the number of similar
    members in a community, not just the ones that are identically matched.

    Idea published as "TCRdiv" in Dash et al. (2017) in the area of T cell receptor
    (TCR) analysis. The community is a TCR repertoire and the members are unique
    TCR clones.

    The diversity index is similar to Simpson's diversity index (SDI) in that it estimates
    the chance of sampling m members from the community that are all more similar
    to each other than some pre-specified threshold.

    Fuzzy diversity is the sum of the probability of n members that match exactly,
    plus the probability of m members matching within the threshold. Or it can be expressed
    as 1 - probability of drawing m members that are not similar (distance > threshold).

    For a confidence interval use some kind of bootstrap approach.

    Parameters
    ----------
    counts : np.ndarray [n]
        Number of observations of each member in the community
    pwmat : np.ndarray [n x n]
        Pairwise distances (not similarities) for all members in the community.
        Consider using 1 to represent members that are not similar, and 0 for others,
        or use the threshold parameter.
    order : int
        Order of the SDI, i.e. the number of members sampled to see if they are
        similar to one another
    threshold : float
        A distance threshold that defines whether two members are similar to one
        another: similar if distance < threshold
    force_sampling : bool
        Force use of the sampling approach, even for order=2
        Only useful for checking the order=2 result.

    Returns
    -------
    div : float
        Estimate of fuzzy diversity"""

    counts = np.asarray(counts)[:, None]

    n = np.sum(counts)
    
    cts = counts * counts.T

    assert cts.shape[0] == pwmat.shape[0]
    assert cts.shape[0] == pwmat.shape[1]

    up = np.triu_indices_from(cts, k=1)

    if order == 2 and not force_sampling:
        """The pwmat represents drawing two members and seeing if they are similar
        so order = 2 is just a simple summary of the pairwise distance matrix:
        the proportion of distance pairs that represent similar members"""
        div = np.sum(cts[up] * (pwmat < threshold)[up]) / np.sum(cts)
    else:
        """Higher orders could be computed using a sampling approach
        since its not scalable to compute all possible triplets, etc.
        to see which are all similar. A triplet consists of 3 pairwise distances
        that must all be within the threshold to be a similar triplet."""
        dvec = (pwmat < threshold)[up]
        prob = cts[up] / np.sum(cts[up])
        samples = np.random.choice(dvec, size=(order, nsamples), replace=True, p=prob)
        
        div = np.mean(np.all(samples, axis=0))

    return 1 - div





