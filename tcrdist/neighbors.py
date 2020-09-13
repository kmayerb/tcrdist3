"""
tcrdist/neigbhors.py

This module contas 
"""

import numpy as np 
import pandas as pd
import parmap

from tcrdist.repertoire import TCRrep
from tcrdist.regex import _index_to_regex_str, _index_to_seqs, _multi_regex 

def compute_ecdf(data, counts=None, thresholds=None):
    """
    Computes the empirical cumulative distribution function at pre-specified
    thresholds. Assumes thresholds is sorted and should be unique.

    data : np.array

    counts : np.array or None 

    threshold : np.array or None
    
    Example
    -------
    >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
    array([0.50083195])
    """
    if thresholds is None:
        thresholds = np.unique(data[:])
    if counts is None:
        counts = np.ones(data.shape)
    
    tot = np.sum(counts)
    # ecdf = np.array([np.sum((data <= t) * counts)/tot for t in thresholds])
    
    """Vectorized and faster, using broadcasting for the <= expression"""
    ecdf = (np.sum((data[:, None] <= thresholds[None, :]) * counts[:, None], axis=0) + 1) / (tot + 1)
    # n_ecdf = (np.sum((data[:, None] <= thresholds[None, :]) * counts[:, None], axis=0) >= n).astype(int)
    return ecdf

  
def _compute_ecdf_rowwise(i, thresholds, data):
    """
    This is a warpper for compute_ecdf using parmap
    """
    ec = compute_ecdf(data = data[i,:], thresholds = thresholds)
    return pd.Series(ec, index = thresholds)


def bkgd_cntl_nn(  tr, 
                   tr_background,
                   ctrl_bkgd = 2*10**-5, 
                   col = 'cdr3_b_aa',
                   add_cols = ['v_b_gene', 'j_b_gene'],
                   ncpus = 4,
                   thresholds = [x for x in range(0,50,2)],
                   generate_regex = True,
                   test_regex = True):   

    """
    bkgd_cntl_nn standfs for background controlled nearest neighbors. 

    tr : tcrdist.repertoire.TCRrep
        TCRrep instance with clone_df of target data
    tr_background : tcrdist.repertoire.TCRrep
        TCRrep instance with clone_df of bulk data
    ctrl_bkgd : float
        Default is 2*10**-5, the acceptable level of background neighbror 
    col : str
        Default is cdr3_b_aa', the column containing CDR3 string
    add_cols  : list
        Extra columns from clone_df to include in centers_df['v_b_gene', 'j_b_gene'],
    ncpus = 4
        passed to pm_processes when using parmap   
    thresholds : list
        Default is [x for x in range(0,50,2)] indicating tcrdist thresholds to compute 
        a ecdf over.
    generate_regex : bool
        if True, generate a regex pattern capturing most if not all of 
        sequences in the neigbhorhood
    test_regex = True
        if True, test regex against both CDR3 in target and background

    Notes
    -----

    The logic behind this function is that we want to find radi that 
    are appropriate for each TCR in a epitope-specific target set 
    that control the frequency of finding neighbors in a background set. 
    For instance, a TCR with common V gene usage and a high probabililty 
    of generation CDR3, will potentially have many nieghbors in a 
    pre-selection background set, whereas a TCR with rare gene usage 
    and a low probability of generation CDR3 can acomodate a larger 
    neighborhood radi. 

    The result dataframe can be sorted to find TCRs that anchor 
    a nieghborhood of TCRs including the most epitope specific
    nieghbors. We hypothesis that these TCRs are more likely to serve as a useful 
    biomarker for searching for biochemically similar TCRs in unlabeled
    datasets. 
    """ 
    thresholds = np.unique(thresholds)
    # for each TCR, we calculate a empirical cummulative 
    # density function along a range of threshold radi
    ecdfs = parmap.map(_compute_ecdf_rowwise, range(0,tr.rw_beta.shape[0]), data = tr.rw_beta, thresholds = thresholds, pm_pbar = True, pm_processes = ncpus)
    # Based on acceptable ctrl_bkgd, we find max acceptable radi from each TCR
    max_radi = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs]
    # THERE IS A NOTABLE BUG IN THE ABOVE LINE!! IF a radius is None (the next line will fail, thus set Nones to 0.
    max_radi = [x if (x is not None) else 0 for x in max_radi]

    # Confirm that this procedure did infact control background hits
    assert np.all([ np.sum(tr.rw_beta[i,:] <= t) < tr.rw_beta.shape[1]*ctrl_bkgd for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))])
    # Tabulate target-set hits. (That is hits within epitope specific set)  
    target_hits = [np.sum(tr.pw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
    background_hits = [np.sum(tr.rw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
    target_neighbors =   [list((tr.pw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
    target_seqs = [tr.clone_df.loc[ar,col].to_list() for ar in target_neighbors]
    background_neighbors = [list((tr.rw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
    background_seqs = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]

    from tcrdist.regex import _index_to_regex_str, _index_to_seqs 

    centers_df = pd.DataFrame({col: tr.clone_df[col].copy(),
                  add_cols[0] : tr.clone_df[add_cols[0]].copy(), 
                  add_cols[1] : tr.clone_df[add_cols[1]].copy(), 
                  'pgen' : tr.clone_df[f'pgen_{col}'],
                  'background_hits': background_hits,
                  'target_hits' : target_hits, 
                  'max_radi': max_radi,
                  'target_neighbors': target_neighbors,
                  'target_seqs': target_seqs,
                  'background_neighbors': background_seqs})

    if generate_regex:
       
        centers_df['regex'] = [_index_to_regex_str(ind = r['target_neighbors'], 
                clone_df = tr.clone_df, 
                pwmat = None, 
                col = col, 
                centroid = tr.clone_df[col][i],
                ntrim = 3,
                ctrim = 2,
                max_ambiguity = 5) for i,r in centers_df.iterrows()]

    if test_regex:
        # Test regex against backgound
        rs = centers_df['regex'].to_list()
        bkgd_cdr3 = tr_background.clone_df[col].to_list()
        
        target_cdr3 = tr.clone_df[col].to_list()
        target_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = target_cdr3 , pm_pbar = True, pm_processes = ncpus)
        target_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in target_regex_hits ]
        centers_df['target_re_hits'] = target_re_hits

        bkgd_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = bkgd_cdr3, pm_pbar = True, pm_processes = ncpus)
        bkgd_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in bkgd_regex_hits]
        centers_df['bkgd_re_hits'] = bkgd_re_hits

    return centers_df
