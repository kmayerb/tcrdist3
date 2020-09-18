"""
tcrdist/neigbhors.py

This module contas 
"""

import numpy as np 
import pandas as pd
import parmap

from tcrdist.repertoire import TCRrep
from tcrdist.regex import _index_to_regex_str, _index_to_seqs, _multi_regex 
from scipy.stats import chi2_contingency


def compute_population_estimate_ecdf(data, weights = None, thresholds=None):
    """
    """
    if thresholds is None:
        thresholds = np.unique(data[:])
  
    if weights is None:
        weights = np.ones(data.shape)
    else:
        weights = np.array(weights)
    
    """Vectorized and faster, using broadcasting for the <= expression"""
    # Step 1 : compare a set of n distances to a set of say 25 threshold 
        # (n dists x 1) broadcast against (1 x 25 thresholds) ==> (n, 25), where each column representd number of hits below the threshold 
    M= 1*(data[:, None] <= thresholds[None, :])
    # Step 2 : Adjust hit be relative weight (i.e. actual_freq_vj / sampled_freq_vj)
        # Broadcast (n, 25) on weights (n, 1) ==> (n, 25)
    M = M * weights[:,None]
    # Step 3: Take row sums, producing the weighted sum at each threshold ==> (25,)
    M =  np.sum(M, axis = 0)
    # Step 4: Divide by the weighted total (i.e., the amount if all seqs were present, which because of enrichment is far greater than 1)
    print(M)   
    ecdf = M / len(weights) 
    return ecdf



def compute_ecdf(data, weights = None, counts=None, thresholds=None):
    """
    Computes the empirical cumulative distribution function at pre-specified
    thresholds. Assumes thresholds is sorted and should be unique.

    data : np.array

    threshold : np.array or None
    
    counts : np.array or None

    weights :  np.array or None
         vector of W_vj  (i.e  actual_freq_vj / sampled_freq_vj)
    
    Example
    -------
    >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
    array([0.5])
    >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), weights =np.array([1,1,1,1,1,1]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
    array([0.5])
    >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), weights =np.array([10,1,1,1,1,1]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
    array([0.8])

    Notes
    -----
    
    Weights are used when we enriche for V,J combinations that are rare. 
    For instance if some V,J pairing 
    occurs at an actual frequency of 1/100 but we provide it at 1/10, 
    we provide an enrichment weight 
    (actual_freq_vj / sampled_freq_v)


    -----
    """

    if thresholds is None:
        thresholds = np.unique(data[:])
    if counts is None:
        counts = np.ones(data.shape)
    
    if weights is None:
        weights = np.ones(data.shape)
    else:
        weights = np.array(weights)
    
    tot = np.sum(counts * weights) 
    """Vectorized and faster, using broadcasting for the <= expression"""
    # Step 1 : compare a set of n distances to a set of say 25 threshold 
        # (n dists x 1) broadcast against (1 x 25 thresholds) ==> (n, 25), where each column representd number of hits below the threshold 
    M= 1*(data[:, None] <= thresholds[None, :])
    # Step 2 : Adjust hit be relative weight (i.e. actual_freq_vj / sampled_freq_vj)
        # Broadcast (n, 25) on weights (n, 1) ==> (n, 25)
    M = M * weights[:,None] * counts[:,None]
    # Step 3: Take row sums, producing the weighted sum at each threshold ==> (25,)
    M =  np.sum(M, axis = 0)
    # Step 4: Divide by the weighted total (i.e., the amount if all seqs were present, which because of enrichment is far greater than 1)
    ecdf = M / tot

    return ecdf


def _compute_pop_estimate_ecdf_rowwise(i, thresholds, data):
    """
    This is a warpper for compute_ecdf using parmap
    """
    ec = compute_population_estimate_ecdf(data = data[i,:], thresholds = thresholds)
    return pd.Series(ec, index = thresholds)

def _compute_weighted_pop_estimate_ecdf_rowwise(i, thresholds, weights, data):
    """
    This is a warpper for compute_ecdf using parmap
    """
    ec = compute_population_estimate_ecdf(data = data[i,:], weights = weights, thresholds = thresholds)
    return pd.Series(ec, index = thresholds)



def _compute_ecdf_rowwise(i, thresholds, data):
    """
    This is a warpper for compute_ecdf using parmap
    """
    ec = compute_ecdf(data = data[i,:], thresholds = thresholds)
    return pd.Series(ec, index = thresholds)

def _compute_weighted_ecdf_rowwise(i, thresholds, weights, data):
    """
    This is a warpper for compute_ecdf using parmap
    """
    ec = compute_ecdf(data = data[i,:], weights = weights, thresholds = thresholds)
    return pd.Series(ec, index = thresholds)


def compute_odds_ratio(pos, neg, bpos, bneg, ps = 1):
    """
    pos : int or float
        number of positve instances in sample
    neg : int or float
        number of negative instances in sample
    bpos : int or float 
        number of positve instances in background population
    bneg : int of float 
        number of negative instances in background population
    ps : int or float
        psuedocount to avoid infinite odds ratio when zeros are possible in the denominator

    Notes
    -----
    An odds ratio (OR) is a statistic that quantifies the strength of the association between two events, A and B. The odds ratio is defined as the ratio of the odds of A in the presence of B and the odds of A in the absence of B

    Examples
    --------    
    >>> compute_odds_ratio(10, 100, 10, 1000, ps = 0)
    10.0
    >>> compute_odds_ratio(10, 100, 10, 1000, ps = 1)
    9.91089108910891
    """
    odds_ratio = ((ps+pos)/(ps+neg))/ ((ps+bpos)/(ps+bneg))
    return odds_ratio

def compute_rate(pos,neg, ps = 1):
    return ((pos + ps) / (pos + neg + 2*ps))

def bkgd_cntl_nn2( tr, 
                   tr_background,
                   weights = None,
                   ctrl_bkgd = 2*10**-5, 
                   col = 'cdr3_b_aa',
                   add_cols = ['v_b_gene', 'j_b_gene'],
                   ncpus = 4,
                   include_seq_info= True,
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
    if weights is None:
        weights = np.one(1, tr.rw_beta.shape[1])
    thresholds = np.unique(thresholds)
    # for each TCR, we calculate a empirical cummulative 
    # density function along a range of threshold radi
    ecdfs = parmap.map(_compute_weighted_ecdf_rowwise, range(0,tr.rw_beta.shape[0]), 
        data = tr.rw_beta, 
        weights = weights,
        thresholds = thresholds, 
        pm_pbar = True, 
        pm_processes = ncpus)

    # Based on acceptable ctrl_bkgd, we find max acceptable radi from each TCR
    max_radi = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs]
    # THERE IS A NOTABLE BUG IN THE ABOVE LINE!! IF a radius is None (the next line will fail, thus set Nones to 0.
    max_radi = [x if (x is not None) else 0 for x in max_radi]

    target_hits    = [np.sum(tr.pw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
    bkgd_hits      = [np.sum(tr.rw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
    bkgd_hits_weighted = [np.sum(1*(tr.rw_beta[i,:] <= t) * (weights)) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
 
    bkgd_total = tr.rw_beta.shape[1]
    bkgd_weighted_total = tr.rw_beta.shape[1] #np.sum(weights)
    ctrl = np.array(bkgd_hits) / bkgd_total 
    ctrl_weighted= np.array(bkgd_hits_weighted)/bkgd_weighted_total 

    centers_df = pd.DataFrame({col: tr.clone_df[col].copy(),
              add_cols[0] : tr.clone_df[add_cols[0]].copy(), 
              add_cols[1] : tr.clone_df[add_cols[1]].copy(), 
              'pgen' : tr.clone_df[f'pgen_{col}'],
              'max_radi': max_radi,
              'target_hits' : target_hits, 
              'bkgd_hits': bkgd_hits,
              'bkgd_hits_weighted': bkgd_hits_weighted,
              'bkgd_total' :  bkgd_total,
              'bkgd_weighted_total' : bkgd_weighted_total,
              'ctrl' : ctrl,
              'ctrl_weighted' : ctrl_weighted})
   

    # Total bakground counts
    n1 = tr.clone_df.shape[0]   
    n2 = bkgd_weighted_total 
    print(f"N1 {n1}")
    print(f"N2 {n2}")  
    centers_df['target_misses'] = centers_df['target_hits'].apply(lambda x : n1-x)
    centers_df['TR'] = [compute_rate(pos=r['target_hits'], neg=(n1-r['target_hits'])) for i,r in centers_df.iterrows()]
    centers_df['TR2'] = [compute_rate(pos=r['target_hits'], neg=r['target_misses']) for i,r in centers_df.iterrows()]
    centers_df['BR_weighted'] = [compute_rate(pos=r['bkgd_hits_weighted'], 
                                    neg=n2-r['bkgd_hits_weighted']) for i,r in centers_df.iterrows()]
    centers_df['RR_weighted'] = centers_df['TR']/centers_df['BR_weighted']
    centers_df['OR_weighted'] =[compute_odds_ratio(pos=r['target_hits'], 
                                       neg=n1-r['target_hits'], 
                                       bpos=r['bkgd_hits_weighted'], 
                                       bneg= n2-r['bkgd_hits_weighted'], ps = 1) for i,r in centers_df.iterrows()]
    centers_df['chi2dist'] = [chi2_contingency(np.array(
            [[1+r['target_hits'], 1+n1-r['target_hits']],[1+r['bkgd_hits_weighted'], 1+n2-r['bkgd_hits_weighted']]]))[0] for _,r in centers_df.iterrows() ]

    if include_seq_info:
        
        target_neighbors = [list((tr.pw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
        centers_df['target_neighbors'] = target_neighbors 
        centers_df['target_seqs'] = [tr.clone_df.loc[ar,col].to_list() for ar in target_neighbors]
        background_neighbors = [list((tr.rw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
        centers_df['background_neighbors'] = background_neighbors
        col = "cdr3_b_aa"
        centers_df['background_seqs'] = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]
        col = "v_b_gene" 
        centers_df['background_v'] = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]
        col = "j_b_gene" 
        centers_df['background_j'] = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]
        col  = 'adj_freq_pVJ'
        centers_df['adj_freq'] = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]
    
    if generate_regex:
        col = "cdr3_b_aa"
        centers_df['regex'] = [_index_to_regex_str(ind = r['target_neighbors'], 
                clone_df = tr.clone_df, 
                pwmat = None, 
                col = col, 
                centroid = tr.clone_df[col][i],
                ntrim = 3,
                ctrim = 2,
                max_ambiguity = 5) for i,r in centers_df.iterrows()]

    if test_regex:
        # THIS INVOLVES TESTING REGEX AGAINST TARGET AND BACKGROUND
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
        bkgd_weighted_re_hits = [np.sum(np.array([1 if (x is not None) else 0 for x in l]) * (weights)) for l in bkgd_regex_hits] 
        # Compute Relative Rates
        centers_df['bkgd_re_weighted_hits'] = bkgd_weighted_re_hits 
        centers_df['TR_re'] = [compute_rate(pos=r['target_re_hits'], neg=n1-r['target_re_hits']) for i,r in centers_df.iterrows()]
        centers_df['BR_re_weighted'] = [compute_rate(pos=r['bkgd_re_weighted_hits'], 
                                     neg=n2-r['bkgd_re_weighted_hits']) for i,r in centers_df.iterrows()]
        centers_df['RR_re_weighted'] = centers_df['TR']/centers_df['BR_re_weighted']
        centers_df['OR_re_weighted'] =[compute_odds_ratio(pos=r['target_re_hits'], 
                                       neg=n1-r['target_re_hits'], 
                                       bpos=r['bkgd_re_weighted_hits'], 
                                       bneg= n2-r['bkgd_re_weighted_hits'], ps = 1) for i,r in centers_df.iterrows()]
        centers_df['chi2re'] = [chi2_contingency(np.array(
            [[1+r['target_re_hits'], 1+n1-r['target_re_hits']],[1+r['bkgd_re_weighted_hits'], 1+n2-r['bkgd_re_weighted_hits']]]))[0] for _,r in centers_df.iterrows() ]
        beta_re = 1
        beta_dist = 1
        centers_df['chi2joint'] = [beta_re  * r['chi2re'] + beta_dist* r['chi2dist'] for _,r in centers_df.iterrows() ]
#
# 
    return centers_df



              












# def bkgd_cntl_nn(  tr, 
#                    tr_background,
#                    weights = None,
#                    ctrl_bkgd = 2*10**-5, 
#                    col = 'cdr3_b_aa',
#                    add_cols = ['v_b_gene', 'j_b_gene'],
#                    ncpus = 4,
#                    thresholds = [x for x in range(0,50,2)],
#                    generate_regex = True,
#                    test_regex = True):   

#     """
#     bkgd_cntl_nn standfs for background controlled nearest neighbors. 

#     tr : tcrdist.repertoire.TCRrep
#         TCRrep instance with clone_df of target data
#     tr_background : tcrdist.repertoire.TCRrep
#         TCRrep instance with clone_df of bulk data
#     ctrl_bkgd : float
#         Default is 2*10**-5, the acceptable level of background neighbror 
#     col : str
#         Default is cdr3_b_aa', the column containing CDR3 string
#     add_cols  : list
#         Extra columns from clone_df to include in centers_df['v_b_gene', 'j_b_gene'],
#     ncpus = 4
#         passed to pm_processes when using parmap   
#     thresholds : list
#         Default is [x for x in range(0,50,2)] indicating tcrdist thresholds to compute 
#         a ecdf over.
#     generate_regex : bool
#         if True, generate a regex pattern capturing most if not all of 
#         sequences in the neigbhorhood
#     test_regex = True
#         if True, test regex against both CDR3 in target and background

#     Notes
#     -----

#     The logic behind this function is that we want to find radi that 
#     are appropriate for each TCR in a epitope-specific target set 
#     that control the frequency of finding neighbors in a background set. 
#     For instance, a TCR with common V gene usage and a high probabililty 
#     of generation CDR3, will potentially have many nieghbors in a 
#     pre-selection background set, whereas a TCR with rare gene usage 
#     and a low probability of generation CDR3 can acomodate a larger 
#     neighborhood radi. 

#     The result dataframe can be sorted to find TCRs that anchor 
#     a nieghborhood of TCRs including the most epitope specific
#     nieghbors. We hypothesis that these TCRs are more likely to serve as a useful 
#     biomarker for searching for biochemically similar TCRs in unlabeled
#     datasets. 
#     """ 
#     thresholds = np.unique(thresholds)
#     # for each TCR, we calculate a empirical cummulative 
#     # density function along a range of threshold radi
#     if weights is None:
#         ecdfs = parmap.map(_compute_ecdf_rowwise, range(0,tr.rw_beta.shape[0]), 
#             data = tr.rw_beta, 
#             thresholds = thresholds, 
#             pm_pbar = True, 
#             pm_processes = ncpus)
#     else: 
#         ecdfs = parmap.map(_compute_weighted_ecdf_rowwise, range(0,tr.rw_beta.shape[0]), 
#             data = tr.rw_beta, 
#             weights = weights,
#             thresholds = thresholds, 
#             pm_pbar = True, 
#             pm_processes = ncpus)

#     # Based on acceptable ctrl_bkgd, we find max acceptable radi from each TCR
#     max_radi = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs]
#     # THERE IS A NOTABLE BUG IN THE ABOVE LINE!! IF a radius is None (the next line will fail, thus set Nones to 0.
#     max_radi = [x if (x is not None) else 0 for x in max_radi]

#     # Confirm that this procedure did infact control background hits
#     # NOTE WE CAN NO LONGER GAURANTEE THIS IS TRUE IF NUMBER OF PERFECT MATCHES EXCEEDS 0
#     # assert np.all([ np.sum(tr.rw_beta[i,:] <= t) < tr.rw_beta.shape[1]*ctrl_bkgd for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))])
#     # Tabulate target-set hits. (That is hits within epitope specific set)  
    
    
#     if weights is None:
#         target_hits = [np.sum(tr.pw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
#         background_hits = [np.sum(tr.rw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
#     else:
#         target_hits              = [np.sum(tr.pw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
#         target_weighted_hits     = [np.sum(tr.pw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
#         background_hits          = [np.sum(tr.rw_beta[i,:] <= t) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
#         background_weighted_hits = [np.sum(1*(tr.rw_beta[i,:] <= t) * (1/weights)) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
    
#     target_neighbors =   [list((tr.pw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
#     target_seqs = [tr.clone_df.loc[ar,col].to_list() for ar in target_neighbors]
#     background_neighbors = [list((tr.rw_beta[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,tr.pw_beta.shape[0]))]
#     background_seqs = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]

    

#     centers_df = pd.DataFrame({col: tr.clone_df[col].copy(),
#                   add_cols[0] : tr.clone_df[add_cols[0]].copy(), 
#                   add_cols[1] : tr.clone_df[add_cols[1]].copy(), 
#                   'pgen' : tr.clone_df[f'pgen_{col}'],
#                   'bkgd_hits': background_hits,
#                   'background_hits': background_hits, # redundant so we don't break anything
#                   'target_hits' : target_hits, 
#                   'max_radi': max_radi,
#                   'target_neighbors': target_neighbors,
#                   'target_seqs': target_seqs,
#                   'background_neighbors': background_seqs})
    
#     if weights is not None:
#         centers_df['target_weighted_hits'] = target_weighted_hits
#         centers_df['bkgd_weighted_hits'] = background_weighted_hits


#     if generate_regex:
       
#         centers_df['regex'] = [_index_to_regex_str(ind = r['target_neighbors'], 
#                 clone_df = tr.clone_df, 
#                 pwmat = None, 
#                 col = col, 
#                 centroid = tr.clone_df[col][i],
#                 ntrim = 3,
#                 ctrim = 2,
#                 max_ambiguity = 5) for i,r in centers_df.iterrows()]

#     if test_regex:
#         from tcrdist.regex import _index_to_regex_str, _index_to_seqs 
#         # Test regex against backgound
#         rs = centers_df['regex'].to_list()
#         bkgd_cdr3 = tr_background.clone_df[col].to_list()
        
#         target_cdr3 = tr.clone_df[col].to_list()
#         target_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = target_cdr3 , pm_pbar = True, pm_processes = ncpus)
#         target_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in target_regex_hits ]
#         centers_df['target_re_hits'] = target_re_hits

#         bkgd_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = bkgd_cdr3, pm_pbar = True, pm_processes = ncpus)
#         bkgd_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in bkgd_regex_hits]
#         centers_df['bkgd_re_hits'] = bkgd_re_hits

#         if weights is not None:       
#             #target_weighted_re_hits = [np.array([1 if (x is not None) else 0 for x in l]) * (1/weights) for l in target_regex_hits] #1*[np.sum(1*(tr.rw_beta[i,:] <= t) * (1/weights)) for t,i in zip(max_radi, range(0,tr.rw_beta.shape[0]))]
#             #centers_df['target_weighted_re_hits'] = target_weighted_re_hits
#             bkgd_weighted_re_hits = [np.sum(np.array([1 if (x is not None) else 0 for x in l]) * (1/weights)) for l in bkgd_regex_hits] 
#             centers_df['bkgd_weighted_re_hits'] = bkgd_weighted_re_hits 

#     n1 = tr.clone_df.shape[0]
#     # Total bakground is just the number of clones on background
#     n2 = tr_background.clone_df.shape[0]
#     centers_df['TR'] =  [compute_rate(pos=r['target_hits'], neg=n1-r['target_hits']) for i,r in centers_df.iterrows()]
#     centers_df['BR'] = [compute_rate(pos=r['bkgd_hits'], neg= n2-r['bkgd_hits']) for i,r in centers_df.iterrows()]
#     centers_df['RR'] = centers_df['TR']/centers_df['BR']
#     centers_df['OR'] = [compute_odds_ratio(pos=r['target_hits'], 
#                                          neg=n1-r['target_hits'], 
#                                          bpos=r['bkgd_hits'], 
#                                          bneg= n2-r['bkgd_hits'], ps = 1) for i,r in centers_df.iterrows()]
#     centers_df['ctrl'] = centers_df['bkgd_hits'] / n2

#     if weights is not None:
#         # Total bakground sum of the inverse of weightes
#         n1 = tr.clone_df.shape[0]
#         n2 = np.sum(1/weights)
#         centers_df['TR_weighted'] = [compute_rate(pos=r['target_weighted_hits'], neg=n1-r['target_weighted_hits']) for i,r in centers_df.iterrows()]
#         centers_df['BR_weighted'] = [compute_rate(pos=r['bkgd_weighted_hits'], neg=n2-r['bkgd_weighted_hits']) for i,r in centers_df.iterrows()]
#         centers_df['RR_weighted'] = centers_df['TR_weighted']/centers_df['BR_weighted']
#         centers_df['OR_weighted'] =[compute_odds_ratio(pos=r['target_hits'], 
#                                                neg=n1-r['target_hits'], 
#                                                bpos=r['bkgd_weighted_hits'], 
#                                                bneg= n2-r['bkgd_weighted_hits'], ps = 1) for i,r in centers_df.iterrows()]
#         centers_df['ctrl_weighted'] = centers_df['bkgd_weighted_hits'] / n2


#     return centers_df
