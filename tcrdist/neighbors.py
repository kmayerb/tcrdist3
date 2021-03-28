"""
tcrdist/neigbhors.py

This module contas 
"""

import numpy as np 
import pandas as pd
import parmap
import scipy.sparse
import warnings
#from tcrdist.repertoire import TCRrep
from tcrdist.regex import _index_to_regex_str, _index_to_seqs, _multi_regex 
from scipy.stats import chi2_contingency



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

def _todense_row(pw, maxd):
    """Make a sparse distance matrix dense, dealing with zeros as > maxd and -1 as true zero"""
    pw = np.asarray(pw.todense())
    pw[pw == 0] = maxd + 1
    pw[pw == -1] = 0
    return pw

#      888                                          888                                           888 d8b          888                                                                    888  .d888 
#      888                                          888                                           888 Y8P          888                                                                    888 d88P"  
#      888                                          888                                           888              888                                                                    888 888    
#  .d88888  .d88b.  88888b.   .d88b.  88888b.   .d88888 .d8888b        .d88b.  88888b.        .d88888 888 .d8888b  888888  8888b.  88888b.   .d8888b .d88b.          .d88b.   .d8888b .d88888 888888 
# d88" 888 d8P  Y8b 888 "88b d8P  Y8b 888 "88b d88" 888 88K           d88""88b 888 "88b      d88" 888 888 88K      888        "88b 888 "88b d88P"   d8P  Y8b        d8P  Y8b d88P"   d88" 888 888    
# 888  888 88888888 888  888 88888888 888  888 888  888 "Y8888b.      888  888 888  888      888  888 888 "Y8888b. 888    .d888888 888  888 888     88888888        88888888 888     888  888 888    
# Y88b 888 Y8b.     888 d88P Y8b.     888  888 Y88b 888      X88      Y88..88P 888  888      Y88b 888 888      X88 Y88b.  888  888 888  888 Y88b.   Y8b.            Y8b.     Y88b.   Y88b 888 888    
#  "Y88888  "Y8888  88888P"   "Y8888  888  888  "Y88888  88888P'       "Y88P"  888  888       "Y88888 888  88888P'  "Y888 "Y888888 888  888  "Y8888P "Y8888 88888888 "Y8888   "Y8888P "Y88888 888    
#                   888                                                                                                                                                                              
#                   888                                                                                                                                                                              
#                   888   

def bkgd_cntl_nn2( tr, 
                   tr_background,
                   weights = None,
                   ctrl_bkgd = 10**-5, 
                   col = 'cdr3_b_aa',
                   add_cols = ['v_b_gene', 'j_b_gene'],
                   pw_mat_str = 'pw_beta',
                   rw_mat_str = 'rw_beta',
                   ncpus = 4,
                   include_seq_info= True,
                   thresholds = [x for x in range(0,50,2)],
                   generate_regex = True,
                   test_regex = True,
                   beta_re = 1,
                   beta_dist = 1,
                   forced_max_radius = None):   

    """
    bkgd_cntl_nn2 stands for background controlled nearest neighbors. 

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
    ncpus : int
        e.g., 4, passed to pm_processes when using parmap   
    include_seq_info : bool 
        If True, returned DataFrame <centers_df> will include ['target_neighbors', 'target_seqs',
        'background_neighbors','background_seqs','background_v', 'background_j']
        as columns in centers_df DataFrame returned by this function. 
        This allows for inspection of sequences found in both the antigen enriched repertoire and supplied
        background.
    thresholds : list
        Default is [x for x in range(0,50,2)] indicating tcrdist thresholds to compute 
        a ecdf over.
    generate_regex : bool
        if True, generate a regex pattern capturing most if not all of 
        sequences in the neigbhorhood
    test_regex = True
        if True, test regex against both CDR3 in target and background
    beta_re : int or float
        for joint chi2, weight for regex based portion
    beta_dist : int of float
        for joint chi2, weight for distance based portion
    forced_max_radius : int or None
        if not None, radius cannot exceed this amount
    
    Returns 
    -------
    centers_df: DataFrame
    
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

        # <cdr_col> usually 'cdr3_b_aa', we save this for later 
    cdr3_col = str(col)
        # <pw_mat> pairwise matrix, typically the network between antigen enriched clones
    pw_mat = getattr(tr, pw_mat_str)# defaults to 'pw_beta',
        # <rw_mat> rectangular matrix, typically the network betwen antigen enriched and bulk clones
    rw_mat = getattr(tr, rw_mat_str)# defaults to 'rw_beta',
        # remind the user if they are using sparse or dense matrices as input
    print(f"PW MATRIX TYPE = {type(pw_mat)}")
    print(f"RW MATRIX TYPE = {type(rw_mat)}")
        # <max_radius> int describing the highest distances, only will be needed for sparse impementation
        # <thresholds> array of thresholds radii to consider 
    thresholds = np.unique(thresholds)
    max_radius = thresholds.max() + 1
        # <weights> inverse probability weights, if none provided assume all are one
    if weights is None:
        warnings.warn("No weights provided trying to use tr_background.clone_df.weights")
        if 'weights' in tr_background.clone_df.columns:
            warnings.warn("USING tr_background.clone_df.weights")
            weights = tr_background.clone_df['weights']
        else:
            warnings.warn("SOMETHING IS PROBABLY WRONG!!!!. You don't have a 'weights' column in your tr_background clone_df, setting to 1, results could be biased if this data is V-J matched")
            weights = np.ones(rw_mat.shape[1])
    thresholds = np.unique(thresholds)
        # for each TCR, we calculate a empirical cummulative 
        # density function along a range of threshold radii



    # 88888b.   .d88b.  888  888  888 
    # 888 "88b d8P  Y8b 888  888  888 
    # 888  888 88888888 888  888  888 
    # 888  888 Y8b.     Y88b 888 d88P 
    # 888  888  "Y8888   "Y8888888P" 

    # KMB: I’ve confirmed that the distance_ecdf(absolute_weight = True) and my function provide same ecdfs. Therefore distance_ecdf can be our ONLY ecdf function (which is nice because it will be used for both plots and results).  My function now relies on it. I have commented out my code but will remove completely in a later version.
    
    from tcrdist.ecdf import distance_ecdf
    thresholds, ecdfs2 = distance_ecdf(pwrect =rw_mat, 
      thresholds = thresholds, 
      weights= weights, 
      pseudo_count=0, 
      skip_diag = False, 
      absolute_weight = True)
    
    ecdfs2 = [pd.Series(x, index = thresholds) for x in ecdfs2]
    max_radi2 = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs2]
    # TO SOLVE NOTABLE BUG IN THE ABOVE LINE!! If a radius is None (the next line will fail, thus set Nones to 0.
    max_radi2 = [x if (x is not None) else 0 for x in max_radi2]

    if forced_max_radius is not None:
        max_radi2 = [min(x,forced_max_radius) for x in max_radi2]

    #          888      888 
    #          888      888 
    #          888      888 
    #  .d88b.  888  .d88888 
    # d88""88b 888 d88" 888 
    # 888  888 888 888  888 
    # Y88..88P 888 Y88b 888 
    #  "Y88P"  888  "Y88888 

    # print("USING : _compute_weighted_pop_estimate_ecdf_rowwise")
    # ecdfs = parmap.map(_compute_weighted_pop_estimate_ecdf_rowwise, range(0,rw_mat.shape[0]), 
    #    data = rw_mat, 
    #    weights = weights,
    #    thresholds = thresholds, 
    #    pm_pbar = True, 
    #    pm_processes = ncpus)
    # #<max_radi> Based on acceptable ctrl_bkgd, we find max acceptable radius from each TCR
    # max_radi = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs]
    # # TO SOLVE NOTABLE BUG IN THE ABOVE LINE!! If a radius is None (the next line will fail, thus set Nones to 0.
    # max_radi = [x if (x is not None) else 0 for x in max_radi]  
    
    # if forced_max_radius is not None:
    #     max_radi = [min(x,forced_max_radius) for x in max_radi]
    # old vs. new methods must yield same results, 
    # TODO: remove old in next version, so there is a commit record
    #import pdb; pdb.set_trace()
    #assert np.all(max_radi2 == max_radi)
    #print("SUCCESS COMPARING METHODS")
    max_radi = max_radi2
        # <target_hits> number of hits within the antigen enriched repertoire at the control radius 
    if scipy.sparse.issparse(pw_mat):
        # sparse implementation, relies on _todense_row (see above in this module)
        target_hits    = [np.sum(_todense_row(pw_mat[i,:], max_radius )[0] <= t) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
    else: 
        target_hits    = [np.sum(pw_mat[i,:] <= t) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
        # <bkgd_hits>  number of hits within the bulk unenriched repertoire at the control radius 
        # <bkgd_hits_weighted> number of hits when inverse probability weights are applied
    if scipy.sparse.issparse(rw_mat):
        bkgd_hits          = [np.sum(_todense_row(rw_mat[i,:], max_radius)[0] <= t) for t,i in zip(max_radi, range(0,rw_mat.shape[0]))]
        bkgd_hits_weighted = [np.sum(1*(_todense_row(rw_mat[i,:], max_radius)[0] <= t) * (weights)) for t,i in zip(max_radi, range(0,rw_mat.shape[0]))]
    else:
        bkgd_hits      = [np.sum(rw_mat[i,:] <= t) for t,i in zip(max_radi, range(0,rw_mat.shape[0]))]
        bkgd_hits_weighted = [np.sum(1*(rw_mat[i,:] <= t) * (weights)) for t,i in zip(max_radi, range(0,rw_mat.shape[0]))]

        # <bkdg_total> number of rows in background
    bkgd_total = rw_mat.shape[1]
        # <bkgd_weighted_total > 
    bkgd_weighted_total = rw_mat.shape[1] # In modern, implemenation demoninator is simply length not np.sum(weights)
        # <ctrl> divide total hits by total posible hits
    ctrl = np.array(bkgd_hits) / bkgd_total 
        # <ctrl_weighted> divide weighted hits byt total possible hits to approximate expected frequency in a real repertoire
    ctrl_weighted= np.array(bkgd_hits_weighted)/bkgd_weighted_total 
        # <center_df> this is the kernal of dataframe
    
        # we typically want pgen in our analysis so
    if f'pgen_{col}' in tr.clone_df.columns:
        centers_df = pd.DataFrame({col: tr.clone_df[col].copy(),
                  add_cols[0] : tr.clone_df[add_cols[0]].copy(), 
                  add_cols[1] : tr.clone_df[add_cols[1]].copy(), 
                  'pgen' : tr.clone_df[f'pgen_{col}'],
                  'radius': max_radi,
                  'target_hits' : target_hits, 
                  'bkgd_hits': bkgd_hits,
                  'bkgd_hits_weighted': bkgd_hits_weighted,
                  'bkgd_total' :  bkgd_total,
                  'ctrl' : ctrl,
                  'ctrl_weighted' : ctrl_weighted})
    else: 
        print("YOU DID NOT PRECALCULATE OPTIONAL PGENS SEE auto_pgen(TCRrep) NEXT TIME")
        centers_df = pd.DataFrame({col: tr.clone_df[col].copy(),
                  add_cols[0] : tr.clone_df[add_cols[0]].copy(), 
                  add_cols[1] : tr.clone_df[add_cols[1]].copy(), 
                  'radius': max_radi,
                  'target_hits' : target_hits, 
                  'bkgd_hits': bkgd_hits,
                  'bkgd_hits_weighted': bkgd_hits_weighted,
                  'bkgd_total' :  bkgd_total,
                  'ctrl' : ctrl,
                  'ctrl_weighted' : ctrl_weighted})

    # Total bakground counts
    n1 = tr.clone_df.shape[0]   
    n2 = bkgd_weighted_total 
    print(f"ANTIGEN-ENRICHED CLONES : {n1}")
    print(f"BULK  BACKGROUND CLONES : {n2}")  
    print("COMPUTING WEIGHTED ODDS RATIO, RELATIVE RATE")
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
    print("COMPUTING CHI-SQUARE STATISTIC BASED ON WEIGHTED DISTANCE HITS <chi2dist>")
    centers_df['chi2dist'] = [chi2_contingency(np.array(
            [[1+r['target_hits'], 1+n1-r['target_hits']],[1+r['bkgd_hits_weighted'], 1+n2-r['bkgd_hits_weighted']]]))[0] for _,r in centers_df.iterrows() ]

    if include_seq_info:
        print("INCLUDING SEQ INFO")
        if scipy.sparse.issparse(pw_mat):
            target_neighbors = [list((_todense_row(pw_mat[i,:], max_radius)[0] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
        else:
            target_neighbors = [list((pw_mat[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
        
        centers_df['target_neighbors'] = target_neighbors 
        centers_df['target_seqs'] = [tr.clone_df.loc[ar,col].to_list() for ar in target_neighbors]
        if scipy.sparse.issparse(rw_mat):
            background_neighbors = [list((_todense_row(rw_mat[i,:], max_radius)[0] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
        else:
            background_neighbors = [list((rw_mat[i,:] <= t).nonzero()[0]) for t,i in zip(max_radi, range(0,pw_mat.shape[0]))]
        

        centers_df['background_neighbors'] = background_neighbors
        #col = "cdr3_b_aa"
        centers_df['background_seqs'] = [tr_background.clone_df.loc[ar,col].to_list() for ar in background_neighbors]
        vcol = add_cols[0] #v_b_gene 
        centers_df['background_v'] = [tr_background.clone_df.loc[ar,vcol].to_list() for ar in background_neighbors]
        jcol = add_cols[1] #j_b_gene
        centers_df['background_j'] = [tr_background.clone_df.loc[ar,jcol].to_list() for ar in background_neighbors]
        adjcol  = 'adj_freq_pVJ'
        if adjcol in tr_background.clone_df.columns:
            centers_df['adj_freq'] = [tr_background.clone_df.loc[ar,adjcol].to_list() for ar in background_neighbors]
    
    if generate_regex:
        print("GENERATING REGEX FOR EACH META-CLONOTYPE")
        #cdr3_col = "cdr3_b_aa"
        centers_df['regex'] = [_index_to_regex_str(ind = r['target_neighbors'], 
                clone_df = tr.clone_df, 
                pwmat = None, 
                col = cdr3_col, 
                centroid = tr.clone_df[cdr3_col][i],
                ntrim = 3,
                ctrim = 2,
                max_ambiguity = 5) for i,r in centers_df.iterrows()]

    if test_regex:
        print("TESTING REGEX FOR EACH META-CLONOTYPE")
        # THIS INVOLVES TESTING REGEX AGAINST TARGET AND BACKGROUND
        # Test regex against backgound
        rs = centers_df['regex'].to_list()
        bkgd_cdr3 = tr_background.clone_df[cdr3_col].to_list()
        target_cdr3 = tr.clone_df[cdr3_col].to_list()
        target_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = target_cdr3 , pm_pbar = True, pm_processes = ncpus)
        target_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in target_regex_hits ]
        centers_df['target_re_hits'] = target_re_hits
        bkgd_regex_hits = parmap.map(_multi_regex, rs, bkgd_cdr3 = bkgd_cdr3, pm_pbar = True, pm_processes = ncpus)
        bkgd_re_hits = [np.sum([1 if (x is not None) else 0 for x in l]) for l in bkgd_regex_hits]
        centers_df['bkgd_re_hits'] = bkgd_re_hits
        bkgd_weighted_re_hits = [np.sum(np.array([1 if (x is not None) else 0 for x in l]) * (weights)) for l in bkgd_regex_hits] 
        # Compute Relative Rates
        print("COMPUTING WEIGHTED ODDS RATIO, RELATIVE RATES FOR EACH REGEX")
        centers_df['bkgd_re_weighted_hits'] = bkgd_weighted_re_hits 
        centers_df['TR_re'] = [compute_rate(pos=r['target_re_hits'], neg=n1-r['target_re_hits']) for i,r in centers_df.iterrows()]
        centers_df['BR_re_weighted'] = [compute_rate(pos=r['bkgd_re_weighted_hits'], 
                                     neg=n2-r['bkgd_re_weighted_hits']) for i,r in centers_df.iterrows()]
        centers_df['RR_re_weighted'] = centers_df['TR']/centers_df['BR_re_weighted']
        centers_df['OR_re_weighted'] =[compute_odds_ratio(pos=r['target_re_hits'], 
                                       neg=n1-r['target_re_hits'], 
                                       bpos=r['bkgd_re_weighted_hits'], 
                                       bneg= n2-r['bkgd_re_weighted_hits'], ps = 1) for i,r in centers_df.iterrows()]
        print("COMPUTING CHI-SQUARE STATISTIC BASED ON WEIGHTED REGEX HITS <chi2re>")
        centers_df['chi2re'] = [chi2_contingency(np.array(
            [[1+r['target_re_hits'], 1+n1-r['target_re_hits']],[1+r['bkgd_re_weighted_hits'], 1+n2-r['bkgd_re_weighted_hits']]]))[0] for _,r in centers_df.iterrows() ]
        print(f"COMPUTING JOINT CHI-SQUARE STATISTIC <chi2joint>BASED ON WEIGHTED REGEX HITS {beta_re} * <chi2re> + {beta_dist} * <chi2dist")
        centers_df['chi2joint'] = [beta_re  * r['chi2re'] + beta_dist* r['chi2dist'] for _,r in centers_df.iterrows() ]



    return centers_df #, thresholds, ecdfs











#          888      888                       888                                            d8b          888                 888 
#          888      888                       888                                            Y8P          888                 888 
#          888      888                       888                                                         888                 888 
#  .d88b.  888  .d88888                   .d88888  .d88b.  88888b.  888d888 .d88b.   .d8888b 888  8888b.  888888 .d88b.   .d88888 
# d88""88b 888 d88" 888                  d88" 888 d8P  Y8b 888 "88b 888P"  d8P  Y8b d88P"    888     "88b 888   d8P  Y8b d88" 888 
# 888  888 888 888  888      888888      888  888 88888888 888  888 888    88888888 888      888 .d888888 888   88888888 888  888 
# Y88..88P 888 Y88b 888                  Y88b 888 Y8b.     888 d88P 888    Y8b.     Y88b.    888 888  888 Y88b. Y8b.     Y88b 888 
#  "Y88P"  888  "Y88888                   "Y88888  "Y8888  88888P"  888     "Y8888   "Y8888P 888 "Y888888  "Y888 "Y8888   "Y88888 
#                                                          888                                                                    
#                                                          888                                                                    
#                                                          888                                                        
# 
# THIS WAS ONCE USED FOR COMPUTING ECDFS BUT WE CAN NOW RELY ONLY ON distance_ecdf

# def compute_population_estimate_ecdf(data, weights = None, thresholds=None):
#     """
#     """
#     if thresholds is None:
#         thresholds = np.unique(data[:])
  
#     if weights is None:
#         weights = np.ones(data.shape[0])
#     else:
#         weights = np.array(weights)
    
#     """Vectorized and faster, using broadcasting for the <= expression"""
#     # Step 1 : compare a set of n distances to a set of say 25 threshold 
#         # (n dists x 1) broadcast against (1 x 25 thresholds) ==> (n, 25), where each column representd number of hits below the threshold 
#     M= 1*(data[:, None] <= thresholds[None, :])
#     # Step 2 : Adjust hit be relative weight (i.e. actual_freq_vj / sampled_freq_vj)
#         # Broadcast (n, 25) on weights (n, 1) ==> (n, 25)
#     M = M * weights[:,None]
#     # Step 3: Take row sums, producing the weighted sum at each threshold ==> (25,)
#     M =  np.sum(M, axis = 0)
#     # Step 4: Divide by the weighted total (i.e., the amount if all seqs were present, which because of enrichment is far greater than 1)
#     # print(M)   
#     ecdf = M / len(weights) 
#     return ecdf


# def _compute_weighted_pop_estimate_ecdf_rowwise(i, thresholds, weights, data, max_radius = 50):
#     """
#     This is a wrapper for compute_ecdf using parmap
#     """
#     if scipy.sparse.issparse(data):
#         row = data[i,:].toarray()[0]
#         row[row==0] = (max_radius + 1)
#         row[row==-1] = 0
#         ec = compute_population_estimate_ecdf(data = row,  weights = weights, thresholds = thresholds)
#     else:
#         ec = compute_population_estimate_ecdf(data = data[i,:],  weights = weights, thresholds = thresholds)
#     return pd.Series(ec, index = thresholds)

# def _compute_pop_estimate_ecdf_rowwise(i, thresholds, data, max_radius = 50):
#     """
#     This is a warpper for compute_ecdf using parmap
#     """
#     if scipy.sparse.issparse(data):
#         row = data[i,:].toarray()[0]
#         row[row==0] = (max_radius + 1)
#         row[row==-1] = 0
#         ec = compute_population_estimate_ecdf(data = row, thresholds = thresholds)
#     else:
#         ec = compute_population_estimate_ecdf(data = data[i,:], thresholds = thresholds)
#     return pd.Series(ec, index = thresholds)


              
#                   888                                           888                          888      888 
#                   888                                           888                          888      888 
#                   888                                           888                          888      888 
#  .d88b.  888  888 888888 888d888 .d88b.  88888b.d88b.   .d88b.  888 888  888         .d88b.  888  .d88888 
# d8P  Y8b `Y8bd8P' 888    888P"  d8P  Y8b 888 "888 "88b d8P  Y8b 888 888  888        d88""88b 888 d88" 888 
# 88888888   X88K   888    888    88888888 888  888  888 88888888 888 888  888 888888 888  888 888 888  888 
# Y8b.     .d8""8b. Y88b.  888    Y8b.     888  888  888 Y8b.     888 Y88b 888        Y88..88P 888 Y88b 888 
#  "Y8888  888  888  "Y888 888     "Y8888  888  888  888  "Y8888  888  "Y88888         "Y88P"  888  "Y88888 
#                                                                          888                              
#                                                                     Y8b d88P                              
#                                                                      "Y88P"                               

# ## SOON TO BE DEPRECIATED BELOW THIS POINT: DO NOT USE

# def compute_ecdf(data, weights = None, counts=None, thresholds=None):
#     """
#     Computes the empirical cumulative distribution function at pre-specified
#     thresholds. Assumes thresholds is sorted and should be unique.

#     data : np.array

#     threshold : np.array or None
    
#     counts : np.array or None

#     weights :  np.array or None
#          vector of W_vj  (i.e  actual_freq_vj / sampled_freq_vj)
    
#     Example
#     -------
#     >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
#     array([0.5])
#     >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), weights =np.array([1,1,1,1,1,1]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
#     array([0.5])
#     >>> compute_ecdf(data = np.array([1,2,3,4,5,6]), weights =np.array([10,1,1,1,1,1]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
#     array([0.8])

#     Notes
#     -----
    
#     Weights are used when we enriche for V,J combinations that are rare. 
#     For instance if some V,J pairing 
#     occurs at an actual frequency of 1/100 but we provide it at 1/10, 
#     we provide an enrichment weight 
#     (actual_freq_vj / sampled_freq_v)


#     -----
#     """

#     if thresholds is None:
#         thresholds = np.unique(data[:])
#     if counts is None:
#         counts = np.ones(data.shape)
    
#     if weights is None:
#         weights = np.ones(data.shape)
#     else:
#         weights = np.array(weights)
    
#     tot = np.sum(counts * weights) 
#     """Vectorized and faster, using broadcasting for the <= expression"""
#     # Step 1 : compare a set of n distances to a set of say 25 threshold 
#         # (n dists x 1) broadcast against (1 x 25 thresholds) ==> (n, 25), where each column representd number of hits below the threshold 
#     M= 1*(data[:, None] <= thresholds[None, :])
#     # Step 2 : Adjust hit be relative weight (i.e. actual_freq_vj / sampled_freq_vj)
#         # Broadcast (n, 25) on weights (n, 1) ==> (n, 25)
#     M = M * weights[:,None] * counts[:,None]
#     # Step 3: Take row sums, producing the weighted sum at each threshold ==> (25,)
#     M =  np.sum(M, axis = 0)
#     # Step 4: Divide by the weighted total (i.e., the amount if all seqs were present, which because of enrichment is far greater than 1)
#     ecdf = M / tot

#     return ecdf



# def _compute_ecdf_rowwise(i, thresholds, data):
#     """
#     This is a warpper for compute_ecdf using parmap
#     """
#     ec = compute_ecdf(data = data[i,:], thresholds = thresholds)
#     return pd.Series(ec, index = thresholds)

# def _compute_weighted_ecdf_rowwise(i, thresholds, weights, data):
#     """
#     This is a warpper for compute_ecdf using parmap
#     """
#     ec = compute_ecdf(data = data[i,:], weights = weights, thresholds = thresholds)
#     return pd.Series(ec, index = thresholds)







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
