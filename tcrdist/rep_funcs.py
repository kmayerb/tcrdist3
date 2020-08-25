"""
tcrdist3: Functional Programming Approach For Higher Memory Demand Use Cases
"""
import pwseqdist as pw
import pandas as pd
import numpy as np

__all__ = ['_pw', '_pws']


def _pws(df, metrics, weights, kargs, df2 = None, cpu = 1, uniquify = True, store = False):
    """    
    _pws performs pairwise distance calculation across a multiple 
    columns of a Pandas DataFrame. This naturally permits calculation of
    a CDR-weighted tcrdistance that incorporates dissimilarity across 
    multiple complimentarity determining regions (see example below):

    Parameters 
    ----------
    df : pd.DataFrame
        Clones DataFrame containing, at a minimum, columns with CDR sequences
    df2 : pd.DataFrame or None
        Second clones DataFrame containing, at a minimum, columns with CDR sequences
    metrics : dict
        Dictionary of functions, specifying the distance metrics to apply to each CDR
    weights : dict
        Weights determining the contributions of each CDR distance to the aggregate distance
    kargs : dict
        Dictionary of Dictionaries
    cpu : int 
        Number of available cpus
    use_numba : bool
        If True, use must use a numba compatible metric 
    store : bool
        If False, only full tcrdist is returned. If True, 
        all component distance matrices are returned.
        
    Returns
    -------
    s : dictionary with tcr_distance.

    Example
    -------
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.rep_funcs import _pw, _pw2
    
    # Define metrics for each region
    metrics = { "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist}

    # Define weights
    weights = { 
                "cdr3_a_aa" : 3,
                "pmhc_a_aa" : 1,
                "cdr2_a_aa" : 1,
                "cdr1_a_aa" : 1,
                "cdr3_b_aa" : 3,
                "pmhc_b_aa" : 1,
                "cdr2_b_aa" : 1,
                "cdr1_b_aa" : 1}

    kargs = {   "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}

    df = pd.DataFrame("dash2.csv")
    _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    """
    metric_keys = list(metrics.keys())
    weight_keys = list(weights.keys())
    assert metric_keys == weight_keys, "metrics and weights keys must be identical"
    
    if kargs is not None:
        kargs_keys  = list(kargs.keys())
        assert metric_keys == kargs_keys,  "metrics and kargs keys must be identical"
    
    tcrdist = None
    s = dict()
    for k in metric_keys:
        if df2 is None:
            pw_mat = _pw(seqs1 = df[k].values, metric = metrics[k], ncpus = cpu, uniqify= uniquify, **kargs[k])
        else:
            pw_mat = _pw(seqs1 = df[k].values, seqs2 = df2[k].values, metric = metrics[k], ncpus = cpu, uniqify= uniquify, **kargs[k])
            
        if store:
           s[k] = pw_mat 
        if tcrdist is None:
            tcrdist = np.zeros(pw_mat.shape, dtype=np.int16)
        tcrdist = tcrdist + (weights[k] * pw_mat)
    
    s['tcrdist'] = tcrdist
    return s


def _pw(metric, seqs1, seqs2=None, ncpus=1, uniqify= True, use_numba = False, **kwargs):
    """
    This is a wrapper for accessing pwseqdist version > 0.2.
    No matter what, it returns squareform results
    """
    pw_mat = pw.apply_pairwise_rect(metric = metric, 
                                    seqs1  = seqs1, 
                                    seqs2  = seqs2, 
                                    ncpus  = ncpus, 
                                    uniqify = uniqify, 
                                    use_numba = use_numba,
                                    **kwargs)
    # if len(pw_mat.shape) == 1:
    #     from scipy.spatial.distance import squareform
    #     pw_mat = squareform(pw_mat)

    return pw_mat

