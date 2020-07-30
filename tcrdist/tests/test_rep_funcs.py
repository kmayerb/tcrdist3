import pandas as pd
import pwseqdist as pw
import numpy as np
from tcrdist.rep_funcs import _pws, _pw

def test_pw_rectangular():
    df = pd.read_csv("dash.csv")
    rmat = _pw(
        metric = pw.metrics.nb_vector_tcrdist,
        seqs1 = df.cdr3_b_aa[0:10],
        seqs2 = df.cdr3_b_aa, 
        ncpus=1, 
        uniqify= True, 
        use_numba = True)
    assert rmat.shape == (10,1924)


def test_pws_rectangular_computation():
    metrics = { "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist}

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

    df = pd.read_csv("dash2.csv")
    df = df.head(10).copy()
    df2 = pd.read_csv("dash2.csv")
    r = _pws(df = df, df2 = df2, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    assert r['tcrdist'].shape == (10,1924)

def test_dash_tcrdist_fixed_gappos_False():
    import pandas as pd
    import pwseqdist as pw
    from tcrdist.rep_funcs import _pws
    
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

    df = pd.read_csv("dash2.csv")
    r = _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    assert r['tcrdist'].shape[0] == 1924
    assert r['tcrdist'].shape[1] == 1924

def test_dash_tcrdist_fixed_gappos_True():
    import pandas as pd
    import pwseqdist as pw
    from tcrdist.rep_funcs import _pws
    
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

    kargs = {   "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':True},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':True},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}

    df = pd.read_csv("dash2.csv")
    r = _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    assert r['tcrdist'].shape[0] == 1924
    assert r['tcrdist'].shape[1] == 1924
    

def test_dash_nw_metric_fixed_gappos_False():
    import pandas as pd
    import pwseqdist as pw
    from tcrdist.rep_funcs import _pws
    
    # Define metrics for each region
    metrics = { "cdr3_a_aa" : pw.metrics.nw_metric,
                "pmhc_a_aa" : pw.metrics.nw_metric,
                "cdr2_a_aa" : pw.metrics.nw_metric,
                "cdr1_a_aa" : pw.metrics.nw_metric,
                "cdr3_b_aa" : pw.metrics.nw_metric,
                "pmhc_b_aa" : pw.metrics.nw_metric,
                "cdr2_b_aa" : pw.metrics.nw_metric,
                "cdr1_b_aa" : pw.metrics.nw_metric}

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

    kargs = {   "cdr3_a_aa" : {'use_numba': False},
                "pmhc_a_aa" : {'use_numba': False},
                "cdr2_a_aa" : {'use_numba': False},
                "cdr1_a_aa" : {'use_numba': False},
                "cdr3_b_aa" : {'use_numba': False},
                "pmhc_b_aa" : {'use_numba': False},
                "cdr2_b_aa" : {'use_numba': False},
                "cdr1_b_aa" : {'use_numba': False}}

    df = pd.read_csv("dash2.csv")
    import multiprocessing
    df = df.head(100)
    r = _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    assert r['tcrdist'].shape[0] == 100
    assert r['tcrdist'].shape[1] == 100



def test_pw():
    df = pd.read_csv("dash2.csv")
    seqs = df['cdr3_b_aa'].values
    r = _pw(metric = pw.metrics.nb_vector_tcrdist, 
            seqs1= seqs,
            seqs2 = None,
            ncpus=1, 
            uniqify=True, 
            use_numba=True, 
            distance_matrix=pw.matrices.tcr_nb_distance_matrix, 
            dist_weight=3, 
            gap_penalty=4, 
            ntrim=3, 
            ctrim=2, 
            fixed_gappos=True)
            
    assert isinstance(r, np.ndarray)

def test_indirect_nw_hamming_metric():
    seqs = ['CASSLDRGEVFF', # Seq1
            'CASSLDRGEVFF', # Seq2 = Seq1 i.e., D(s1,s2) = 0
            'CASSLDAGEVFF', # Seq3 = Seq1 (subsitition A for R) i.e., D(s1,s3) = 1
            'CASSLRGEVFF']  # Seq4 = Seq1 (Delete D for R) i.e., D(s1,s4) = 1 and D(s3,s4) = 2

    result = _pw(metric = pw.metrics.nw_hamming_metric, 
                  seqs1 = seqs, 
                  seqs2 = None, 
                  ncpus=1, 
                  uniqify= True)
    expectation = np.array([[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 2],[1, 1, 2, 0]])
    assert np.all(result == expectation )

def test_indirect_nw_metric():

    seqs = ['CASSLDRGEVFF', # Seq1
            'CASSLDRGEVFF', # Seq2 = Seq1 i.e., D(s1,s2) = 0
            'CASSLDAGEVFF', # Seq3 = Seq1 (subsitition A for R) i.e., D(s1,s3) = 1
            'CASSLRGEVFF' ] # Seq4 = Seq1 (Delete D for R) i.e., D(s1,s4) = 1 and D(s3,s4) = 2

    result = _pw(metric  = pw.metrics.nw_metric, 
                  seqs1   = seqs, 
                  seqs2   = None, 
                  ncpus   = 1, 
                  uniqify = True)

    assert isinstance(result, np.ndarray)

def test_direct_nw_hamming_metric():
    """
    This example uses nw_hamming_metric. It aligns two sequences and then calculates number of
    mismatched positions. 
    """
    
    seqs = ['CASSLDRGEVFF', # Seq1
            'CASSLDRGEVFF', # Seq2 = Seq1 i.e., D(s1,s2) = 0
            'CASSLDAGEVFF', # Seq3 = Seq1 (subsitition A for R) i.e., D(s1,s3) = 1
            'CASSLRGEVFF']  # Seq4 = Seq1 (Delete D for R) i.e., D(s1,s4) = 1 and D(s3,s4) = 2
    
    result      = pw.apply_pairwise_rect(metric = pw.metrics.nw_hamming_metric, seqs1 = seqs, uniqify= False, ncpus=1)
    expectation = np.array([[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 2],[1, 1, 2, 0]])
    assert np.all(result == expectation )


def test_indirect_nw_metric_uniquify_returns_full_dimension():

    seqs = ['CASSLDRGEVFF', # Seq1
            'CASSLDRGEVFF', # Seq2 = Seq1 i.e., D(s1,s2) = 0
            'CASSLDAGEVFF', # Seq3 = Seq1 (subsitition A for R) i.e., D(s1,s3) = 1
            'CASSLRGEVFF' ]  # Seq4 = Seq1 (Delete D for R) i.e., D(s1,s4) = 1 and D(s3,s4) = 2

    result = _pw(metric  = pw.metrics.nw_metric, 
                  seqs1   = seqs, 
                  seqs2   = None, 
                  ncpus   = 1, 
                  uniqify = True)

    assert result.shape[1] == 4
    assert isinstance(result, np.ndarray)
