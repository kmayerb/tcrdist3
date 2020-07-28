import pytest



def test_example_4():
    """
    If you want "tcrdistances" using a different metric.

    Here we illustrate the use a metric that uses the 
    Needleman-Wunsch algorithm to align sequences and then 
    calculate the number of mismatching positions (pw.metrics.nw_hamming_metric)

    This method doesn't rely on Numba so it can run faster using multiple cpus.
    """
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    import multiprocessing

    df = pd.read_csv("dash.csv")
    df = df.head(100) # for faster testing
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                use_defaults=False,
                compute_distances = False,
                cpus = 1,
                db_file = 'alphabeta_gammadelta_db.tsv')

    metrics_a = {
        "cdr3_a_aa" : pw.metrics.nw_hamming_metric ,
        "pmhc_a_aa" : pw.metrics.nw_hamming_metric ,
        "cdr2_a_aa" : pw.metrics.nw_hamming_metric ,
        "cdr1_a_aa" : pw.metrics.nw_hamming_metric }

    metrics_b = {
        "cdr3_b_aa" : pw.metrics.nw_hamming_metric ,
        "pmhc_b_aa" : pw.metrics.nw_hamming_metric ,
        "cdr2_b_aa" : pw.metrics.nw_hamming_metric ,
        "cdr1_b_aa" : pw.metrics.nw_hamming_metric  }

    weights_a = { 
        "cdr3_a_aa" : 1,
        "pmhc_a_aa" : 1,
        "cdr2_a_aa" : 1,
        "cdr1_a_aa" : 1}

    weights_b = { 
        "cdr3_b_aa" : 1,
        "pmhc_b_aa" : 1,
        "cdr2_b_aa" : 1,
        "cdr1_b_aa" : 1}

    kargs_a = {  
        'cdr3_a_aa' : 
            {'use_numba': False},
        'pmhc_a_aa' : {
            'use_numba': False},
        'cdr2_a_aa' : {
            'use_numba': False},
        'cdr1_a_aa' : {
            'use_numba': False}
        }
    kargs_b = {  
        'cdr3_b_aa' : 
            {'use_numba': False},
        'pmhc_b_aa' : {
            'use_numba': False},
        'cdr2_b_aa' : {
            'use_numba': False},
        'cdr1_b_aa' : {
            'use_numba': False}
        }

    tr.metrics_a = metrics_a
    tr.metrics_b = metrics_b

    tr.weights_a = weights_a
    tr.weights_b = weights_b

    tr.kargs_a = kargs_a 
    tr.kargs_b = kargs_b

    tr.compute_distances()

    tr.pw_cdr3_b_aa
    tr.pw_beta
