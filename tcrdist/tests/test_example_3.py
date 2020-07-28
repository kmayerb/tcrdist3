import pytest



def test_example_3():
    """
    If want a 'tcrdistances' AND you want control over EVERY parameter.
    """
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    
    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                compute_distances = False,
                db_file = 'alphabeta_gammadelta_db.tsv')

    metrics_a = {
        "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
        "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist}

    metrics_b = {
        "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
        "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
        "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist }

    weights_a= { 
        "cdr3_a_aa" : 3,
        "pmhc_a_aa" : 1,
        "cdr2_a_aa" : 1,
        "cdr1_a_aa" : 1}

    weights_b = { 
        "cdr3_b_aa" : 3,
        "pmhc_b_aa" : 1,
        "cdr2_b_aa" : 1,
        "cdr1_b_aa" : 1}

    kargs_a = {  
        'cdr3_a_aa' : 
            {'use_numba': True, 
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 
            'dist_weight': 1, 
            'gap_penalty':4, 
            'ntrim':3, 
            'ctrim':2, 
            'fixed_gappos': False},
        'pmhc_a_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True},
        'cdr2_a_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight': 1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True},
        'cdr1_a_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True}
        }
    kargs_b= {  
        'cdr3_b_aa' : 
            {'use_numba': True, 
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 
            'dist_weight': 1, 
            'gap_penalty':4, 
            'ntrim':3, 
            'ctrim':2, 
            'fixed_gappos': False},
        'pmhc_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight': 1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True},
        'cdr2_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True},
        'cdr1_b_aa' : {
            'use_numba': True,
            'distance_matrix': pw.matrices.tcr_nb_distance_matrix,
            'dist_weight':1,
            'gap_penalty':4,
            'ntrim':0,
            'ctrim':0,
            'fixed_gappos':True}
        }   

    tr.metrics_a = metrics_a
    tr.metrics_b = metrics_b

    tr.weights_a = weights_a
    tr.weights_b = weights_b

    tr.kargs_a = kargs_a 
    tr.kargs_b = kargs_b

    tr.compute_distances()
    tr.pw_alpha
    tr.pw_beta