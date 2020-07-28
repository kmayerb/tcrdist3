import pytest



def test_example_8():
    """
    You want a 'tcrdistance' but you don't want to bother with the tcrdist3 framework. 
    
    Note that the columns names are completely arbitrary under this 
    framework, so one can directly compute a tcrdist on a 
    AIRR, MIXCR, VDJTools, or other formated file without any
    reformatting.
    """ 
    import multiprocessing
    import pandas as pd
    import pwseqdist as pw
    from tcrdist.rep_funcs import _pws, _pw  
    
    df_airr = pd.read_csv("dash_beta_airr.csv")

    # Choose the metrics you want to apply to each CDR
    metrics = { 'cdr3_aa' : pw.metrics.nb_vector_tcrdist,
                'cdr2_aa' : pw.metrics.nb_vector_tcrdist,
                'cdr1_aa' : pw.metrics.nb_vector_tcrdist}
    
    # Choose the weights that are right for you.
    weights = { 'cdr3_aa' : 3,
                'cdr2_aa' : 1,
                'cdr1_aa' : 1}

    # Provide arguments for the distance metrics 
    kargs = {   'cdr3_aa' : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                'cdr2_aa' : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                'cdr1_aa' : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
                
    # Here are your distance matrices
    from tcrdist.rep_funcs import _pws
    
    dmats = _pws(df = df_airr,
             metrics = metrics, 
             weights= weights, 
             kargs=kargs, 
             cpu = 1, 
             store = True)

    dmats['tcrdist']
