import pytest 



def test_example_2():
    """  
    If you want 'tcrdistances' with changes over some parameters.
    For instance you want to change the gap penalty on CDR3s to 5. 
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

    tr.kargs_a['cdr3_a_aa']['gap_penalty'] = 5 
    tr.kargs_b['cdr3_b_aa']['gap_penalty'] = 5 
    
    tr.compute_distances()
    
    tr.pw_alpha
    tr.pw_beta