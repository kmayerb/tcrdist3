import pytest 



def test_example_1():
    """
    If you just want a 'tcrdistances' using pre-set default setting.
    
        You can access distance matrices:
            tr.pw_alpha     - alpha chain pairwise distance matrix
            tr.pw_beta      - alpha chain pairwise distance matrix
            tr.pw_cdr3_a_aa - cdr3 alpha chain distance matrix
            tr.pw_cdr3_b_aa - cdr3 beta chain distance matrix
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    
    tr.pw_alpha
    tr.pw_beta
    tr.pw_cdr3_a_aa
    tr.pw_cdr3_b_aa


