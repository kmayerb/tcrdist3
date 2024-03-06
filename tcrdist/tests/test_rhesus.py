import pytest 



def test_example_rhesus():
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

    df = pd.read_csv("rhesus.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'rhesus', 
                chains = ['alpha','beta'], 
                db_file = 'combo_xcr_2023-12-30.tsv')
    
    tr.pw_alpha
    tr.pw_beta
    tr.pw_cdr3_a_aa
    tr.pw_cdr3_b_aa


