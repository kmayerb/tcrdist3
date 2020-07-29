import pytest



def test_example_9():
    """
    If you already have a clones file and want 
    to compute 'tcrdistances' on a DataFrame with 
    custom columns names.
    
    Set:
    1. Assign TCRrep.clone_df
    2. set infer_cdrs = False,
    3. compute_distances = False
    4. deduplicate = False
    5. customize the keys for metrics, weights, and kargs with the lambda
        customize = lambda d : {new_cols[k]:v for k,v in d.items()} 
    6. call .calculate_distances()
    """
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    
    new_cols = {'cdr3_a_aa':'c3a', 'pmhc_a_aa':'pa', 'cdr2_a_aa':'c2a','cdr1_a_aa':'c1a',
                'cdr3_b_aa':'c3b', 'pmhc_b_aa':'pb', 'cdr2_b_aa':'c2b','cdr1_b_aa':'c1b'}
    
    df = pd.read_csv("dash2.csv").rename(columns = new_cols) 
    
    tr = TCRrep(
            cell_df = df,
            clone_df = df,              #(1)
            organism = 'mouse', 
            chains = ['alpha','beta'],
            infer_all_genes = True, 
            infer_cdrs = False,         #(2)s
            compute_distances = False,  #(3)
            deduplicate=False,          #(4)
            db_file = 'alphabeta_gammadelta_db.tsv')

    customize = lambda d : {new_cols[k]:v for k,v in d.items()} #(5)
    tr.metrics_a = customize(tr.metrics_a)
    tr.metrics_b = customize(tr.metrics_b)
    tr.weights_a = customize(tr.weights_a)
    tr.weights_b = customize(tr.weights_b)
    tr.kargs_a = customize(tr.kargs_a)
    tr.kargs_b = customize(tr.kargs_b)

    tr.compute_distances() #(6)
    
    # Notice that pairwise results now have custom names 
    tr.pw_c3b
    tr.pw_c3a
    tr.pw_alpha
    tr.pw_beta
    