import pytest 



def test_workflow_2():
    """
    Load all the TCRs associated with a particular epitope in 
    the Adaptive Biotechnology COVID19 Data Release 2
    """
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.adpt_funcs import get_basic_centroids

    path = os.path.join('tcrdist', 'data', 'covid19')
    file = 'mira_epitope_16_1683_QYIKWPWYI_YEQYIKWPW_YEQYIKWPWY.tcrdist3.csv'
    filename = os.path.join(path,file)
    
    df = pd.read_csv(filename, sep = ",")
   
    df = df[['cell_type','subject','v_b_gene','j_b_gene','cdr3_b_aa',
            'epitope', 'age', 'sex','cohort']]
            
    df['count'] = 1
    
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['beta'])
 
    tr = get_basic_centroids(tr, max_dist = 200)
    
    tr.centroids_df
