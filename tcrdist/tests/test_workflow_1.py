import pytest 



def test_workflow_1():
    """
    Load all the TCRs associated with a particular epitope in 
    the Adaptive Biotechnology COVID19 Data Release 2
    """
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    path = os.path.join('tcrdist', 'data', 'covid19')
    file = 'mira_epitope_9_2477_FLQSINFVR_FLQSINFVRI_FLYLYALVYF_GLEAPFLYLY_INFVRIIMR_LQSINFVRI_LQSINFVRII_QSINFVRII_SINFVRIIMR_VYFLQSINF_VYFLQSINFV_YFLQSINFVR_YLYALVYFL.tcrdist3.csv'
    filename = os.path.join(path,file)
    
    df = pd.read_csv(filename, sep = ",")
   
    df = df[['cell_type','subject','v_b_gene','j_b_gene','cdr3_b_aa',
            'epitope', 'age', 'sex', 'race','cohort',
            'hla-a', 'hla-a_1', 'hla-b', 'hla-b_1', 'hla-c','hla-c_1', 
            'dpa1_1', 'dpb1', 'dpb1_1', 
            'dqa1', 'dqa1_1', 'dqb1','dqb1_1', 'drb1', 'drb1_1', 'drb3']]
    
    df['count'] = 1
    
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['beta'])
