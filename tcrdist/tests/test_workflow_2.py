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
    
    tr.clone_df['covid'] = ['healthy' if x.find("Healthy") != -1 else "covid" for x in tr.clone_df.cohort]
    
    from tcrdist.rep_diff import neighborhood_diff, hcluster_diff, member_summ
    import hierdiff
    #tr.clone_df['covid'] = ['healthy' if x.find("Healthy") != -1 else "covid" for x in tr.clone_df.cohort]
    #nd = neighborhood_diff(tr.clone_df, tr.pw_beta, x_cols = ['covid'], count_col = 'count')

    tr.clone_df['covid'] = ['healthy' if x.find("Healthy") != -1 else "covid" for x in tr.clone_df.cohort]
    res, Z= hcluster_diff(tr.clone_df, tr.pw_beta, x_cols = ['covid'], count_col = 'count')
    
    res_summary = member_summ(res_df = res, clone_df = tr.clone_df, addl_cols=['cohort','subject'])
    
    res_detailed = pd.concat([res, res_summary], axis = 1)
    
    html = hierdiff.plot_hclust_props(Z,
                title='PA Epitope Example',
                res=res_detailed,
                tooltip_cols=['cdr3_b_aa','v_b_gene', 'j_b_gene','cohort','subject'],
                alpha=0.05, 
                alpha_col='pvalue')
    
    with open('hierdiff_example.html', 'w') as fh:
        fh.write(html)
