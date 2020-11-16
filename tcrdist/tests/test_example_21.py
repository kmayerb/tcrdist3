



def test_TCRpublic_with_neighborhood_dif():
    """
    Use values from neighborhood_diff
    """
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic  
    fn = os.path.join('tcrdist','data','covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.radius.csv')
    df = pd.read_csv(fn)
    tr = TCRrep(cell_df = df[['cohort','subject','v_b_gene', 'j_b_gene','cdr3_b_aa', 'radius']], 
                organism = "human", 
                chains = ["beta"])
    
    from tcrdist.rep_diff import neighborhood_diff
    ndif = neighborhood_diff(   clone_df= tr.clone_df, 
                                    pwmat = tr.pw_beta, 
                                    count_col = 'count', 
                                    x_cols = ['cohort'], 
                                    knn_radius = 25, 
                                    test_method = "chi2")
    # Add neighbors and other columns of interest 
    # from neighbor_diff result to the clone_df
    tr.clone_df = pd.concat([tr.clone_df, ndif[['neighbors', 'K_neighbors','val_0','ct_0','pvalue']] ], axis = 1)
    # Because neighors and K_neighbor are already added to the clone_df 
    # TCRpublic.report() uses those instead of finding new ones.
    tp = TCRpublic(
        tcrrep = tr, 
        output_html_name = "quasi_public_clones_with_ndif.html")
    # Add any columns neighbor_diff columns 
    #that you want to display in the final report
    tp.labels.append('val_0')
    tp.labels.append('ct_0')
    tp.labels.append('pvalue')
    # chagne sort to be pvalue not publicity
    tp.sort_columns = ['pvalue']
    # because you are sorting by pvalue, change to True
    tp.sort_ascending = True
    tp.report()

