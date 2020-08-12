import pytest



def test_introduction_6():
    """
    Basic Specificity Neighborhoods based on a Hierarchical Clustering
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta','alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    
    from tcrdist.rep_diff import hcluster_diff, member_summ
    from hierdiff import plot_hclust_props

    # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
    tr.clone_df['PA'] = ['PA' if x == 'PA' else 'X' for x in tr.clone_df.epitope]
    
    res, Z= hcluster_diff(tr.clone_df, tr.pw_beta, x_cols = ['PA'], count_col = 'count')
    
    res_summary = member_summ(res_df = res, clone_df = tr.clone_df, addl_cols=['epitope'])
    
    res_detailed = pd.concat([res, res_summary], axis = 1)
    
    html = plot_hclust_props(Z,
                title='PA Epitope Example',
                res=res_detailed,
                tooltip_cols=['cdr3_b_aa','v_b_gene', 'j_b_gene','epitope'],
                alpha=0.00001, colors = ['blue','gray'],
                alpha_col='pvalue')

    with open('hierdiff_example.html', 'w') as fh:
        fh.write(html)