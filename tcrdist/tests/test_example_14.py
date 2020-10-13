



def test_example_tree_args():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash.csv").sample(100, random_state=1).reset_index(drop = True)

    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    tcrtree = TCRtree(tcrrep = tr, 
          html_name = 'dash.mouse.b.tree.html')

    tcrtree.default_hcluster_diff_kwargs['x_cols'] = ['epitope']

    tcrtree.default_member_summ_kwargs['addl_cols'] : ['subject', 'epitope']

    tcrtree.default_plot_hclust_props['alpha_col'] = 'pvalue'
    tcrtree.default_plot_hclust_props['alpha'] = 1.0

    tcrtree.build_tree()