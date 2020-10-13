



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

    tcrtree.default_hcluster_diff_kwargs = \
            {'clone_df': None,
             'pwmat': None,
             'x_cols': ['epitope'],
             'Z': None,
             'count_col': 'count',
             'subset_ind': None,
             'hclust_method': 'complete',
             'optimal_ordering': True,
             'test_method': 'fishers'}

    tcrtree.default_member_summ_kwargs = \
            {'key_col': 'neighbors_i',
            'count_col': 'count',
            'addl_cols': ['subject'],
            'addl_n': 1}

    tcrtree.default_plot_hclust_props = \
            {'title': '',
            'alpha_col': 'pvalue',
            'alpha': 0.05,
            'tooltip_cols': ['subject',
            'mean_dist',
            'pct_dist_75',
            'pct_dist_50',
            'pct_dist_25',
            'fuzzy_simpson_diversity_75',
            'fuzzy_simpson_diversity_50',
            'fuzzy_simpson_diversity_25',
            'cdr3_b_aa',
            'v_b_gene',
            'j_b_gene',
            'svg_beta',
            'svg_raw_beta',
            'ref_size_beta',
            'ref_unique_beta',
            'percent_missing_beta']}

    tcrtree.build_tree()

