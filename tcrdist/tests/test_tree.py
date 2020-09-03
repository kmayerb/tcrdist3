"""
tcrdist3 hierahical clustering trees alpha/beta fill tests including
generating SVGs these tests take some time.
"""

def test_init_tree():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.mouse.ab.tree.html')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_size_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_unique_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('percent_missing_olga_beta')
    tcrtree.build_tree()

    assert os.path.isfile('dash.mouse.ab.tree.html')

def test_init_tree_beta():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.mouse.b.tree.html')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_size_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_unique_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('percent_missing_olga_beta')
    tcrtree.build_tree()
    assert os.path.isfile('dash.mouse.b.tree.html')

def test_init_tree_alpha():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.mouse.a.tree.html')



    tcrtree.build_tree()
    assert os.path.isfile('dash.mouse.a.tree.html')


def test_init_tree_human_beta():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash_human.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.human.b.tree.html')
    tcrtree.combine_olga = True
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_size_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_unique_olga_beta')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('percent_missing_olga_beta')
    tcrtree.build_tree( )
    assert os.path.isfile('dash.human.b.tree.html')


def test_init_tree_human_alpha():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash_human.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.human.a.tree.html')  
    tcrtree.combine_olga = True
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_size_olga_alpha')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('ref_unique_olga_alpha')
    tcrtree.default_plot_hclust_props['tooltip_cols'].append('percent_missing_olga_alpha')
    tcrtree.build_tree( )

    assert os.path.isfile('dash.human.a.tree.html')


def test_init_tree_human_alpha_beta():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash_human.csv").sample(10).reset_index(drop = True)
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['alpha', 'beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.human.ab.tree.html')
    tcrtree.build_tree()
    assert os.path.isfile('dash.human.ab.tree.html')


