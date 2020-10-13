



def test_example_tree_args():
    import os
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.tree import TCRtree

    df = pd.read_csv("dash.csv").sample(100).reset_index(drop = True)

    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    tcrtree = TCRtree(tcrrep = tr, html_name = 'dash.mouse.ab.tree.html')

    tcrtree.build_tree()

    assert os.path.isfile('dash.mouse.ab.tree.html')