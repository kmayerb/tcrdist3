



def test_example_tree():
    """
    An example showing how to create an interactive
    tree from a sample of mouse TCRs 
    """
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

    tcrtree.build_tree()

    assert os.path.isfile('dash.mouse.b.tree.html')