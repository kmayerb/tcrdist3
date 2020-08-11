import pytest 



def test_introduction_1():
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False)