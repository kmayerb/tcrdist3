import pytest



def test_introduction_4():
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    """
    tr.pw_beta

    array([[  0., 162.,  90., ..., 120.,  75.,  75.],
           [162.,   0., 123., ...,  96., 126., 138.],
           [ 90., 123.,   0., ...,  81., 120., 120.],
           ...,
           [120.,  96.,  81., ...,   0.,  69.,  69.],
           [ 75., 126., 120., ...,  69.,   0.,  12.],
           [ 75., 138., 120., ...,  69.,  12.,   0.]])
    """