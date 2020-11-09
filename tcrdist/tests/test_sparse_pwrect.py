import pytest

def test_example_10_sparse():
    """
    If you just want a 'tcrdistances' of some target seqs against another set computed sparsely

    (1) cell_df is asigned the first 10 cells in dash.csv
    (2) compute tcrdistances with default settings.
    (3) compute rectangular distance between clone_df and df2.
    (4) compute rectangular distance between clone_df and any 
    arbtirary df3, which need not be associated with the TCRrep object.
    (5) compute rectangular distance with only a subset of the TCRrep.clone_df
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import pw2dense
    import numpy as np

    df = pd.read_csv("dash.csv")
    df2 = pd.read_csv("dash2.csv")
    df = df.head(10)                        #(1)
    tr = TCRrep(cell_df = df,               #(2)
                df2 = df2, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    tr.compute_rect_distances(df = tr.clone_df, df2 = df2)
    assert tr.rw_alpha.shape == (10,1924) 
    assert tr.rw_beta.shape  == (10,1924)

    rw_alpha = tr.rw_alpha.copy()
    rw_beta = tr.rw_beta.copy()

    radius = 100
    tr.cpus = 1
    tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = df2, radius = radius)
    d = pw2dense(tr.rw_alpha, radius)
    print(rw_alpha[:2, :10])
    print(tr.rw_alpha.todense()[:2, :10])
    print(d[:2, :10])
    assert np.all(rw_alpha[rw_alpha <= radius] == d[d <= radius])

    d = pw2dense(tr.rw_beta, radius)
    assert np.all(rw_beta[rw_beta <= radius] == d[d <= radius])

    
    radius = 5000
    tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = df2, radius = radius)
    d = pw2dense(tr.rw_alpha, radius)
    print(rw_alpha[:2, :10])
    print(tr.rw_alpha.todense()[:2, :10])
    assert np.all(rw_alpha == d)

    d = pw2dense(tr.rw_beta, radius)
    assert np.all(rw_beta == d)
    

def test_example_10_sparse_multiprocessing():
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import pw2dense
    import numpy as np

    df2 = pd.read_csv("dash2.csv")
    tr = TCRrep(cell_df = df2,               #(2)
                df2 = df2, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    tr.compute_rect_distances(df = tr.clone_df, df2 = df2)
    assert tr.rw_alpha.shape == (1924,1924) 
    assert tr.rw_beta.shape  == (1924,1924)

    rw_alpha = tr.rw_alpha.copy()
    rw_beta = tr.rw_beta.copy()

    radius = 150
    tr.cpus = 2
    tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = df2, radius = radius)
    d = pw2dense(tr.rw_alpha, radius)
    assert np.all(rw_alpha[rw_alpha <= radius] == d[d <= radius])

    d = pw2dense(tr.rw_beta, radius)
    assert np.all(rw_beta[rw_beta <= radius] == d[d <= radius])