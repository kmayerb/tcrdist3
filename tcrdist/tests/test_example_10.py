import pytest 



def test_example_10():
    """
    If you just want a 'tcrdistances' of some target seqs against another set.

    (1) cell_df is asigned the first 10 cells in dash.csv
    (2) compute tcrdistances with default settings.
    (3) compute rectangular distance between clone_df and df2.
    (4) compute rectangular distance between clone_df and any 
    arbtirary df3, which need not be associated with the TCRrep object.
    (5) compute rectangular distance with only a subset of the TCRrep.clone_df
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    df2 = pd.read_csv("dash2.csv")
    df = df.head(10)                        #(1)
    tr = TCRrep(cell_df = df,               #(2)
                df2 = df2, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    
    assert tr.pw_alpha.shape == (10,10) 
    assert tr.pw_beta.shape  == (10,10)

    tr.compute_rect_distances()             # (3) 
    assert tr.rw_alpha.shape == (10,1924) 
    assert tr.rw_beta.shape  == (10,1924)

    df3 = df2.head(100)

    tr.compute_rect_distances(df = tr.clone_df, df2 = df3)  # (4) 
    assert tr.rw_alpha.shape == (10,100) 
    assert tr.rw_beta.shape  == (10,100)

    tr.compute_rect_distances(  df = tr.clone_df.iloc[0:2,], # (5)
                                df2 = df3)  
    assert tr.rw_alpha.shape == (2,100) 
    assert tr.rw_beta.shape  == (2,100)
