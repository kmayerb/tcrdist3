import pytest 

def test_example_10():
    """
    If you just want a 'tcrdistances' of some target seqs against another set.

    (1) cell_df is asigned the first 10 cells in dash.csv
    (2) TCRrep with default settings.
    (3) compute rectangular distance between clone_df and df2.
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv").head(10)   #(1)
    tr = TCRrep(cell_df = df,               #(2)
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                compute_distances = False)
    df2 = pd.read_csv("dash.csv")
    tr2 = TCRrep(cell_df = df2,               #(2)
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                compute_distances = False)
  
    tr.compute_rect_distances( df = tr.clone_df, df2 = tr2.clone_df) #(3)

    assert tr.rw_alpha.shape == (10,1920) 
    assert tr.rw_beta.shape  == (10,1920)
