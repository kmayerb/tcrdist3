import pandas as pd 
from tcrdist.repertoire import TCRrep


def test_example_11():
    """
    Example of comparing two clone_df from separate repertoires of different sizes
    
    1. Load a bulk beta-chain dataset with > 70 Clones.
    2. Prepare a beta chain clones_df DataFrame with cdr3_b_aa, pmhc_b_aa, cdr2_b_aa, cdr1_b_aa columns.
    3. Extract copy of processed dataframe.
    4. Load data of interest into a TCRrep instance.
    5. Compare the first 100 clones in primary dataset (100 Clones) against the the entire bulk dataset (77,538 Clones)
    """
    df_bulk = pd.read_csv('bulk.csv', sep = ",")        #(1)
    tr_bulk = TCRrep(cell_df = df_bulk,                 #(2)
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False)
    df_search = tr_bulk.clone_df.copy()                 #(3)

    df = pd.read_csv("dash.csv").head(100)              

    tr = TCRrep(cell_df = df,                           #(4)
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    assert tr.pw_beta.shape == (100,100)
    assert df_search.shape  == (77538, 8)
    tr.compute_rect_distances(df = tr.clone_df, df2 = df_search) #(5)
    assert tr.rw_cdr3_b_aa.shape == (100, 77538)
