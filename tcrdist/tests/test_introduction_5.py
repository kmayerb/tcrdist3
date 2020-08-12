import pytest



def test_introduction_5():
    """
    Basic Specificity Neighborhoods based on a Fixed Radius
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta','alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    
    from tcrdist.rep_diff import neighborhood_diff
    # diff testing is pasted on binary comparison, so all epitope not 'PA' are set to 'X'
    tr.clone_df['PA'] = ['PA' if x == 'PA' else 'X' for x in tr.clone_df.epitope]
    
    # Larger Radius
    tr.nn_df = neighborhood_diff(clone_df= tr.clone_df, 
        pwmat = tr.pw_beta + tr.pw_alpha, 
        count_col = 'count', 
        x_cols = ['PA'], 
        knn_radius = 150)

    tr.nn_df[['K_neighbors', 'val_0', 'ct_0', 'val_2', 'ct_2', 'RR','OR', 'pvalue', 'FWERp','FDRq']].\
        sort_values(['FDRq']).sort_values(['OR','ct_0'], ascending = False)

    # Smaller Radius
    tr.nn_df = neighborhood_diff(clone_df= tr.clone_df, 
        pwmat = tr.pw_beta + tr.pw_alpha, 
        count_col = 'count', 
        x_cols = ['PA'], 
        knn_radius = 75)

    tr.nn_df[['K_neighbors', 'val_0', 'ct_0', 'val_2', 'ct_2', 'RR','OR', 'pvalue', 'FWERp','FDRq']].\
        sort_values(['FDRq']).sort_values(['OR','ct_0'], ascending = False)