import pytest

def test_alpha_beta():
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
    from tcrdist.rep_funcs import  compute_n_tally_out_of_memory2
    from hierdiff.association_testing import cluster_association_test

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df.sample(100, random_state = 1), 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv', 
                compute_distances = True,
                store_all_cdr = False)

    check_beta = tr.pw_beta.copy(); check_beta[check_beta == 0] = 1
    check_alpha = tr.pw_alpha.copy(); check_alpha[check_alpha == 0] = 1
    check_alpha_beta = check_beta + check_alpha
    

    S, fragments = compute_pw_sparse_out_of_memory2(    tr = tr,
                                                        row_size      = 50,
                                                        pm_processes  = 1,
                                                        pm_pbar       = True,
                                                        max_distance  = 1000,
                                                        reassemble    = True,
                                                        cleanup       = False,
                                                        assign        = True)
    
    assert np.all(tr.pw_beta == check_beta)
    assert np.all(tr.pw_alpha == check_alpha)

    ndif1 = compute_n_tally_out_of_memory2(fragments, 
                                         to_file = False, 
                                         to_memory = True,
                                         pm_processes = 2, 
                                         x_cols = ['epitope'],
                                         count_col='count',
                                         knn_neighbors= None,
                                         knn_radius =100)

    from hierdiff.association_testing import cluster_association_test
    ndif1 = cluster_association_test(res = ndif1, y_col='cmember', method='chi2')


    from tcrdist.rep_diff import neighborhood_diff
    ndif2 = neighborhood_diff(clone_df= tr.clone_df, 
        pwmat = np.array(tr.pw_beta.todense() + tr.pw_alpha.todense()),
        count_col = 'count', 
        x_cols = ['epitope'], 
        knn_radius = 100, 
        test_method = "chi2")

    assert ndif1.shape == ndif2.shape
    np.all(ndif2['FDRq'].to_list() == ndif2['FDRq'].to_list())





def test_beta():
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
    from tcrdist.rep_funcs import  compute_n_tally_out_of_memory2
    from hierdiff.association_testing import cluster_association_test

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df.sample(100, random_state = 1), 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv', 
                compute_distances = True,
                store_all_cdr = False)

    check_beta = tr.pw_beta.copy(); check_beta[check_beta == 0] = 1

    S, fragments = compute_pw_sparse_out_of_memory2(    tr = tr,
                                                        row_size      = 50,
                                                        pm_processes  = 2,
                                                        pm_pbar       = True,
                                                        max_distance  = 1000,
                                                        reassemble    = True,
                                                        cleanup       = False,
                                                        assign        = True)
    import numpy as np
    assert np.all(tr.pw_beta == check_beta)
    assert 'pw_alpha' not in tr.__dict__.keys()


    ndif1 = compute_n_tally_out_of_memory2(fragments, 
                                         to_file = False, 
                                         to_memory = True,
                                         pm_processes = 2, 
                                         x_cols = ['epitope'],
                                         count_col='count',
                                         knn_neighbors= None,
                                         knn_radius =100)

    from hierdiff.association_testing import cluster_association_test
    ndif1 = cluster_association_test(res = ndif1, y_col='cmember', method='chi2')


    from tcrdist.rep_diff import neighborhood_diff
    ndif2 = neighborhood_diff(clone_df= tr.clone_df, 
        pwmat = np.array(tr.pw_beta.todense()),
        count_col = 'count', 
        x_cols = ['epitope'], 
        knn_radius = 100, 
        test_method = "chi2")

    assert ndif1.shape == ndif2.shape
    np.all(ndif2['FDRq'].to_list() == ndif2['FDRq'].to_list())



def test_alpha():
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
    from tcrdist.rep_funcs import  compute_n_tally_out_of_memory2
    from hierdiff.association_testing import cluster_association_test

    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df.sample(100, random_state = 1), 
                organism = 'mouse', 
                chains = ['alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv', 
                compute_distances = True,
                store_all_cdr = False)

    check_alpha = tr.pw_alpha.copy(); check_alpha[check_alpha == 0] = 1

    S, fragments = compute_pw_sparse_out_of_memory2(    tr = tr,
                                                        row_size      = 50,
                                                        pm_processes  = 2,
                                                        pm_pbar       = True,
                                                        max_distance  = 1000,
                                                        reassemble    = True,
                                                        cleanup       = False,
                                                        assign        = True)
    import numpy as np
    assert np.all(tr.pw_alpha == check_alpha)
    assert 'pw_beta' not in tr.__dict__.keys()

    ndif1 = compute_n_tally_out_of_memory2(fragments, 
                                         to_file = False, 
                                         to_memory = True,
                                         pm_processes = 2, 
                                         x_cols = ['epitope'],
                                         count_col='count',
                                         knn_neighbors= None,
                                         knn_radius =100)

    from hierdiff.association_testing import cluster_association_test
    ndif1 = cluster_association_test(res = ndif1, y_col='cmember', method='chi2')


    from tcrdist.rep_diff import neighborhood_diff
    ndif2 = neighborhood_diff(clone_df= tr.clone_df, 
        pwmat = np.array(tr.pw_alpha.todense()),
        count_col = 'count', 
        x_cols = ['epitope'], 
        knn_radius = 100, 
        test_method = "chi2")

    assert ndif1.shape == ndif2.shape
    np.all(ndif2['FDRq'].to_list() == ndif2['FDRq'].to_list())

