



def test_neighbors_and_publicity_directly():
    """
    Instead of enforcing a fixed radius, 
    use a radius specific to each
    centroid, specified in an additional 
    column.
    """
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic  
    fn = os.path.join('tcrdist','data','covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.radius.csv')
    df = pd.read_csv(fn)
    tr = TCRrep(cell_df = df[['cohort','subject','v_b_gene', 'j_b_gene','cdr3_b_aa', 'radius']], 
                organism = "human", 
                chains = ["beta"])
    
    # NEIGHBORS
    from tcrdist.public import _neighbors_fixed_radius
    from tcrdist.public import _neighbors_variable_radius
    # returns lists of lists of all neighbors at fixed of variable radii
    _neighbors_fixed_radius(pwmat = tr.pw_beta, radius = 18)
    _neighbors_variable_radius(pwmat = tr.pw_beta, radius_list = tr.clone_df.radius)

    # returns the number (K) neighbors at fixed or vriable radii
    from tcrdist.public import _K_neighbors_fixed_radius
    from tcrdist.public import _K_neighbors_variable_radius
    _K_neighbors_fixed_radius(pwmat = tr.pw_beta, radius = 18)
    _K_neighbors_variable_radius(pwmat = tr.pw_beta, radius_list = tr.clone_df.radius)

    # First find neighbors by your favorite method 
    tr.clone_df['neighbors'] = _neighbors_variable_radius(
        pwmat = tr.pw_beta, 
        radius_list = tr.clone_df.radius)
    # Once neighbors are added to a clone_df you can easily determine publicity. 
    tr.clone_df['nsubject']   = tr.clone_df['neighbors'].\
        apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
    tr.clone_df['qpublic']   = tr.clone_df['nsubject'].\
        apply(lambda x: x > 1)
