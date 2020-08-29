import pytest
import pandas as pd
import numpy as np
import os

def test__default_tcrsampler_human_beta():
    from tcrdist.automate import _default_tcrsampler_human_beta
    from tcrsampler.sampler import TCRsampler
    r = _default_tcrsampler_human_beta()
    assert isinstance(r, TCRsampler)

def test__default_tcrsampler_human_alpha():
    from tcrdist.automate import _default_tcrsampler_human_alpha
    from tcrsampler.sampler import TCRsampler
    r = _default_tcrsampler_human_alpha()
    assert isinstance(r, TCRsampler)

def test__default_tcrsampler_mouse_beta():
    from tcrdist.automate import _default_tcrsampler_mouse_beta
    from tcrsampler.sampler import TCRsampler
    r = _default_tcrsampler_mouse_beta()
    assert isinstance(r, TCRsampler)

def test__default_tcrsampler_mouse_alpha():
    from tcrdist.automate import _default_tcrsampler_mouse_alpha
    from tcrsampler.sampler import TCRsampler
    r = _default_tcrsampler_mouse_alpha()
    assert isinstance(r, TCRsampler)

"""
Test Automatic Pgen Esimator Functions
"""
def test__auto_pgen_mouse_beta():
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv").sample(10, random_state = 3)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    from tcrdist.automate import _auto_pgen
    tr = _auto_pgen(tr, chain = 'beta', organism = 'mouse')
    assert isinstance(tr.clone_df.pgen_cdr3_b_aa, pd.Series)

    
def test_auto_pgen_mouse_beta():
    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv").sample(10, random_state = 3)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    from tcrdist.automate import auto_pgen
    tr = auto_pgen(tr)
    assert isinstance(tr.clone_df.pgen_cdr3_b_aa, pd.Series)

def test_auto_pgen_mouse_alpha_beta_ValueError():
    """ If auto_pgen called on TCRrep with organism == mouse and 
    chain including 'alpha', Raises ValueError"""

    import pandas as pd
    from tcrdist.repertoire import TCRrep

    df = pd.read_csv("dash.csv").sample(10, random_state = 3)
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    from tcrdist.automate import auto_pgen
    with pytest.raises(ValueError):
        tr = auto_pgen(tr)
        
def test__auto_pgen_human_alpha_beta():
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    df = pd.read_csv("dash_human.csv").sample(10, random_state = 3)
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
   
    from tcrdist.automate import _auto_pgen
    tr = _auto_pgen(tr, chain = 'beta', organism = 'human')
    assert isinstance(tr.clone_df.pgen_cdr3_b_aa, pd.Series)
    assert isinstance(tr.clone_df.pgen_cdr3_b_aa[0], float)
    tr = _auto_pgen(tr, chain = 'alpha', organism = 'human')
    assert isinstance(tr.clone_df.pgen_cdr3_a_aa, pd.Series)
    assert isinstance(tr.clone_df.pgen_cdr3_a_aa[0], float)

def test_auto_pgen_human_alpha_beta():
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    df = pd.read_csv("dash_human.csv").sample(10, random_state = 3)
    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')
    from tcrdist.automate import auto_pgen
    tr = auto_pgen(tr)
    assert isinstance(tr.clone_df.pgen_cdr3_b_aa, pd.Series)



