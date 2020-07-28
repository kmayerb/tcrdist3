import pytest
from tcrdist.repertoire import TCRrep
import pandas as pd
import numpy as np


"""
This testing modul is for tcrdist.repertoire
"""

"""
INTEGRATION TESTS
"""
def test_gamma_delta_all_steps():
    """
    Test that correctly formated gamma delta data can be 
    fully processed during TCRrep initialization
    """
    df = pd.read_csv("sant.csv")
    tr = TCRrep(
    cell_df = df, 
    organism = "human", 
    chains = ['gamma', "delta"],
    infer_cdrs     = True,
    deduplicate    = True,
    use_defaults   = True,
    store_all_cdr  = True,
    compute_distances = True,
    cpus    = 1, 
    db_file = 'alphabeta_gammadelta_db.tsv')
    assert isinstance(tr.pw_gamma, np.ndarray)
    assert isinstance(tr.pw_delta, np.ndarray)

def test_alpha_beta_all_steps():
    """
    Test that correctly formated alpha beta data can be 
    fully processed during TCRrep initialization
    """
    df = pd.read_csv("dash.csv")
    tr = TCRrep(
    cell_df = df, 
    organism = "mouse", 
    chains = ['alpha', "beta"],
    infer_cdrs     = True,
    deduplicate    = True,
    use_defaults   = True,
    store_all_cdr  = True,
    compute_distances = True,
    cpus    = 1, 
    db_file = 'alphabeta_gammadelta_db.tsv')

    assert isinstance(tr.pw_alpha, np.ndarray)
    assert isinstance(tr.pw_beta, np.ndarray)

def test_alpha_beta_all_steps_imgt_aligned_FALSE():
    """
    Test that correctly formated alpha beta data can be 
    fully processed during TCRrep initialization
    """
    df = pd.read_csv("dash.csv")
    tr = TCRrep(
    cell_df = df, 
    organism = "mouse", 
    chains = ['alpha', "beta"],
    infer_cdrs     = True,
    imgt_aligned   = False,
    deduplicate    = True,
    use_defaults   = True,
    store_all_cdr  = True,
    compute_distances = True,
    cpus    = 1, 
    db_file = 'alphabeta_gammadelta_db.tsv')

    assert isinstance(tr.pw_alpha, np.ndarray)
    assert isinstance(tr.pw_beta, np.ndarray)

def test_gamma_delta_manually_step_by_step():
    """
    Test that the user can go through the automatic 
    steps, manually, one-by-one for well formatted
    gamma delta data.
    """
    df = pd.read_csv("sant.csv")
    tr = TCRrep(
    cell_df = df, 
    organism = "human", 
    chains = ['gamma', "delta"],
    imgt_aligned   = False,
    infer_cdrs     = False,
    infer_index_cols = False,
    deduplicate    = False,
    use_defaults   = False,
    store_all_cdr  = False,
    compute_distances = False,
    cpus    = 1, 
    db_file = 'alphabeta_gammadelta_db.tsv')

    tr.infer_cdrs_from_v_gene(chain = "gamma")
    tr.infer_cdrs_from_v_gene(chain = "delta")
    tr.infer_index_cols()
    tr.show_incomplete()
    tr.deduplicate()
    tr._initialize_chain_specific_attributes()
    tr.stora_all_cdr = True
    tr.compute_distances()
    tr.pw_gamma

def test_alpha_beta_manually_step_by_step():
    """
    Test that the user can go through the automatic 
    steps, manually, one-by-one for well formatted
    alpha beta data.
    """
    df = pd.read_csv("dash.csv")
    tr = TCRrep(
    cell_df = df, 
    organism = "mouse", 
    chains = ['alpha', "beta"],
    imgt_aligned   = False,
    infer_cdrs     = False,
    infer_index_cols = False,
    deduplicate    = False,
    use_defaults   = False,
    store_all_cdr  = False,
    compute_distances = False,
    cpus    = 1, 
    db_file = 'alphabeta_gammadelta_db.tsv')

    tr.infer_cdrs_from_v_gene(chain = "alpha")
    tr.infer_cdrs_from_v_gene(chain = "beta")
    tr.infer_index_cols()
    
    tr.deduplicate()
    tr._initialize_chain_specific_attributes()
    tr.stora_all_cdr = True
    tr.compute_distances()
    tr.pw_beta
   

"""
NULL INPUTS - Test that object can be created without any real inputs
"""
def test_null_import():
    """
    Null Import - Test that we can import with no inputs and everything set to false
    """
    with pytest.warns(UserWarning):
        tr = TCRrep(
            organism          = "mouse",
            chains            = ["alpha", "beta"],
            cell_df           = None,
            clone_df          = None,
            imgt_aligned      = False,
            infer_cdrs        = False,
            infer_index_cols  = False,
            deduplicate       = False,
            use_defaults      = False,
            store_all_cdr     = False,
            compute_distances = False,
            cpus              = 1, 
            db_file           = 'alphabeta_gammadelta_db.tsv')
    assert isinstance(tr.chains, list) 

def test_null_import_defaults_left_blank():
    with pytest.warns(UserWarning):
        tr = TCRrep(
            cell_df           = None,
            clone_df          = None,
            imgt_aligned      = False,
            infer_cdrs        = False,
            infer_index_cols  = False,
            deduplicate       = False,
            use_defaults      = False,
            store_all_cdr     = False,
            compute_distances = False)
    assert isinstance(tr.chains, list) 


def test_null_import_gamma_delta():
    with pytest.warns(UserWarning):
        tr = TCRrep(
            organism          = "human",
            chains            = ["gamma", "delta"],
            cell_df           = None,
            clone_df          = None,
            imgt_aligned      = False,
            infer_cdrs        = False,
            infer_index_cols  = False,
            deduplicate       = False,
            use_defaults      = False,
            store_all_cdr     = False,
            compute_distances = False,
            cpus              = 1, 
            db_file           = 'alphabeta_gammadelta_db.tsv')
    assert isinstance(tr.chains, list) 



"""
VALIDATION TESTS
"""
def test_validate_organism():
    """Test that incorrect organisms raises ValueError"""
    with pytest.raises(ValueError) as info:
        tr = TCRrep(
            organism          = "hamster",
            chains            = ["alpha", "beta"])
    assert str(info.value) == "organism must be 'mouse' or 'human'"

def test_validate_chains():
    """Test that incorrect chain raise causes ValueError"""
    with pytest.raises(ValueError) as info:
        tr = TCRrep(
            organism          = "mouse",
            chains            = ["ALPHA", "beta"])
    assert str(info.value) == "TCRrep chains arg can be one or more of the following ['alpha', 'beta', 'gamma', 'delta'] case-sensitive"


def test_validate_imgt_aligned():
    """Test that incorrect chain raise causes ValueError"""
    with pytest.raises(ValueError) as info:
        tr = TCRrep(
            organism          = "mouse",
            chains            = ["alpha", "beta"], 
            imgt_aligned      = "align")
    assert str(info.value) == "TCRrep imgt_aligned argument must be a boolean"

def test_validate_cell_df_is_pandas_dataframe():
    with pytest.warns(UserWarning):
        tr = TCRrep(
            cell_df           = np.zeros(1),
            organism          = "mouse",
            chains            = ["alpha", "beta"],
            clone_df          = None,
            imgt_aligned      = False,
            infer_cdrs        = False,
            infer_index_cols  = False,
            deduplicate       = False,
            use_defaults      = False,
            store_all_cdr     = False,
            compute_distances = False,
            cpus              = 1, 
            db_file           = 'alphabeta_gammadelta_db.tsv')

def test_warn_on_bad_db_file():
    """Test that incorrect chain raise causes ValueError"""
    with pytest.warns(UserWarning):
        tr = TCRrep(
            organism          = "human",
            chains            = ["gamma", "delta"],
            cell_df           = None,
            clone_df          = None,
            imgt_aligned      = False,
            infer_all_genes   = False,
            infer_cdrs        = False,
            infer_index_cols  = False,
            deduplicate       = False,
            use_defaults      = False,
            store_all_cdr     = False,
            compute_distances = False,
            cpus              = 1, 
            db_file           = 'WRONGalphabeta_gammadelta_db.tsv')