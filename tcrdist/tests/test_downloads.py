import pytest 



def test_downloads():
    from tcrdist.setup_tests import list_available_zip_files  
    list_available_zip_files()
    """
    ['dash.zip',
    'human_T_alpha_beta_sim200K.zip',
    'vdjDB_PMID28636592.zip',
    'sant.csv.zip',
    'bulk.csv.zip',
    'wiraninha_sampler.zip',
    'ruggiero_mouse_sampler.zip',
    'ruggiero_human_sampler.zip',
    'britanova_human_beta_t_cb.tsv.sampler.tsv.zip',
    'emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip',
    'ruggiero_human_alpha_t.tsv.sampler.tsv.zip',
    'ruggiero_human_beta_t.tsv.sampler.tsv.zip']
    """
    from tcrdist.setup_tests import download_and_extract_zip_file
    """
    Get dash.zip contents
    """
    download_and_extract_zip_file( 'dash.zip', source = "dropbox", dest = ".")
    """
    Assert that file has been downloaded
    """
    import os
    assert os.path.isfile("dash.zip")
    """
    Assert that individual files of successively been inflated
    """
    assert os.path.isfile("dash.csv")
    assert os.path.isfile("dash2.csv")
    assert os.path.isfile("dash_human.csv")
    assert os.path.isfile("dash_beta_airr.csv")