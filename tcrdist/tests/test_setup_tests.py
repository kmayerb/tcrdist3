import pytest 
from tcrdist import setup_tests #get_url, download_and_extract_zip_file
from tcrdist import paths
import os

def test_get_url():
    url = setup_tests.get_url("bulk.csv.zip", source = 'dropbox')
    assert url == "https://www.dropbox.com/s/g6k2h1ed5d5sabz/bulk.csv.zip?dl=1"

def test_download_and_extract_zip_file():
    setup_tests.download_and_extract_zip_file(zipfile = "sant.csv.zip", source = 'dropbox', dest = paths.path_to_base)
    # assert download succeed
    assert os.path.isfile(os.path.join(paths.path_to_base, "sant.csv.zip"))
    # assert download + unzip succeed
    assert os.path.isfile(os.path.join(paths.path_to_base, "sant.csv"))
    
