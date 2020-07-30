from zipfile import ZipFile
import requests
from tcrdist import paths

__all__ = ['download_and_extract_zip_file']

"""
python -c "from tcrdist.setup_tests import *; download_and_extract_zip_file('bulk.csv.zip')"
python -c "from tcrdist.setup_tests import *; download_and_extract_zip_file('bulk.csv.zip')"
"""

# public facing data files, L looksup url based for dropbox and aws
L = {"dash.zip":
        {'dropbox': {
            'url' : "https://www.dropbox.com/s/pce3f9816ntzjki/dash.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "sant.csv.zip":
        {'dropbox': {
            'url' : "https://www.dropbox.com/s/8p3djrdd270ad0n/sant.csv.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "bulk.csv.zip":
        {'dropbox': {
            'url' : "https://www.dropbox.com/s/g6k2h1ed5d5sabz/bulk.csv.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "wiraninha_sampler.zip":
        {'dropbox': {
            'url' : "https://www.dropbox.com/s/ily0td3tn1uc7bi/wiraninha_sampler.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "ruggiero_mouse_sampler.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/yz8v1c1gf2eyzxk/ruggiero_mouse_sampler.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "ruggiero_human_sampler.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/jda6qtemk65zlfk/ruggiero_human_sampler.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "britanova_human_beta_t_cb.tsv.sampler.tsv.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/87n5v2by80xhy1q/britanova_human_beta_t_cb.tsv.sampler.tsv.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/04mxrzw7f5wkg1x/emerson_human_beta_t_cmvneg.tsv.sampler.tsv.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "ruggiero_human_alpha_t.tsv.sampler.tsv.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/9h84bzhd0asfym7/ruggiero_human_alpha_t.tsv.sampler.tsv.zip?dl=1"},
        'aws': { 
            'url' : None}
    },
    "ruggiero_human_beta_t.tsv.sampler.tsv.zip":
        {'dropbox':{
            'url' : "https://www.dropbox.com/s/onr5lntmlm4fivi/ruggiero_human_beta_t.tsv.sampler.tsv.zip?dl=1"},
        'aws': { 
            'url' : None}
        }
    }

def get_url(zipfile, source = 'dropbox'):
    url = L[zipfile][source]['url']
    return url

def download_and_extract_zip_file(zipfile, source = 'dropbox', dest = paths.path_to_base):
    """ Downloads and extract a zip file to destinateion folder """
    url = get_url(zipfile, source = source)
    r = requests.get(url) 
    with open(zipfile,'wb') as f:
        f.write(r.content)
    with ZipFile(zipfile, 'r') as zipObj:
        # Extract all the contents of zip file in different directory
        zipObj.extractall(dest)

