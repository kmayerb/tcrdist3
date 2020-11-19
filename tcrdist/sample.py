import os
import re
import pandas as pd
import numpy as np
import warnings
from tcrsampler.sampler import TCRsampler

def _default_sampler(organism = 'human', chain = 'beta'):
    assert organism in ['human', 'mouse']
    assert chain in ['beta','alpha']

    default_tcrsampler_generator = {
        ('human','beta'): 
            _default_tcrsampler_human_beta,
        ('human','alpha'): 
            _default_tcrsampler_human_alpha,
        ('mouse','beta'): 
            _default_tcrsampler_mouse_beta,
        ('mouse','alpha'): 
            _default_tcrsampler_mouse_alpha, 
        }[(organism, chain)]
    
    return default_tcrsampler_generator

def _default_sampler_olga(organism = 'human', chain = 'beta'):
    assert organism in ['human', 'mouse']
    assert chain in ['beta','alpha']
    if organism == "mouse" and chain == "alpha":
        raise KeyError("No currenlty available mouse-alpha chain background form olga")

    default_tcrsampler_generator = {
        ('human','beta'): 
            _default_tcrsampler_olga_human_beta,
        ('human','alpha'): 
            _default_tcrsampler_olga_human_alpha,
        ('mouse','beta'): 
            _default_tcrsampler_olga_mouse_beta
        }[(organism, chain)]
    
    return default_tcrsampler_generator


def _default_tcrsampler_olga_human_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_human_beta_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_olga_human_alpha(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_human_alpha_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_olga_mouse_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background =  'olga_mouse_beta_t.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='olga_sampler.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t


def _default_tcrsampler_human_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default human beta sampler 'britanova_human_beta_t_cb.tsv.sampler.tsv'

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv'
        
    if default_background_if_missing is None:
        default_background_if_missing ='britanova_human_beta_t_cb.tsv.sampler.tsv.zip'
    
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_human_alpha(default_background = None, default_background_if_missing=None ):
    """
    Responsible for providing the default human alpha sampler 'ruggiero_human_alpha_t.tsv.sampler.tsv'
    """
    from tcrsampler.sampler import TCRsampler
    if default_background is None:
        default_background = 'ruggiero_human_alpha_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_human_alpha_t.tsv.sampler.tsv.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_mouse_beta(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default mouse beta sampler

    Returns
    -------
    t : tcrsampler.sampler.TCRsampler 
    """
    from tcrsampler.sampler import TCRsampler

    if default_background is None:
        default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_mouse_sampler.zip'

    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def _default_tcrsampler_mouse_alpha(default_background = None, default_background_if_missing=None):
    """
    Responsible for providing the default mouse alpha sampler
    """
    from tcrsampler.sampler import TCRsampler
    
    if default_background is None:
        default_background = 'ruggiero_mouse_alpha_t.tsv.sampler.tsv'
    if default_background_if_missing is None:
        default_background_if_missing =  'ruggiero_mouse_sampler.zip'
    
    print(default_background)

    try: 
        t = TCRsampler(default_background=default_background)
    except OSError:
        t = TCRsampler()
        t.download_background_file(default_background_if_missing)
        t = TCRsampler(default_background=default_background)
    return t

def allele_01(genename):
    """
    >>> allele_01('TRBV19*02')
    'TRBV19*01'
    """
    g,a = re.match(pattern = '(.*)([0-9])$', string= genename).groups()
    return(f"{g}1")