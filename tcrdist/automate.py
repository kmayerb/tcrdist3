"""
Functions that automate expressive functions of tcrdist3. 

These use sensible defaults but users should consider 
flexibile elements of these functions that might better 
suit their specific data or research questions. 
"""
import os
import pandas as pd
import numpy as np
import warnings
from progress.bar import IncrementalBar

"""
Test Automatic Pgen Esimator Functions
"""
def auto_pgen(tcrrep):
	"""
	Automate a pgen estimation of cdr3s alpha/beta given a tcrrep with a clones_df attribute

	Parameters
	----------
	tcrrep : tcrdist.repertoire.TCRrep
		TCRrep instance with a clone_df

	Returns
	-------
	tcrrep : tcrdist.repertoire.TCRrep

	tcrrep : tcrdist.repertoire.TCRrep
		TCRrep instance with a clone_df
	"""
	if tcrrep.organism == "mouse" and 'alpha' in tcrrep.chains:
		raise ValueError("UNFORTUNATELY OLGA DOES NOT YET SUPPORT MOUSE ALPHA PGEN ESTIMATES")

	for chain in tcrrep.chains:
		tcrrep = _auto_pgen(tcrrep=tcrrep, organism = tcrrep.organism, chain = chain)
	return tcrrep


def _auto_pgen(tcrrep = None, organism = 'human', chain = 'beta'):
	"""
	Automate a pgen estimation of cdr3s alpha/beta given a tcrrep with a clones_df attribute

	Parameters
	----------
	tcrrep : tcrdist.repertoire.TCRrep
		TCRrep instance with a clone_df
	organism : str
		'human' or 'mouse'
	chain : str
		'beta' or 'alpha'

	Returns
	-------
	tcrrep : tcrdist.repertoire.TCRrep
	"""
	import tcrdist
	import parmap
	import pandas as pd
	from tcrdist.pgen import OlgaModel
	assert organism in ['human', 'mouse']
	assert chain in ['beta', 'alpha']
	assert isinstance(tcrrep, tcrdist.repertoire.TCRrep)
	assert isinstance(tcrrep.clone_df, pd.DataFrame)

	cdr3_col = {'alpha': 'cdr3_a_aa', 'beta': 'cdr3_b_aa'}[chain]
	cdr3s = tcrrep.clone_df[cdr3_col]

	olga_models = { 
		('human', 'beta') : 
			OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ"),
		('human', 'alpha') :
			OlgaModel(chain_folder = "human_T_alpha", recomb_type="VJ"),
		('mouse', 'beta') : 
			OlgaModel(chain_folder = "mouse_T_beta", recomb_type="VDJ")}
	
	olga_model = olga_models[(organism, chain)]

	pgens = parmap.map(olga_model.compute_aa_cdr3_pgen, cdr3s, pm_pbar = True)
	tcrrep.clone_df[f"pgen_{cdr3_col}"] = pgens

	return tcrrep
	


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
