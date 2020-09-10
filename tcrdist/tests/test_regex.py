from tcrdist.regex import _index_to_regex_str, _index_to_seqs, _matrix_to_regex, _index_to_matrix, _list_to_regex_component
import pytest
import pandas as pd
import os
from tcrdist.repertoire import TCRrep	
from tcrdist.regex import _index_to_matrix

"""
Fixture
"""
fn = os.path.join('dash.csv')
df = pd.read_csv(fn , sep = ",")
df = df[(df.epitope == 'NP')].reset_index()
tr = TCRrep( cell_df = df,
        organism = 'mouse',
        chains = ['beta'],
        db_file = 'alphabeta_gammadelta_db.tsv')
"""
unit and coverage tests for the /tcrdist/regex.py 
"""
def test_list_to_regex_component_sinlge():
	assert _list_to_regex_component(['A'],max_ambiguity = 3) == 'A'

def test_list_to_regex_component_group():
	assert _list_to_regex_component(['A','T','C'],max_ambiguity = 3) == '[ATC]'

def test_list_to_regex_component_group_more_than_max_ambiguity():
	assert _list_to_regex_component(['A','T','C','G'],max_ambiguity = 3) == '.'

def test_list_to_regex_component_change_ambiguity_threshold():
	assert _list_to_regex_component(['A','T','C','G'],max_ambiguity = 4) == '[ATCG]'

def test_list_to_regex_component_gap_makes_group_optional():
	assert _list_to_regex_component(['A','-','C']) == '[AC]?'

def test__matrix_to_regex():
	import numpy as np
	import pandas as pd
	arr = np.array([   [ 0. ,  3.3,  0. ,  0. ,  0.4,  0. ,  0. ,  0. ],
					   [ 0. ,  0. ,  0. ,  0. ,  0.  ,  0. ,  3.6,  4.5],
					   [ 0. ,  0. ,  0. ,  0. ,  0.  ,  0. ,  0. ,  0. ],
					   [ 0. ,  0. ,  0. ,  0. ,  0.  ,  4.5,  0. ,  0. ],
					   [ 4.5,  0. ,  0.4,  0. ,  0.  ,  0. ,  0. ,  0. ],
					   [ 0. ,  0. ,  0. ,  0. ,  2.67,  0. ,  0. ,  0. ],
					   [ 0. ,  0. ,  0. ,  0. ,  0.  ,  0. ,  0. ,  0. ],
					   [ 0. ,  0. ,  0. ,  0.1,  0.  ,  0. ,  0. ,  0. ]])
	matrix = pd.DataFrame(arr, 
		index = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G'], 
		columns = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G'])
	assert _matrix_to_regex(matrix) == '(G[AQ]D)'

def test_index_to_matrix(tr =tr ):
	r = _index_to_matrix(ind = [1,2,3], clone_df=tr.clone_df, pwmat = tr.pw_beta, col = 'cdr3_b_aa', centroid = None )
	assert isinstance(r, pd.DataFrame)
	assert r.shape == (24,14)
	assert r.columns.to_list() == ['C', 'A', 'S', 'S', 'G', 'K', 'Q', 'N', 'T', 'S', 'Q', 'F', 'T', 'F']

def test_index_to_seqs(ind = [1,2,3], clone_df = tr.clone_df, col = 'cdr3_b_aa'):
	dfnode   = clone_df.iloc[ind,].copy()
	seqs = dfnode[col].to_list()
	assert seqs == ['CASSGGANTGQLYF', 'CASSGKQNTSQFTF', 'CASSLKQGRSDYTF']

def test_index_to_regex_str():
	"""
	ind : list 
		iloc row index in clone_df
	clone_df : pd.DataFrae
		DataFrame with cdr3 sequence
	pwmat: np.ndarray
		Matrix with pairwise inofrmation
	col : str
		'cdr3_b_aa', 
	ntrim : int
		default 3, 
	ctrim : int
		default 2,  
	max_ambiguity : int
		default 3, the maximum number of amino acids at a position before a '.' wildcard is applied
	"""
	regex_pattern = _index_to_regex_str(ind = [1,2,3], 
										clone_df = tr.clone_df, 
										pwmat = tr.pw_beta, 
										centroid = None,
										col = 'cdr3_b_aa', 
										ntrim = 3, 
										ctrim = 2,  
										max_ambiguity = 3)

	assert regex_pattern == '(S[GL][GK][AQ][NG][RT][GS][DQ][LFY])'

