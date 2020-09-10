"""
tcrdist/regex.py 

Regex tools for defining regex patterns from a list of seququences, 
aligned using computepal_motif.
"""

import numpy as np
import pandas as pd
from palmotif import compute_pal_motif
import re

def _list_to_regex_component(l, max_ambiguity = 3):
	""" 
	list of str to regex 

	Parameters
	----------

	l : list
		list of strings
	max_ambiguity : int
		default is 3, more than max results in '.' rather than a set []
	
	Example
	-------
	>>> _list_to_regex_component(['A'])
	'A'
	>>> _list_to_regex_component(['A','T','C'])
	'[ATC]'
	>>> _list_to_regex_component(['A','T','C','G'])
	'.'
	>>> _list_to_regex_component(['A','-','C'])
	'[AC]?'
	"""
	if len(l) < 1:
		s = '.?'
	elif len(l) < 2:
		s = l[0]
	elif len(l) <= max_ambiguity:
		s = f"[{''.join(l)}]"
	else:
		s = '.'
	
	if s == '-':
		s = s.replace('-', '.?') 
	elif s.find('-') != -1:
		s = s.replace('-', '') 
		s = f"{s}?"
	return s


def _index_to_matrix(ind, clone_df, pwmat = None, col = 'cdr3_b_aa', centroid = None ):
	"""
	Example
	-------


	"""
	dfnode   = clone_df.iloc[ind,].copy()

	seqs = dfnode[col].to_list()
	if centroid is None:
		pwnode   = pwmat[ind,:][:,ind].copy()
		iloc_idx = pwnode.sum(axis = 0).argmin() 
		centroid = dfnode[col].to_list()[iloc_idx]
	
	matrix, stats = compute_pal_motif(seqs = seqs, centroid = centroid)
	return matrix


def _matrix_to_regex(matrix, ntrim = 3, ctrim = 2, max_ambiguity = 3):
	"""
	Example
	-------
	>>> import numpy as np
	>>> import pandas as pd
	>>> arr = np.array([[ 0.        ,  3.32910223,  0.        ,  0.        ,  0.41764668,
	         0.        ,  0.        ,  0.        ],
	       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
	         0.        ,  3.60160675,  4.5849625 ],
	       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
	         0.        ,  0.        ,  0.        ],
	       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
	         4.5849625 ,  0.        ,  0.        ],
	       [ 4.5849625 ,  0.        ,  0.41764668,  0.        ,  0.        ,
	         0.        ,  0.        ,  0.        ],
	       [ 0.        ,  0.        ,  0.        ,  0.        ,  2.66666667,
	         0.        ,  0.        ,  0.        ],
	       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
	         0.        ,  0.        ,  0.        ],
	       [ 0.        ,  0.        ,  0.        , -0.01922274,  0.        ,
	         0.        ,  0.        ,  0.        ]])
	>>> matrix = pd.DataFrame(arr, index = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G'], columns = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G'])
	>>> matrix_to_regex(matrix)
	'(G[AQ]D)'
	"""
	regex_parts =[matrix.iloc[:,][matrix.iloc[:,i] != 0].index.to_list() for i in range(ntrim ,matrix.shape[1] - ctrim)]
	regex_parts_str = [_list_to_regex_component(x, max_ambiguity = max_ambiguity) for x in regex_parts] 
	regex_pattern = f"({''.join(regex_parts_str)})"
	return regex_pattern


def _index_to_regex_str(ind, 
						clone_df, 
						pwmat, 
						centroid = None,
						col = 'cdr3_b_aa', 
						ntrim = 3, 
						ctrim = 2,  
						max_ambiguity = 3):
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
	
	mat = _index_to_matrix(ind = ind, clone_df= clone_df, pwmat = pwmat, col = col, centroid = centroid )
	
	regex_str = _matrix_to_regex(matrix = mat, 
		ntrim = ntrim, 
		ctrim = ctrim, 
		max_ambiguity =  max_ambiguity)

	return(regex_str)


def _multi_regex(regex , bkgd_cdr3):
	"""
	Search a regex pattern in a list of string 
	"""
	result = [re.search(string = s, pattern = regex ) for s in bkgd_cdr3] 
	result = [1 if (x is not None) else None for x in result]
	return result

def _index_to_seqs(ind, clone_df, col):
	dfnode   = clone_df.iloc[ind,].copy()
	seqs = dfnode[col].to_list()
	return seqs 
