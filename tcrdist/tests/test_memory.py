import pytest
import numpy as np
import pandas as pd
import os
import numpy as np
from scipy import sparse

from tcrdist.memory import _partition
from tcrdist.memory import collapse_csrs
from tcrdist.memory import gen_sparse_rw_on_fragment
from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory

def test_partion_size_2_on_6():
	"""
	UNIT TEST _partition
	"""
	assert _partition([1,2,3,4,5,6], 2) == [[1, 2], [3, 4], [5, 6]]

def test_partion_size_3_on_6():
	assert _partition([1,2,3,4,5,6], 3) == 	[[1, 2, 3], [4, 5, 6]]

def test_partion_size_3_on_5():
	assert _partition([1,2,3,4,5], 3) == 	[[1, 2, 3], [4, 5]]


def test_collapse_csrs_axis_0():
	"""
	UNIT TEST _collapse_csrs (axis = 0)
	"""
	A = np.array([[1, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 1], [0, 0, 0, 2, 0, 0]])
	B = np.array([[3, 3, 3, 1, 0, 0], [5, 5, 2, 0, 0, 1], [6, 7, 0, 2, 0, 0]])
	A_s = sparse.csr_matrix(A)
	B_s = sparse.csr_matrix(B)
	sparse.save_npz("A.npz", A_s)
	sparse.save_npz("B.npz", B_s)
	AB_s = collapse_csrs(["A.npz","B.npz"], axis = 0)
	AB = AB_s.todense()
	assert np.all(AB == np.concatenate([A,B], axis = 0))


def test_collapse_csrs_axis_1():
	"""
	UNIT TEST _collapse_csrs (axis = 1)
	"""
	A = np.array([[1, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 1], [0, 0, 0, 2, 0, 0]])
	B = np.array([[3, 3, 3, 1, 0, 0], [5, 5, 2, 0, 0, 1], [6, 7, 0, 2, 0, 0]])
	A_s = sparse.csr_matrix(A)
	B_s = sparse.csr_matrix(B)
	sparse.save_npz("A.npz", A_s)
	sparse.save_npz("B.npz", B_s)
	AB_s = collapse_csrs(["A.npz","B.npz"], axis = 1)
	AB = AB_s.todense()
	assert np.all(AB == np.concatenate([A,B], axis = 1))

def test_gen_sparse_rw_on_fragment():
	"""
	UNIT TEST ON _gen_sparse_rw_on_fragment, 
		TESTS THAT A SPARSE REPRESENTATION MATRIX IS COMPUTED, SAVED TO .npz 
		AND RETURNED AS CORRECT DIMENSIONS
	"""
	df = pd.read_csv("dash.csv")
	tr = TCRrep(cell_df = df,               #(2)
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
	            compute_distances = False)
	r = gen_sparse_rw_on_fragment(tcrrep = tr, ind = list(range(0,10)), outfile = "TESTONLY.csr", matrix_name = "rw_beta", max_distance = 1000)
	assert r == True
	assert os.path.isfile("TESTONLY.csr.npz")
	sparse_matrix_fragment = sparse.load_npz("TESTONLY.csr.npz")
	assert sparse_matrix_fragment.todense().shape[0] == 10
	assert sparse_matrix_fragment.todense().shape[1] == tr.clone_df.shape[0]









