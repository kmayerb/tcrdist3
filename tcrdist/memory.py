import numpy as np
import pandas as pd
import copy
import os
import parmap 
import secrets
from progress.bar import IncrementalBar
from hierdiff.tally import neighborhood_tally
from scipy import sparse
"""
Large datasets can consume more memory than is typically available. 

The memory module solves this problem with functions
that allow for TCRdistance computation to 
be broken down into fragments that can be stored on disk
and reassembled in a sparser format. 
"""

def _partition(l, n):
    """
    _partition breaks up a list in parts of size n

    Parameters
    ----------
    l : list
    
    n : int

    >>> _partition([1,2,3,4,5,6], 2)
    [[1, 2], [3, 4], [5, 6]]
    >>> _partition([1,2,3,4,5,6], 3)
    [[1, 2, 3], [4, 5, 6]]
        """
    n = int(n)
    r = [l[i:i + n] for i in range(0, len(l), n)]
    return r


def gen_sparse_rw_on_fragment(tcrrep, ind, outfile, matrix_name = 'rw_beta', max_distance=50):
	"""
	gen_sparse_rw_on_fragment generates a sparse matrix of distaances
	on a fragment of overall comparisons. 

	Suppose a clone_df has m unique clones. A full matrix is m x m. 
	this computes distances of the rectangular fragment 
	m[ind,] x m. The fragment is converted to a sparse compressed 
	format. All true zeros are set to 1, and all values above <max_distance>
	are set to one, for a massive space savings.

	tcrrep : TCRrep
		TCRrep instance 
	ind : list
		line of row index position specifying the rows of the clone_df
	outfile: str
		string name of the fragment
	matrix_name : str
		the matrix attribute to be stored ('rw_beta')

	Example
	-------
	See usage in tcrdist.rep_funcs.compute_pw_sparse_out_of_memory

	"""
	tr = copy.deepcopy(tcrrep)
	tr.compute_rect_distances(df = tr.clone_df.iloc[ind,], df2 = tr.clone_df, store = False)
	M = getattr(tr, matrix_name)
	# Set all true 0s to 1s
	M[M == 0] = 1
	# Set all values greater than max_distance to zero to create massive scarcity
	M[M > max_distance] = 0 
	# convert to scr matrix spa
	M = sparse.csr_matrix(M)
	# Write to outfile
	sparse.save_npz(outfile, M)
	del tr
	del M
	return True


def collapse_csrs(list_of_csr_filenames, axis = 0):
	"""
	Given a list of filenames refering to the sparse
	repressentation of ordered row-wise (axis = 0) chunks
	of a matrix, collapse chunks into a single 
	matrix in sparse format.

	Parameters
	----------
	list_of_csr_filenames : list of strings
		list of string pointing to filenames of .npz sparse matrices
	axis : int
		if 0, combine with vstack, else if 1

	Example
	-------
	>>> import numpy as np
	>>> from scipy import sparse
	>>> A = np.array([[1, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 1], [0, 0, 0, 2, 0, 0]])
	>>> B = np.array([[3, 3, 3, 1, 0, 0], [5, 5, 2, 0, 0, 1], [6, 7, 0, 2, 0, 0]])
	>>> A_s = sparse.csr_matrix(A)
	>>> B_s = sparse.csr_matrix(B)
	>>> sparse.save_npz("A.npz", A_s)
	>>> sparse.save_npz("B.npz", B_s)
	>>> AB_s = collapse_csrs(["A.npz","B.npz"], axis = 0)
	>>> AB = AB_s.todense()
	>>> assert np.all(AB == np.concatenate([A,B], axis = 0))
	"""
	assert axis in [0,1], "axis must be 0 (for hstack on rows) or 1 (for vstack on columns)."
	if axis == 0:
		full_spase_matrix = sparse.vstack([sparse.load_npz(f_npz) for f_npz in list_of_csr_filenames])
	if axis == 1:
		full_spase_matrix = sparse.hstack([sparse.load_npz(f_npz) for f_npz in list_of_csr_filenames])
	return full_spase_matrix





# def gen_rw_distance_on_fragment(tcrrep, ind, outfile, rw = 'rw_beta'):
# 	"""
# 	gen_rw_distance_on_fragment
# 	"""
# 	import copy
# 	tr = copy.deepcopy(tcrrep)
# 	tr.compute_rect_distances(df = tr.clone_df.iloc[ind,], df2 = tr.clone_df)
# 	with open(outfile, 'wb') as frag:
# 		np.save(frag, getattr(tr, rw), allow_pickle = False)
# 	del tr
# 	return outfile


# def gen_ndif_on_fragment(tcrrep, infile, ind, outfile):
# 	import copy
# 	tr = copy.deepcopy(tcrrep)
# 	with open(infile ,"rb") as frag:
# 		rwmat = np.load(frag)
# 	ndif = neighborhood_tally(df_pop = tr.clone_df, pwmat = rwmat, x_cols = ['epitope'], df_centroids=tr.clone_df.iloc[ind,], count_col='count', knn_neighbors=None, knn_radius=50)
# 	ndif.to_csv(outfile, index = False)
# 	del ndif
# 	return outfile

# def gen_ndif_on_fragment2(tcrrep, infile, ind, outfile):
# 	import copy
# 	tr = copy.deepcopy(tcrrep)
# 	with open(infile ,"rb") as frag:
# 		rwmat = np.load(frag)
# 	ndif = neighborhood_tally(df_pop = tr.clone_df, pwmat = rwmat, x_cols = ['vzv'], df_centroids=tr.clone_df.iloc[ind,], count_col='single', knn_neighbors=None, knn_radius=50)
# 	ndif.to_csv(outfile, index = False)
# 	del ndif
# 	return outfile

# def concat(dest,fragments):
# 	from progress.bar import IncrementalBar
# 	bar = IncrementalBar(f'Processing Files', max = len(fragments), suffix='%(percent)d%%')
# 	mysep = ","
# 	mycatfile = f"{dest}/all.ndif25.csv"
# 	counter = 0
# 	for f in fragments:
# 		df = pd.read_csv(os.path.join(f), sep = ',')
# 		if counter == 0:
# 			df.to_csv(mycatfile, header=True, index = False, sep = mysep)
# 		if counter > 0:
# 			df.to_csv(mycatfile, header=False, index = False, sep = mysep, mode = "a")
# 		del df
# 		counter =+1 
# 		bar.next()
# 	bar.next(); bar.finish()