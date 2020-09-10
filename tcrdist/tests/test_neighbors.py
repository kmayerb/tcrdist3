"""
test nieghbors
"""


def test_example():
	"""
	The purpose of this example is to show the use of 
	chosing thresholds based on background discovery rate
	"""
	import pandas as pd
	import numpy as np 
	from tcrdist.repertoire import TCRrep
	from tcrdist.neighbors import compute_ecdf, bkgd_cntl_nn
	from tcrdist.automate import auto_pgen
	from tcrdist.regex import _index_to_regex_str, _index_to_seqs 
	# Set the acceptable fraction of background that may 
	# be in the nieghborhood of a TCR


	df_background = pd.read_csv("m60_bkgd_test_input.csv")
	
	tr_background = TCRrep(cell_df = df_background, 
					organism = "human", 
					chains= ['beta'], 
					compute_distances = False)

	df = pd.read_csv("m60_test_input.csv")
	
	tr = TCRrep(cell_df = df, 
	            organism = 'human', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv')
	
	auto_pgen(tr)

	tr.compute_rect_distances(df = tr.clone_df, 
							  df2 = tr_background.clone_df, 
							  store = False)

	assert tr.rw_beta.shape[0] == tr.clone_df.shape[0]

	centers_df = bkgd_cntl_nn(	tr = tr, 
								tr_background = tr_background,
								ctrl_bkgd = 2*10**-5, 
								col = 'cdr3_b_aa',
								ncpus = 6,
								thresholds = [x for x in range(0,50,2)])
	
	centers_df.sort_values(['target_hits'], ascending = False)


def test_compute_ecdf():
	"""
	Test import, and ecdf with some or all optional arguments
	"""
	from tcrdist.neighbors import compute_ecdf
	import numpy as np
	assert isinstance(compute_ecdf(data = np.array([1,2,3,4,5,6])), np.ndarray)
	assert isinstance(compute_ecdf(data = np.array([1,2,3,4,5,6]), thresholds = np.array([3])), np.ndarray)
	r = compute_ecdf(data = np.array([1,2,3,4,5,6]), counts =np.array([100,100,100,100,100,100]), thresholds = np.array([3]))
	assert np.all( np.isclose(r, np.array([0.50083195]) ) )


