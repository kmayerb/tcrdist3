"""
test nieghbors
"""


def test_example():
	"""
	The purpose of this example is to show the use of 
	chosing thresholds based on background discovery rate
	"""
	import os
	import pandas as pd
	import numpy as np 
	from tcrdist.repertoire import TCRrep
	from tcrdist.neighbors import compute_ecdf, bkgd_cntl_nn2
	from tcrdist.automate import auto_pgen
	from tcrdist.regex import _index_to_regex_str, _index_to_seqs 
	from tcrdist.summarize import _summ, _dist_summ, _select, filter_gt, filter_is, test_for_subsets, test_for_almost_subsets

	fn = os.path.join('tcrdist', 'data', 'covid19', "m60_bkgd_test_input.csv")
	df_background = pd.read_csv(fn)
	
	tr_background = TCRrep(cell_df = df_background, 
					organism = "human", 
					chains= ['beta'], 
					compute_distances = False)

	fn = os.path.join('tcrdist', 'data', 'covid19', "m60_test_input.csv")
	df = pd.read_csv(fn)
	
	tr = TCRrep(cell_df = df, 
	            organism = 'human', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv')
	
	auto_pgen(tr)

	tr.compute_rect_distances(df = tr.clone_df, 
							  df2 = tr_background.clone_df, 
							  store = False)

	assert tr.rw_beta.shape[0] == tr.clone_df.shape[0]

	centers_df = bkgd_cntl_nn2(	tr = tr, 
								tr_background = tr_background,
								ctrl_bkgd = 2*10**-5, 
								col = 'cdr3_b_aa',
								ncpus = 2,
								thresholds = [x for x in range(0,50,2)])
	
	centers_df.sort_values(['target_hits'], ascending = False)
	
	# Let's add some useful summary info about publicity and number of unique seqs
	def tabulate_subjects(neighbor_df, clone_df, col_nn ='target_neighbors', col_seq= 'target_seqs'):
		# Tabulate the number of unique subjects at each node
		neighbor_df['nsubject'] = neighbor_df[col_nn].apply( lambda x: len(set(_select(clone_df, iloc_rows =x, col = 'subject'))))
		# Tabulate the number of unique sequences at each node
		neighbor_df['useq']    = neighbor_df[col_seq].apply( lambda x: len(set(x)))
		return neighbor_df

	tabulate_subjects(centers_df, tr.clone_df)
	
	# Let's add some useful info about Probability of generation
	centers_df['pgen_cdr3_b_aa'] = tr.clone_df.pgen_cdr3_b_aa.copy()
	centers_df['pgen_dist'] = _summ(df = tr.clone_df, 
		 indices = centers_df['target_neighbors'], 
		 column = 'pgen_cdr3_b_aa', 
		 f=_dist_summ)

	from scipy.stats import chi2_contingency
	n1 = 436
	n2 = 100000
	beta_re = 1 
	beta_dist = 1 
	centers_df['chi2re'] = [chi2_contingency(np.array( [[r['target_re_hits'], n1-r['target_re_hits']],[r['bkgd_re_hits'], n2-r['bkgd_re_hits']]]))[0] for _,r in centers_df.iterrows() ]
	centers_df['chi2dist'] = [chi2_contingency(np.array( [[r['target_hits'], n1-r['target_hits']],[r['background_hits'], n2-r['background_hits']]]))[0] for _,r in centers_df.iterrows() ]
	centers_df['chi2joint'] = [beta_re  * r['chi2re'] + beta_dist* r['chi2dist'] for _,r in centers_df.iterrows() ]

	# Rank and select non-redundant
	sorted_centers_df = centers_df.sort_values(['chi2joint'], ascending = False).copy()
	sorted_centers_df['novel'] = test_for_almost_subsets(sorted_centers_df['target_neighbors'], 5)
	sorted_filtered_centers_df = filter_is(sorted_centers_df, 'novel', 1).copy()
	#sorted_filtered_centers_df.to_csv("m60_centers_df.csv")
	





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


