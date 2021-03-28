"""
centers 

Module contains functions for evaluating TCRs as center(oids) of meta-clonotypes.

find_center 
"""
import warnings 
import numpy as np
from tcrdist.ecdf import distance_ecdf

def calc_radii(tr, tr_bkgd, chain = 'beta', ctrl_bkgd = 10**-5, use_sparse = True, max_radius=50, chunk_size=100, **kwargs):
	"""
	Simply find maximum radii based on an antigen enriched repertoires <tr> 
	and a background antigen-naive background. 

	IMPORTANT: This function will work with a precomputed TCRrep.rw_chain matrix if it 
	is supplied and matches row dimensions of tr_bkgd. This way the user does 
	not have to recompute it each time they with to change ctrl_bkgd. It 
	also allows them to customize the distance prior to this stage, although 
	**kwargs can be passed to the compute_sparse_rect_distances or 
	compute_rect_distances function

	Parameters
	----------
	tr : tcrdist.repertoire.TCRrep
		An antigen enriched repertoires TCRrep object
	tr_bkgd : tcrdist.repertoire.TCRrep
		A background antigen-naive background TCRrep object
	chain : str
		e.g, 'beta', ctrl_bkgd = 10**-5, use_sparse = True, max_radius=50, chunk_size=50
	ctrl_bkgd : float
		e.g., 10**-5, 
	use_sparse : bool 
		If True, uses a sparse implementation, 
	max_radius : int 
		Values beyond this max_radius are set to 0 in sparcification
	chunk_size : int
		How many rows to process at a time, based on memory available (100 is a good default, but for more see notes)
	**kwargs will be passed to either compute_sparse_rect_distances or compute_rect_distances
	Returns
	-------
	max_radii : list
		List of length equal to the numer of rows in <tr> clond_df. 
		These are radius at which the number of expected background 
		sequences are 'controlled' at a rate of <ctlr_bkgd> 

	Notes
	-----
	TODO: Discuss chunk size and memory
	"""
	assert chain in ['alpha','beta','gamma', 'delta']
	if 'weights' not in tr_bkgd.clone_df.columns:
		warnings.warn("No weights provided in background repertoire, setting to 1")
		tr_bkgd.clone_df['weights'] = 1

	# USER MAY HAVE ALRRADY COMPUTED TCRrep.rw_
	if getattr(tr, f"rw_{chain}", None) is not None:
		if getattr(tr, f"rw_{chain}").shape[1] == tr_bkgd.clone_df.shape[0]:
			print(f"IT APPEARS THAT (TCRrep.rw_{chain}) HAS ALREDY BEEN COMPUTED AND MATCHES BACKGROUND TCRrep SIZE")
			print(f"USING EXISTING (TCRrep.rw_{chain}). SET TCRrep.rw_{chain} = None IF YOU WANT TO RECOMPUTE IT.")
	else:
		if use_sparse:
			print(f"COMPUTING SPARSE RECT MATRIX TO FIND RADIUS: (TCRrep.rw_{chain})")
			print(f"USING {tr.cpus} CPUS")
			tr.compute_sparse_rect_distances(df = tr.clone_df, 
											 df2 = tr_bkgd.clone_df,
											 radius=max_radius,
											 chunk_size=chunk_size,
											 **kwargs)
		else:
			print(f"COMPUTING FULL RECT MATRIX TO FIND RADIUS, (TCRrep.rw_{chain})")
			tr.compute_rect_distances(df = tr.clone_df, 
									  df2 = tr_bkgd.clone_df, 
									  store = False,
									  **kwargs)
	print(f"COMPUTING ECDFS PER TCR, TO FIND APPROPRIATE MAX RADII AT {ctrl_bkgd}")
	thresholds, ecdfs = distance_ecdf(pwrect = getattr(tr, f"rw_{chain}"), 
						  thresholds = np.array(range(0,max_radius, 2)), 
						  weights= tr_bkgd.clone_df.weights, 
						  pseudo_count=0, 
						  skip_diag = False, 
						  absolute_weight = True)
	# Based on acceptable ctrl_bkgd, we find max acceptable radi from each TCR
	#import pdb; import pdb; pdb.set_trace()
	all_radii = [pd.Series(x, index = thresholds) for x in ecdfs]
	max_radii = [s[s<=ctrl_bkgd].last_valid_index() for s in all_radii]
		# WARNING: There is a potential BUG in the output of the above line. 
		# That is, iF a radius is None (the next line will fail, thus set Nones to 0.
	max_radii = [i if (i is not None) else 0 for i in max_radii]
	print(f"RETURNING LIST OF MAX RADII")
	
	return max_radii, thresholds, ecdfs 







def find_centers_beta(	
	background_filename,
	target_filename,
	ncpus,
	min_nsubject,
	ctrl_bkgd = 10**-5, 
	prefilter = False):
	import os
	import pandas as pd
	import numpy as np 
	from tcrdist.repertoire import TCRrep
	from tcrdist.neighbors import compute_ecdf, bkgd_cntl_nn2
	from tcrdist.automate import auto_pgen
	from tcrdist.rep_diff import neighborhood_diff
	from tcrdist.summarize import _summ, _dist_summ, _select, filter_gt, filter_is, test_for_subsets, test_for_almost_subsets
	import scipy.sparse
	
	df_background       =   pd.read_csv(background_filename)
	print(df_background)
	tr_background       =   TCRrep( cell_df = df_background.copy(), 
									organism = "human", 
									chains= ['beta'], 
									compute_distances = False)

	df_mira             =   pd.read_csv(target_filename)
	df_mira             =   df_mira[['subject','cell_type','v_b_gene', 'j_b_gene', 'cdr3_b_aa']]
	print(df_mira)
	tr                  =   TCRrep(	cell_df = df_mira.copy(), 
									organism = 'human', 
									chains = ['beta'], 
									db_file = 'alphabeta_gammadelta_db.tsv',
									store_all_cdr = False,
									compute_distances = True)
	
	if prefilter:
		# We can greatly cut down on the number of searches if we drop centroids without minimum publicicity
		nn_df = neighborhood_diff(clone_df= tr.clone_df,
								  pwmat = tr.pw_beta,
								  count_col = 'count',
								  x_cols = ['cell_type'],
								  knn_radius = 37)
		def tabulate_publicity(neighbor_df, clone_df, col_nn ='neighbors'):
			# Tabulate the number of unique subjects at each node
			neighbor_df['nsubject'] = neighbor_df[col_nn].apply( lambda x: len(set(_select(clone_df, iloc_rows =x, col = 'subject'))))
			return neighbor_df
		print(f"TABULATING PUBLIC CLUSTERS")
		nn_df = tabulate_publicity(nn_df, tr.clone_df)
		nn_df = filter_gt(nn_df, 'nsubject' , min_nsubject)
	
		if nn_df.shape[0] == 0:
			centers_df = pd.DataFrame({}, columns =  ['cdr3_b_aa','v_b_gene','j_b_gene','pgen','max_radi','target_hits','bkgd_hits','bkgd_hits_weighted','bkgd_total','ctrl','ctrl_weighted','target_misses','TR','TR2','BR_weighted','RR_weighted','OR_weighted','chi2dist','target_neighbors','target_seqs','background_neighbors','background_seqs','background_v','background_j','regex','target_re_hits','bkgd_re_hits','bkgd_re_weighted_hits','TR_re','BR_re_weighted','RR_re_weighted','OR_re_weighted','chi2re','chi2joint','nsubject'])
			tr.pw_beta[tr.pw_beta == 0]  = 1 # set true zeros to 1
			tr.pw_beta[tr.pw_beta > 50]  = 0 # ignores everything less than 100
			pw_beta_sparse = scipy.sparse.csr_matrix(tr.pw_beta)	
			return centers_df, pw_beta_sparse
	
		tr.clone_df = tr.clone_df.loc[nn_df.index, :].reset_index(drop = True)
		del nn_df
		# Compute pairwise again with filtered set
		tr.compute_distances()
		# compute pgens automatically, currently parmap will max out cpus on this step 
		
		

	print("COMPUTING PROBABILITY OF GENERATION")
	auto_pgen(tr)
	print(f"COMPUTING RECT DIST {tr.clone_df.shape[0]}x{tr_background.clone_df.shape[0]}")
	tr.compute_rect_distances(df = tr.clone_df, 
							  df2 = tr_background.clone_df, 
							  store = False)

	assert tr.rw_beta.shape[0] == tr.clone_df.shape[0]

	centers_df = bkgd_cntl_nn2(	tr = tr, 
								tr_background = tr_background,
								ctrl_bkgd = ctrl_bkgd,  #ctrl_bkgd = 2*10**-5
								weights =tr_background.clone_df.weights,
								col = 'cdr3_b_aa',
								ncpus = ncpus,
								thresholds = [x for x in range(0,38,2)], # Settign 38 as the max radius
								generate_regex = True, 
								test_regex = True)

	def tabulate_publicity(neighbor_df, clone_df, col_nn ='neighbors'):
		# Tabulate the number of unique subjects at each node
		neighbor_df['nsubject'] = neighbor_df[col_nn].apply( lambda x: len(set(_select(clone_df, iloc_rows =x, col = 'subject'))))
		return neighbor_df

	centers_df = tabulate_publicity(neighbor_df = centers_df, clone_df = tr.clone_df, col_nn ='target_neighbors')
		
	tr.rw_beta[tr.rw_beta == 0]  = 1 # set true zeros to 1
	tr.rw_beta[tr.rw_beta > 50]  = 0 # ignores everything less than 100
	rw_beta_sparse = scipy.sparse.csr_matrix(tr.rw_beta)
	#scipy.sparse.save_npzz(output_matrix_filename, rw_beta_sparse)
	return centers_df, rw_beta_sparse


def rank_centers(centers_filename = None, centers_df = None, rank_column = 'chi2joint', min_nsubject = 2, min_nr = 1):
	"""
	This function takes the output of tcrdist.neighbors.bkgd_cntl_nn2(), 
	a set of scored metaclonotypes (centers - TCRs + radius) 
	and ranks them by chi2 statistics, 
	prioritizing those that include lots of target sequences 
	while minimizing inclusion of background sequences. 
	
	Parameters
	----------
	centers_filename : str or None
		User can only provide centers_df or centers_filename but not both
		The filepath to a file containing metaclonotype centers information, generally produced with 
		tcrdist.neighbors.bkgd_cntl_nn2()
	centers_df : DataFrame or None
		User can only provide centers_df or centers_filename but not both.
		The Pandas DataFraem containing metaclonotype centers information, generally produced with 
		tcrdist.neighbors.bkgd_cntl_nn2()
	rank_column : str
		Default : 'chi2joint' (or 'chi2joint' (radius+motif averaged) or chi2re'(using motif only), 'chi2dist' (using radius only) 
	min_nsubject : int
		Default 2, (minimum publicity of the meta-clonotype). 
		That is, the minimum number of unique subjects contributing TCRs 
		among a group of biochemically TCRs to form a meta-clonotype. 
	min_nr : int
		Default 1, (minimum non-redundancy). Once the metaclonotypes are ranked, 
		the function requires that lower ranked meta-clonotypes to have a minimum number
		<min_nr> of new sequences not already spanned by a higher ranked meta-clonotype. 
	
	Returns
	-------
	df : DataFrame
	
	"""
	import pandas as pd
	import ast
	from tcrdist.summarize import filter_gt, filter_is, test_for_subsets, test_for_almost_subsets
	
	if centers_filename is not None and centers_df is not None:
		raise ValueError("rank centers can use <centers_filename> or <centers_df> but not both")
	if centers_df is None:
		df = pd.read_csv(centers_filename)
	else: 
		df = centers_df.copy()

	# VERY IMPORTANT NOTE, pandas reads lists as strings '[1,2]'; so we use ast.literal_eval to convert back t a list 
	if not isinstance(df['target_neighbors'][0], list):
		df['target_neighbors'] = df['target_neighbors'].apply(lambda s: list(ast.literal_eval(s)))
	df = df.sort_values(rank_column, ascending = False)
	df['novel'] = test_for_almost_subsets(df['target_neighbors'], min_nr)
	df = filter_gt(df, 'nsubject', min_nsubject).copy()
	df = filter_is(df, 'novel', 1).copy()
	return df


import re
import os
import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep
from tcrdist.adpt_funcs import _valid_cdr3
import scipy.sparse
import ast

def check_tsv_csv(filename, check_column_names = ['cdr3_b_aa', 'v_b_gene', 'j_b_gene']):
	"""
	To avoid problems with .tsv or .csv, check for appropriate seperator based on expected columns
	"""
	#shape_csv = pd.read_csv(filename, sep = ",").shape
	columns_csv = pd.read_csv(filename, sep = ",").columns
	#shape_tsv = pd.read_csv(filename, sep = "\t").shape
	columns_tsv = pd.read_csv(filename, sep = "\t").columns
	
	if set(check_column_names) - set(columns_csv) == set():
		sep = ","
	elif set(check_column_names) - set(columns_tsv) == set():
		sep = "\t"
	else: 
		sep = False
		raise IOError("File provided does not appear to be either a .csv or .tsv or lacks required columns {check_column_names}")
	return sep


def centers_v_bulk(search_filename, bulk_filename, sep_search_filename = "\t"):
	"""
	"""
	# Get appropriate seperator (in this case we expect "\t", trust than verify)
	sep_bulk_filename = check_tsv_csv(filename = bulk_filename, check_column_names = ['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'templates', 'productive_frequency','count'])
	bulk_df = pd.read_csv(bulk_filename, sep = sep_bulk_filename)

	# Adaptive uses the terms templates, which are synonymous with counts in our nomenclature
	bulk_df['count'] = bulk_df['templates'].copy()
	# Ensure that all bulk sequences have valid cdr3s
	v = bulk_df.cdr3_b_aa.apply(lambda x: _valid_cdr3(x))
	# Ensure length CDR3 > 5
	ls = bulk_df.cdr3_b_aa.apply(lambda x: len(x) > 5)
	bulk_df = bulk_df[(v) & (ls)]

	# Select only important columns
	bulk_df = bulk_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'templates', 'productive_frequency','count']]
	# Assign a cid for tracking purposes
	bulk_df['cid'] = [f"cid{i}" for i in bulk_df.index]
	# Load clean bulk data to assign CDR1,2,3
	tr_bulk = TCRrep(	cell_df = bulk_df, 
						organism = 'human', 
						chains = ['beta'], 
						compute_distances= False)

	# Here we expect comma or tab, but we check file for correct seperator (TAB, TWO MANY COMMAS CAUSES A FAILURE)
	#sep_search_filename = check_tsv_csv(
	#	filename = search_filename, 
	#	check_column_names = ['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'max_radi' ,'pgen','regex'])
	
		# Read search file
	search_df = pd.read_csv(search_filename, sep = sep_search_filename )
	# Assign a count column of 1, for purposes of loading data into TCRrep
	search_df['count'] = 1

	# If source and index columns are missing we provide them, select only the relevant columns
	try:
		search_df = search_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'count', 'max_radi' ,'pgen','index','source','regex','target_hits', 'target_seqs']].copy()
	except KeyError:
		search_df['source'] = search_filename
		search_df['index'] = search_df.index.to_list()
		search_df = search_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'count', 'max_radi' ,'pgen','index','source','regex','target_hits', 'target_seqs']].copy()

	# Load the search file into TCRrep instance, getting CDR1, 2, 2.5 by vgene name infrerence
	tr_search = TCRrep( cell_df = search_df.copy(),
						organism = 'human', 
						chains = ['beta'], 
						compute_distances= False)

	# Compute Rectangular distance (search on rows, bulk on columns)
	tr_search.compute_rect_distances(df = tr_search.clone_df, 
									 df2 = tr_bulk.clone_df, 
									 store = False)

	# Each search sequence has a max radius 
	max_radi = tr_search.clone_df.max_radi.values
	# Convert to a 2D array to permit broadcasting (n,1) with (n,m) rw_matrix
	max_radi = max_radi.reshape(len(max_radi), 1)

	# Assert that [n,1] aligns with [n,m] 
	assert max_radi.shape[0] == tr_search.rw_beta.shape[0]
	# Get index of True, where 
	ij = tr_search.rw_beta < max_radi
	# Get the column index of each bulk sequences within the row specific radius
	icol = [np.where(x) for x in ij]
	bulk_hits = np.sum(ij, axis = 1)
	assert np.all(bulk_hits == [len(x[0]) for x in icol])
	# Retrieve the distance of sequences with the max radius for each row
	idist =[tr_search.rw_beta[i,j] for i,j in enumerate(icol)]

	# Retrieve sequences from the bulk clone_df
	iseqs = [tr_bulk.clone_df['cdr3_b_aa'].iloc[x].to_list() for x in icol]
	ivgenes = [tr_bulk.clone_df['v_b_gene'].iloc[x].to_list() for x in icol]
	ijgenes = [tr_bulk.clone_df['j_b_gene'].iloc[x].to_list() for x in icol]
	
	# Retrieve abundances from the bulk clone df
	itemplates = [tr_bulk.clone_df['templates'].iloc[x].to_list()  for x in icol]
	icounts = [tr_bulk.clone_df['count'].iloc[x].to_list() for x in icol]
	ifreqs = [tr_bulk.clone_df['productive_frequency'].iloc[x].to_list()  for x in icol]

	assert [np.sum(x) for x in itemplates] == [np.sum(x) for x in icounts]
	isumtemplates = [np.sum(x) for x in itemplates]
	isumcounts    = [np.sum(x) for x in icounts]
	isumfreqs     = [np.sum(x) for x in ifreqs ]

	df_summ = pd.DataFrame({ 'bulk_sum_freq'     : isumfreqs,
							 'bulk_sum_templates': isumtemplates,
							 'bulk_sum_counts'   : isumcounts,
							 'bulk_seqs'      	 : iseqs,
							 'bulk_v_genes'   	 : ivgenes,
							 'bulk_j_genes'   	 : ijgenes,
							 'bulk_distances' 	 : idist,
							 'bulk_templates' 	 : itemplates,
							 'bulk_counts'    	 : icounts, 
							 'bulk_freqs'     	 : ifreqs})

	assert tr_search.clone_df.shape[0] == search_df.shape[0]

	#search_df = pd.read_csv(search_filename, sep = sep)
	###### <<>>>>>
	result_df = pd.concat([tr_search.clone_df, df_summ], axis = 1)
	result_df['sourcefile'] = bulk_filename
	result_df['searchfile'] = search_filename
	
	regex_on_bulk = [[re.search(pattern = r['regex'], string =  s) for s in r['bulk_seqs']] for _,r in result_df.iterrows()]
	bulk_regex_match = [[True if (x is not None) else False for x in sublist] for sublist in regex_on_bulk]
	result_df['bulk_regex_match'] = bulk_regex_match
	result_df['bulk_sum_freqs_regex_adj'] = [pd.Series(r['bulk_freqs'], dtype = "float64")[pd.Series(r['bulk_regex_match'], dtype = 'bool')].sum() for i,r in result_df.iterrows()]
	result_df['bulk_sum_templates_regex_adj'] = [pd.Series(r['bulk_templates'], dtype = "int32")[pd.Series(r['bulk_regex_match'], dtype = 'bool')].sum() for i,r in result_df.iterrows()]
	result_df['bulk_sum_counts_regex_adj'] = [pd.Series(r['bulk_counts'], dtype = "int32")[pd.Series(r['bulk_regex_match'],  dtype = 'bool')].sum() for i,r in result_df.iterrows()]
	
	# Tabulating Tcrdist == 0
	ij = tr_search.rw_beta == 0
	icol = [np.where(x) for x in ij]
	itemplates = [tr_bulk.clone_df['templates'].iloc[x].to_list()  for x in icol]
	icounts = [tr_bulk.clone_df['count'].iloc[x].to_list() for x in icol]
	ifreqs = [tr_bulk.clone_df['productive_frequency'].iloc[x].to_list()  for x in icol]
	assert [np.sum(x) for x in itemplates] == [np.sum(x) for x in icounts]
	isumcounts0 = [np.sum(x) for x in icounts]
	isumtemplates0 = [np.sum(x) for x in itemplates]
	isumfreqs0     = [np.sum(x) for x in ifreqs ]
	result_df['bulk_sum_freqs_tcrdist0'] = isumfreqs0
	result_df['bulk_sum_templates_tcrdist0'] = isumtemplates0
	result_df['bulk_sum_counts_tcrdist0'] = isumcounts0


	# Tabulate Perfect Match 
	result_df['cdr3_exact_match'] = [[s == r['cdr3_b_aa'] for s in r['bulk_seqs']] for _,r in result_df.iterrows()]
	result_df['v_gene_exact_match'] = [[s == r['v_b_gene'] for s in r['bulk_v_genes']] for _,r in result_df.iterrows()]
	def safely_compare_boolan_lists(a,b):
		if len(a) < 1:
			r = list()
		else:
			r = np.array(a) & np.array(b)
			r = list(r)
		return r
	result_df['vcdr3_exact_match'] = [safely_compare_boolan_lists(r['cdr3_exact_match'], r['v_gene_exact_match']) for _,r in result_df.iterrows()]
	# Summarize Exact Matches
	result_df['bulk_sum_freqs_vcdr3match'] = [np.sum(pd.Series(r['bulk_freqs'],dtype = 'float64')[r['vcdr3_exact_match']]) for _,r in result_df.iterrows()]
	result_df['bulk_sum_templates_vcdr3match'] = [np.sum(pd.Series(r['bulk_templates'],dtype = 'int32')[r['vcdr3_exact_match']]) for _,r in result_df.iterrows()]
	result_df['bulk_sum_counts_vcdr3match'] = [np.sum(pd.Series(r['bulk_counts'],dtype = 'int32')[r['vcdr3_exact_match']]) for _,r in result_df.iterrows()]

	# Tabulate Perfect Matches to any of the Target Seqs
	result_df['target_seqs'] = [ast.literal_eval(x) for x in result_df['target_seqs']]
	result_df['bulk_seq_within_targetseqs'] = [[s in r['target_seqs'] for s in r['bulk_seqs']] for _,r in result_df.iterrows()]
	result_df['bulk_sum_freqs_within_targetset'] = [np.sum(pd.Series(r['bulk_freqs'],dtype = 'float64')[r['bulk_seq_within_targetseqs']]) for _,r in result_df.iterrows()]
	result_df['bulk_sum_templates_within_targetset'] = [np.sum(pd.Series(r['bulk_templates'],dtype = 'int32')[r['bulk_seq_within_targetseqs']]) for _,r in result_df.iterrows()]
	result_df['bulk_sum_counts_within_targetset'] = [np.sum(pd.Series(r['bulk_counts'],dtype = 'int32')[r['bulk_seq_within_targetseqs']]) for _,r in result_df.iterrows()]


	desired_output_column_order = [  'cdr3_b_aa',
									 'v_b_gene',
									 'j_b_gene',
									 'max_radi',
									 'pgen',
									 'index',
									 'source',
									 'target_seqs',
									 'regex',
									 'cdr1_b_aa',
									 'cdr2_b_aa',
									 'pmhc_b_aa',
									 'count',
									 'clone_id',
									 'sourcefile',
									 'searchfile',
									 
									 'bulk_distances',
									 'bulk_templates',
									 'bulk_counts',
									 'bulk_freqs',
									 'bulk_seqs',
									 'bulk_v_genes',
									 'bulk_j_genes',
									 'bulk_regex_match',

									 'bulk_sum_freq',
									 'bulk_sum_templates',
									 'bulk_sum_counts',
									  
									 'bulk_sum_freqs_regex_adj',
									 'bulk_sum_templates_regex_adj',
									 'bulk_sum_counts_regex_adj',
									 
									 'bulk_sum_freqs_tcrdist0',
									 'bulk_sum_templates_tcrdist0',
									 'bulk_sum_counts_tcrdist0',
									 #'cdr3_exact_match',
									 #'v_gene_exact_match',
									 #'vcdr3_exact_match',
									 'bulk_sum_freqs_vcdr3match',
									 'bulk_sum_templates_vcdr3match',
									 'bulk_sum_counts_vcdr3match',

									 'bulk_sum_freqs_within_targetset',
									 'bulk_sum_templates_within_targetset',
									 'bulk_sum_counts_within_targetset']

	result_df = result_df[desired_output_column_order] 

	# COMPRESS SPARSE MATRIX FOR LATER REFRENCE
	tr_search.rw_beta[tr_search.rw_beta == 0]  = 1 # set true zeros to 1
	tr_search.rw_beta[tr_search.rw_beta > 50]  = 0 # ignores everything less than 100
	rw_beta_sparse = scipy.sparse.csr_matrix(tr_search.rw_beta)
	tr_search.rw_beta = rw_beta_sparse

	return result_df, rw_beta_sparse, tr_search
