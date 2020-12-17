# tabulation 
import os
import pandas as pd
import numpy as np
import re
import scipy.sparse

def tabulate(clone_df1, clone_df2, pwmat = None, cdr3_name = 'cdr3_b_aa', v_gene_name = 'v_b_gene', j_gene_name = 'j_b_gene'):
	"""
	Tabulates frequency and counts of biochemically 
	similar clones in a bulk clone DataFrame <clone_df2>
	to special cloens o interest in <clone_df1>, most likely 
	derived from antigen-enriched repertoires (MIRA or tetramer-based enrichments). 

	It is REQUIRED that clone_df1 contain a 'radius' and 'regex' column, 
	specifiying the maximum search radius per clone and motif constraint.

	Results are provided for 
		(1)             TCRDIST <=  radius 
		(2) (regex_adj) TCRDIST <= radius + CDR3 matches regex
		(3) (tcrdist0)  TCRDIST == 0, (in almost all cases a direct match of both TRV Gene and CDR3 AA)

	Parameters
	----------
	clone_df1 : pd.DataFrame
		contains the metaclonotype definitions (at a minimum: cdr3,cdr2,cdr1,phmc,vgene,radius,regex)
	clone_df2 : pd.DataFrame
		contains the bulk clones (at a minimum: cdr3,cdr2,cdr1,phmc,vgene)
	pwmat : np.ndarray or scipy.sparse.csr.csr_matrix
		
	max_radius : int 
		This refers to the maximum radius used overall in the sparse compression. 
		row-wise all zeros will be converted to this value + 1.
	cdr3_name : str
		must describe CDR3 in both clone_df1 and clone_df2 e.g., 'cdr3_b_aa'
	v_gene_name : str
		must describe TRV gene =in both clone_df1 and clone_df2e.g., 'v_b_gene'
	j_gene_name : str
		must describe TRJ gene in both clone_df1 and clone_df2 e.g., 'j_b_gene'
	
	Returns 
	-------
	result_df : pd.DataFrame

	"""
	assert 'radius' in clone_df1.columns, '<radius> must be a column of clone_df1'
	assert 'regex' in clone_df1.columns, '<radius> must be a column of clone_df1'
	assert isinstance(clone_df1, pd.DataFrame), 'clone_df1 must be a Pandas DataFrame'
	assert isinstance(clone_df2, pd.DataFrame), 'clone_df1 must be a Pandas DataFrame'
	assert isinstance(pwmat,np.ndarray) or isinstance(pwmat, scipy.sparse.csr.csr_matrix)
	assert pwmat.shape[0] == clone_df1.shape[0]
	assert pwmat.shape[1] == clone_df2.shape[0]
	print(f"MATRIX SHAPE {pwmat.shape}")
	print(f"MATRIX TYPE  {type(pwmat)}")
	
	max_radii = clone_df1['radius'].values
	max_radius = clone_df1['radius'].values.max() + 1

	# Each search sequence has a max radius 
	icol = list()
	icol_0 = list()
	idist = list()
	for i,radius in zip(range(pwmat.shape[0]), max_radii):
		# convert a row of sparse matrix to an array, replace zeros with max_radius, replace -1 with 0
		if isinstance(pwmat, scipy.sparse.csr.csr_matrix):
			row_distances = pwmat[i,].toarray()[0] #)[0]
			row_distances[row_distances == 0] = max_radius
			row_distances[row_distances == -1] = 0
		elif isinstance(pwmat, np.ndarray):
			row_distances = pwmat[i,] #)[0]
		else:
			raise TypeError('<pwma> must be a numpy.ndarray or scipy.sparse.csr.csr_matrix')

		columns_within_radius = list(np.where(row_distances <= radius)[0])
		columns_within_radius_0 = list(np.where(row_distances == 0)[0])
		distances_within_radius = list(row_distances[columns_within_radius])
		icol.append(columns_within_radius)
		idist.append(distances_within_radius)
		icol_0.append(columns_within_radius_0)
	
	# Retrieve sequences from the bulk clone_df
	iseqs   = [clone_df2[cdr3_name].iloc[x].to_list() for x in icol]
	ivgenes = [clone_df2[v_gene_name].iloc[x].to_list()  for x in icol]
	ijgenes = [clone_df2[j_gene_name].iloc[x].to_list()  for x in icol]
	# Retrieve abundances from the bulk clone df
	icounts    = [clone_df2['count'].iloc[x].to_list()                 for x in icol]
	ifreqs     = [clone_df2['productive_frequency'].iloc[x].to_list()  for x in icol]
	
	isumcounts    = [np.sum(x) for x in icounts]
	isumfreqs     = [np.sum(x) for x in ifreqs ]

	df_summ = pd.DataFrame({ 'bulk_sum_freq'     : isumfreqs,
							 'bulk_sum_counts'   : isumcounts,
							 'bulk_seqs'      	 : iseqs,
							 'bulk_v_genes'   	 : ivgenes,
							 'bulk_j_genes'   	 : ijgenes,
							 'bulk_distances' 	 : idist,
							 'bulk_counts'    	 : icounts, 
							 'bulk_freqs'     	 : ifreqs})
	result_df = pd.concat([clone_df1.copy(), df_summ], axis = 1)

	regex_on_bulk = [[re.search(pattern = r['regex'], string =  s) for s in r['bulk_seqs']] for _,r in result_df.iterrows()]
	bulk_regex_match = [[True if (x is not None) else False for x in sublist] for sublist in regex_on_bulk]
	result_df['bulk_regex_match'] = bulk_regex_match
	result_df['bulk_sum_freqs_regex_adj'] = [pd.Series(r['bulk_freqs'], dtype = "float64")[pd.Series(r['bulk_regex_match'], dtype = 'bool')].sum() for i,r in result_df.iterrows()]
	result_df['bulk_sum_counts_regex_adj'] = [pd.Series(r['bulk_counts'], dtype = "int32")[pd.Series(r['bulk_regex_match'],  dtype = 'bool')].sum() for i,r in result_df.iterrows()]

	# TCRDIST == 0
	icounts_0 = [clone_df2['count'].iloc[x].to_list() for x in icol_0]
	ifreqs_0 = [clone_df2['productive_frequency'].iloc[x].to_list()  for x in icol_0]
	isumcounts_0 = [np.sum(x) for x in icounts_0]
	isumfreqs_0     = [np.sum(x) for x in ifreqs_0 ]
	result_df['bulk_sum_freqs_tcrdist0'] = isumfreqs_0
	result_df['bulk_sum_counts_tcrdist0'] = isumcounts_0

	return result_df
