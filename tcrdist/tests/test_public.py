"""
test_public.py
"""
import pandas as pd
import numpy as np
import os
import pytest 
from tcrdist.repertoire import TCRrep
from tcrsampler.sampler import TCRsampler
from tcrdist.public import *
from tcrdist.tree import _default_tcrsampler_mouse_beta, _default_tcrsampler_mouse_alpha
import scipy

# Fixture Mouse Beta
df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index().copy()
tr = TCRrep(cell_df = df.head(300).copy(), 
            organism = 'mouse', 
            chains = ['beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv', 
            compute_distances = True)
tr.clone_df['radius'] = 23
tsampler_beta = _default_tcrsampler_mouse_beta()

# Fixture Mouse Alpha
tra = TCRrep(cell_df = df.head(300).copy(), 
            organism = 'mouse', 
            chains = ['alpha'], 
            db_file = 'alphabeta_gammadelta_db.tsv', 
            compute_distances = True)
tra.clone_df['radius'] = 23
tsampler_alpha = _default_tcrsampler_mouse_alpha()


def test_TCRpublic():
	df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index().copy()
	tr = TCRrep(cell_df = df.head(200).copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv', 
	            compute_distances = True)
	tr.clone_df['radius'] = 40
	from tcrdist.public import TCRpublic
	tp = TCRpublic(tcrrep = tr, organism = 'mouse', chain = 'beta')
	tp.report()

def test__neighbors_fixed_radius():
	""" Returns the list of neighbor column indices if within the fixed radius """
	result = _neighbors_fixed_radius(pwmat = tr.pw_beta, radius = 20) 
	assert isinstance(result, list)

def test__K_neighbors_fixed_radius():
	result = _K_neighbors_fixed_radius(pwmat = tr.pw_beta, radius = 20) 
	assert isinstance(result, list)

def test__neighbors_variable_radius():
	""" Returns the list of neighbor column indices if within the fixed radius """
	result = _neighbors_variable_radius(pwmat = tr.pw_beta, radius_list = tr.clone_df['radius'].to_list())
	assert isinstance(result, list)

def test__K_neighbors_variable_radius():
	""" Returns the number of neighbors (self-inclusive) if within the fixed radius"""
	result = _K_neighbors_variable_radius(pwmat = tr.pw_beta, radius_list = tr.clone_df['radius'].to_list())
	assert isinstance(result, list)

def test_neighbors_sparse_fixed_radius():
	M = np.array([[-1,0,1,4],[1,-1,0,2],[0,10,-1,3],[0,0,2,-1]])
	S = scipy.sparse.csr_matrix(M)
	NN = _neighbors_sparse_fixed_radius(csrmat = S, radius = 4)
	assert NN == [[0, 2, 3], [0, 1, 3], [2, 3], [2, 3]]

def test_neighbors_sparse_variable_radius():
	M = np.array([[-1,1,1,4],[1,-1,1,2],[1,10,-1,3],[5,1,2,-1]])
	S = scipy.sparse.csr_matrix(M)
	N = _neighbors_sparse_variable_radius(csrmat = S, radius_list = [1,1,20,1])
	assert N ==  [[0, 1, 2], [0, 1, 2], [0, 1, 2, 3], [1,3]]

def test_make_motif_logo():
	
	svg, svg_raw = make_motif_logo(
			tcrsampler = tsampler_beta,
			clone_df = tr.clone_df,
			pwmat = tr.pw_beta,
			centroid = 'CASSDFDWGGDAETLYF',
			v_gene = 'TRBV13-1*01',
			radius = 24,
			pwmat_str = 'pw_beta',
			cdr3_name = 'cdr3_b_aa',
			v_name = 'v_b_gene',
			gene_names = ['v_b_gene','j_b_gene'])
	assert isinstance(svg, str)
	assert svg.startswith('<svg')
	assert isinstance(svg_raw, str)
	assert svg_raw.startswith('<svg')

# TODO: Test with alpha,beta,gamma,delta # Test with nndif and hdif
def test_quasi_public_meta_clonotypes_beta():
	RB = _quasi_public_meta_clonotypes(clone_df = tr.clone_df.copy(), 
								  pwmat = tr.pw_beta.copy(),
								  tcrsampler = tsampler_beta, 
								  cdr3_name = 'cdr3_b_aa',
								  v_gene_name = 'v_b_gene',
								  nr_filter = True,
								  output_html_name = "quasi_public_clones_beta_testonly.html",
								  labels = ['clone_id',
								  			'cdr3_b_aa', 
								  			'v_b_gene',
								  			'j_b_gene',
								  			'radius',
								  			'neighbors',
								  			'K_neighbors',
								  			#'cdr3s',
								  			'nsubject',
								  			'qpublic',
								  			'cdr3_b_aa.summary',
								  			'v_b_gene.summary',
								  			'j_b_gene.summary',
								  			'subject.summary'],
								  fixed_radius = False,
								  radius = 24,
								  query_str = 'qpublic == True & K_neighbors > 3',
								  kargs_member_summ = {
								 	'key_col'   : 'neighbors', 
								 	'count_col' : 'count',
								 	'addl_cols' : ['subject'],
								 	'addl_n'    : 4},
								  kargs_motif = {
								  	'pwmat_str'  : 'pw_beta',
									'cdr3_name'  : 'cdr3_b_aa',
									'v_name'     : 'v_b_gene',
									'gene_names' : ['v_b_gene','j_b_gene']})



def test_quasi_public_meta_clonotypes_beta_from_nieghbor_diff():
	from tcrdist.rep_diff import neighborhood_diff
	ndif = neighborhood_diff(   clone_df= tr.clone_df, 
								pwmat = tr.pw_beta, 
								count_col = 'count', 
								x_cols = ['epitope'], 
								knn_radius = 25)

	nn_clone_df = pd.concat([tr.clone_df, ndif[['neighbors', 'K_neighbors','val_0','ct_0']] ], axis = 1)
	RB = _quasi_public_meta_clonotypes(clone_df = nn_clone_df, 
							  pwmat = tr.pw_beta.copy(),
							  tcrsampler = tsampler_beta, 
							  cdr3_name = 'cdr3_b_aa',
							  v_gene_name = 'v_b_gene',
							  nr_filter = True,
							  output_html_name = "quasi_public_clones_beta_testonly_from_NN.html",
							  labels = ['clone_id',
							  			'val_0','ct_0', # <- NNDIF COLUMNS
							  			'cdr3_b_aa', 
							  			'v_b_gene',
							  			'j_b_gene',
							  			'radius',
							  			'neighbors',
							  			'K_neighbors',
							  			'nsubject',
							  			'qpublic',
							  			'cdr3_b_aa.summary',
							  			'v_b_gene.summary',
							  			'j_b_gene.summary',
							  			'subject.summary'],
							  fixed_radius = False,
							  radius = 24,
							  query_str = 'qpublic == True & K_neighbors > 3',
							  kargs_member_summ = {
							 	'key_col'   : 'neighbors', 
							 	'count_col' : 'count',
							 	'addl_cols' : ['subject'],
							 	'addl_n'    : 4},
							  kargs_motif = {
							  	'pwmat_str'  : 'pw_beta',
								'cdr3_name'  : 'cdr3_b_aa',
								'v_name'     : 'v_b_gene',
								'gene_names' : ['v_b_gene','j_b_gene']})



def test_quasi_public_meta_clonotypes_alpha():
	RA = _quasi_public_meta_clonotypes(clone_df = tra.clone_df.copy(), 
							  pwmat = tra.pw_alpha.copy(),
							  tcrsampler = tsampler_alpha, 
							  cdr3_name = 'cdr3_a_aa',
							  v_gene_name = 'v_a_gene',
							  nr_filter = True,
							  output_html_name = "quasi_public_clones_alpha_testonly.html",
							  labels = ['clone_id',
							  			'cdr3_a_aa', 
							  			'v_a_gene',
							  			'j_a_gene',
							  			'radius',
							  			'neighbors',
							  			'K_neighbors',
							  			#'cdr3s',
							  			'nsubject',
							  			'qpublic',
							  			'cdr3_a_aa.summary',
							  			'v_a_gene.summary',
							  			'j_a_gene.summary',
							  			'subject.summary'],
							  fixed_radius = False,
							  radius = 24,
							  query_str = 'qpublic == True & K_neighbors > 3',
							  kargs_member_summ = {
							 	'key_col'   : 'neighbors', 
							 	'count_col' : 'count',
							 	'addl_cols' : ['subject'],
							 	'addl_n'    : 4},
							  kargs_motif = {
							  	'pwmat_str'  : 'pw_alpha',
								'cdr3_name'  : 'cdr3_a_aa',
								'v_name'     : 'v_a_gene',
								'gene_names' : ['v_a_gene','j_a_gene']})
	
		

