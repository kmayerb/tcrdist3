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

# Fixture Mouse Beta
df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index().copy()

tr = TCRrep(cell_df = df.head(300).copy(), 
            organism = 'mouse', 
            chains = ['beta'], 
            db_file = 'alphabeta_gammadelta_db.tsv', 
            compute_distances = True)
tr.clone_df['radius'] = 23
tsampler_beta = TCRsampler(default_background='olga_mouse_beta_t.sampler.tsv')

# Fixture Mouse Alpha
tra = TCRrep(cell_df = df.head(300).copy(), 
            organism = 'mouse', 
            chains = ['alpha'], 
            db_file = 'alphabeta_gammadelta_db.tsv', 
            compute_distances = True)
tra.clone_df['radius'] = 23
tsampler_alpha = TCRsampler(default_background='ruggiero_mouse_alpha_t.tsv.sampler.tsv')


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
								knn_radius = 50)

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
	
		

