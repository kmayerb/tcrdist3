import os
import random
import re
import parmap
import timeit
import pandas as pd
import numpy as np 
from tcrdist.repertoire import TCRrep
import scipy.sparse
from tcrdist.regex import _index_to_regex_str, _index_to_seqs 
from tcrdist.neighbors import _todense_row, compute_rate, compute_odds_ratio, bkgd_cntl_nn2, bkgd_cntl_nn3
from scipy.stats import chi2_contingency
from inspect import getfile

def test_bkgd_cntl_nn2_setup():
	random.seed(1)
	
	tcrdist_path        =   os.path.split(getfile(TCRrep))[0]
	
	fn_background       =   'm60_bkgd_test_input.csv'
	fn_background       =   os.path.join(tcrdist_path, 'data', 'covid19', fn_background)
	df_background       =   pd.read_csv(fn_background)
	
	#if there is no weights column, generate some random weights
	if 'weights' not in df_background.columns:
		df_background['weights'] = [random.random() for i in df_background.iloc[:,1]]
	
	tr_background       =   TCRrep( cell_df = df_background.copy(), 
									organism = "human", 
									chains= ['beta'], 
									compute_distances = False)
	
	fn_test             =  'm60_test_input.csv'
	fn_test             =   os.path.join(tcrdist_path, 'data', 'covid19', fn_test)
	df_test             =   pd.read_csv(fn_test)
	df_test             =   df_test[['v_b_gene', 'j_b_gene', 'cdr3_b_aa']]
	tr                  =   TCRrep(	cell_df = df_test.copy(), 
									organism = 'human', 
									chains = ['beta'], 
									db_file = 'alphabeta_gammadelta_db.tsv',
									store_all_cdr = False,
									compute_distances = True)
	
	tr.compute_rect_distances(df = tr.clone_df, 
							  df2 = tr_background.clone_df, 
							  store = False)
	
	assert tr.rw_beta.shape[0] == tr.clone_df.shape[0]
	
	return(tr,tr_background)

def test_bkgd_cntl_nn2(tr,tr_background,ncpus=2):
	centers_df = bkgd_cntl_nn2(	tr = tr, 
								tr_background = tr_background,
								ctrl_bkgd = 10**-6, 
								weights = tr_background.clone_df.weights,
								col = 'cdr3_b_aa',
								ncpus = ncpus,
								thresholds = [x for x in range(0,50,2)], 
								generate_regex = True, 
								test_regex = True)
	return(centers_df)

def test_bkgd_cntl_nn3(tr,tr_background,ncpus=2):	
	centers_df = bkgd_cntl_nn3(	tr = tr, 
								tr_background = tr_background,
								ctrl_bkgd = 10**-6, 
								weights = tr_background.clone_df.weights,
								col = 'cdr3_b_aa',
								ncpus = ncpus,
								thresholds = [x for x in range(0,50,2)], 
								generate_regex = True, 
								test_regex = True)
	return(centers_df)

tr,tr_background = test_bkgd_cntl_nn2_setup()
centers_df_nn2 = test_bkgd_cntl_nn2(tr,tr_background)
centers_df_nn3 = test_bkgd_cntl_nn3(tr,tr_background)

#now compare the output - AssertionError if they don't match
assert(centers_df_nn2.round(6).equals(centers_df_nn3.round(6))) #rounding to dodge floating point error



