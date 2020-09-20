"""
This is an integration test of the major features in background.py
"""
import os
import pytest

__all__ = ['test_background_generation_in_mira_60'] 

def test_background_generation_toy_example():
	import sys
	import os
	import numpy as np
	import pandas as pd
	from tcrsampler.sampler import TCRsampler
	from tcrdist.background import make_gene_usage_counter, get_gene_frequencies, calculate_adjustment, make_gene_usage_counter
	from tcrdist.background import make_vj_matched_background, make_flat_vj_background
	from tcrdist.background import get_stratified_gene_usage_frequency
	from tcrdist.background import sample_britanova
	"""
	SUPPOSE WE HAVE SOME REPERTOIRE WITH THE FOLLOWING GENE USAGE SPECIFIED BY ix
	< df_target > For testing we will use a set of 25 TCRs generated from rare and semi-rare V,J pairings. We use 25 only 
	because we will be comuting distances against 4.6 Million seqs.
		1. TCRsampler, replacing gene occurance frequencies with subject tratified estimates
		NOTE: with replace = True .vj_occur_freq will now be the stratified value
		2. Make V,J gene usage matched backgound to match usage in df_target
		3. Use a subject-stratifeid random draw from the Britanova Chord Blood Samples
		4. Make V,J gene usage matched backgound to match usage in df_target
	"""
	ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv') 	# 1
	ts = get_stratified_gene_usage_frequency(ts = ts, replace = True)

	ix =[['TRBV19*01', 'TRBJ2-5*01', 3],['TRBV24-1*01', 'TRBJ2-4*01', 3],['TRBV25-1*01', 'TRBJ2-4*01', 3],['TRBV30*01', 'TRBJ2-3*01', 2],['TRBV5-4*01', 'TRBJ2-3*01', 2],['TRBV11-2*01', 'TRBJ2-2*01', 2],['TRBV2*01', 'TRBJ1-5*01', 1],['TRBV12-5*01', 'TRBJ2-7*01', 1],['TRBV4-1*01', 'TRBJ1-6*01', 1],['TRBV6-5*01', 'TRBJ1-6*01', 1],['TRBV13*01', 'TRBJ2-3*01', 1],['TRBV18*01', 'TRBJ2-3*01', 1],['TRBV14*01', 'TRBJ2-7*01', 1],['TRBV6-6*01', 'TRBJ2-7*01', 1],['TRBV10-3*01', 'TRBJ2-3*01', 1],['TRBV7-2*01', 'TRBJ2-1*01', 1],['TRBV5-1*01', 'TRBJ2-1*01', 1]]
	flatten = lambda l: [item for sublist in l for item in sublist]
	df_target = pd.concat([pd.DataFrame({'cdr3_b_aa' : flatten(ts.sample([[x[0], x[1], x[2]]])) , 'v_b_gene':x[0], 'j_b_gene':x[1]}) for x in ix]).reset_index(drop = True)

	gene_usage_counter = make_gene_usage_counter(df_target)  							# 2
	df_vj_bkgd = make_vj_matched_background(ts = ts,
		gene_usage_counter = gene_usage_counter, 	
		size = 101000, # Ask for a few extra as Olga can return none if it makes too many non-productive CDR3s
		recomb_type="VDJ", 
		chain_folder = "human_T_beta",
		cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
	df_vj_bkgd = df_vj_bkgd.sample(100000).reset_index(drop = True)
	df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ")
	df_vj_bkgd['source'] = "vj_matched"

	df_britanova_100K = sample_britanova(size = 100000)									# 3
	df_britanova_100K = get_gene_frequencies(ts = ts, df = df_britanova_100K) 
	df_britanova_100K['weights'] = 1
	df_britanova_100K['source'] = "stratified_random"
	df_bkgd = pd.concat([df_vj_bkgd, df_britanova_100K], axis = 0).\
		reset_index(drop = True) 														# 4 

	assert df_bkgd.shape[0] == 200000
	"""
	Visually inspect the gene_usage between target seqs and vj-matched background
	"""
	df_check_match = pd.concat([df_vj_bkgd.groupby(['v_b_gene', 'j_b_gene']).size()/df_vj_bkgd.shape[0], df_target.groupby(['v_b_gene', 'j_b_gene']).size()/df_target.shape[0]], axis = 1)
	assert np.all(abs(df_check_match[0] - df_check_match[1]) < 0.001)
	return df_bkgd


@pytest.mark.skip(reason="This test documents an example, but is redundant. Skipped in the interest of CI time.")
def test_background_generation_in_mira_60(fn = os.path.join('tcrdist','data','covid19','mira_epitope_60_436_MWSFNPETNI_SFNPETNIL_SMWSFNPET.tcrdist3.csv')):
	import sys
	import os
	import numpy as np
	import pandas as pd
	from tcrsampler.sampler import TCRsampler
	from tcrdist.background import make_gene_usage_counter, get_gene_frequencies, calculate_adjustment, make_gene_usage_counter
	from tcrdist.background import make_vj_matched_background, make_flat_vj_background
	from tcrdist.background import get_stratified_gene_usage_frequency
	from tcrdist.background import sample_britanova
	"""
	SUPPOSE WE HAVE SOME REPERTOIRE WITH THE FOLLOWING GENE USAGE SPECIFIED BY ix
	< df_target > For testing we will use a set of 25 TCRs generated from rare and semi-rare V,J pairings. We use 25 only 
	because we will be comuting distances against 4.6 Million seqs.
		1. TCRsampler, replacing gene occurance frequencies with subject tratified estimates
		NOTE: with replace = True .vj_occur_freq will now be the stratified value
		2. Make V,J gene usage matched backgound to match usage in df_target
		3. Use a subject-stratifeid random draw from the Britanova Chord Blood Samples
		4. Make V,J gene usage matched backgound to match usage in df_target
	"""
	ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv') 	# 1
	ts = get_stratified_gene_usage_frequency(ts = ts, replace = True)

	df_target = pd.read_csv(fn)
	df_target =	df_target[['v_b_gene','j_b_gene','cdr3_b_aa']]

	gene_usage_counter = make_gene_usage_counter(df_target)  							# 2

	df_vj_bkgd = make_vj_matched_background(ts = ts,
		gene_usage_counter = gene_usage_counter, 	
		size = 150000, # Ask for a few extra as Olga can return none if it makes too many non-productive CDR3s
		recomb_type="VDJ", 
		chain_folder = "human_T_beta",
		cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
	df_vj_bkgd = df_vj_bkgd.sample(100000).reset_index(drop = True)
	df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ")
	df_vj_bkgd['source'] = "vj_matched"

	df_britanova_100K = sample_britanova(size = 100000)									# 3
	df_britanova_100K = get_gene_frequencies(ts = ts, df = df_britanova_100K) 
	df_britanova_100K['weights'] = 1
	df_britanova_100K['source'] = "stratified_random"
	df_bkgd = pd.concat([df_vj_bkgd, df_britanova_100K], axis = 0).\
		reset_index(drop = True) 														# 4 

	assert df_bkgd.shape[0] == 200000
	#df_bkgd.
	return df_bkgd



