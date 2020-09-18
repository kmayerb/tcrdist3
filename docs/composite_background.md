"""
composite_backgrounds.md
Effectiveness of Backgrounds.md
"""



```python
import numpy as np
import os
import pandas as pd
import sys
from tcrsampler.sampler import TCRsampler
from tcrdist.repertoire import TCRrep
from tcrdist.automate import auto_pgen
from tcrdist.background import make_gene_usage_counter, make_vj_matched_background, make_flat_vj_background, get_gene_frequencies, calculate_adjustment
```

```python
ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
ts.build_background(stratify_by_subject = True, use_frequency = False)

# Functions for making a background set
ix =[['TRBV19*01', 'TRBJ2-5*01', 3],
['TRBV24-1*01', 'TRBJ2-4*01', 3],
['TRBV25-1*01', 'TRBJ2-4*01', 3],
['TRBV30*01', 'TRBJ2-3*01', 2],
['TRBV5-4*01', 'TRBJ2-3*01', 2],
['TRBV11-2*01', 'TRBJ2-2*01', 2],
['TRBV2*01', 'TRBJ1-5*01', 1],
['TRBV12-5*01', 'TRBJ2-7*01', 1],
['TRBV4-1*01', 'TRBJ1-6*01', 1],
['TRBV6-5*01', 'TRBJ1-6*01', 1],
['TRBV13*01', 'TRBJ2-3*01', 1],
['TRBV18*01', 'TRBJ2-3*01', 1],
['TRBV14*01', 'TRBJ2-7*01', 1],
['TRBV6-6*01', 'TRBJ2-7*01', 1],
['TRBV10-3*01', 'TRBJ2-3*01', 1],
['TRBV7-2*01', 'TRBJ2-1*01', 1],
['TRBV5-1*01', 'TRBJ2-1*01', 1]]

"""
< df_rare > 

For testing we will use a set of 25 TCRs generated from 
rare and semi-rare V,J pairings. We use 25 only 
because we will be comuting distances against 4.6 Million 
seqs.
"""
flatten = lambda l: [item for sublist in l for item in sublist]
df_rare= pd.concat([pd.DataFrame({'cdr3_b_aa' : flatten(ts.sample([[x[0], x[1], x[2]]])) , 'v_b_gene':x[0], 'j_b_gene':x[1]}) for x in ix]).reset_index(drop = True)


""" Make 100,000 a VJ-Matched Background """
gene_usage_counter = make_gene_usage_counter(df_rare)
df_vj_bkgd = make_vj_matched_background(ts = ts,
	gene_usage_counter = gene_usage_counter, 	
	size = 100000, 
	recomb_type="VDJ", 
	chain_folder = "human_T_beta",
	cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])

""" Make 100,000 seqs background with flat uniform distribution of all V,J-gene usages"""
df_vj_flat = make_flat_vj_background(ts = ts, size = 100000)


""" Take 100,000 real seqs at random """
df_real = ts.ref_df.sample(100000, random_state = 1).rename(columns = {'v_reps':'v_b_gene','j_reps' :'j_b_gene', 'cdr3' : 'cdr3_b_aa'})[['v_b_gene','j_b_gene', 'cdr3_b_aa']]
df_real = get_gene_frequencies(ts = ts, df = df_real) 
cols = ['v_b_gene','j_b_gene']; adjcol = 'pVJ'
xdf = df_real[cols].groupby(cols).size().reset_index()
xdf.columns = ['v_b_gene','j_b_gene', 'n']
# For each V,J pairing compute frequency in this reference
xdf['ref_freq'] = xdf['n'] / xdf['n'].sum()
df_real = df_real.merge(xdf, how = "left", on = cols).reset_index()
""" Assign Weights to Real Sequences"""
df_real['weights'] =  pd.Series(df_real[adjcol] / df_real['ref_freq'])




"""
df_bkgd0 VJUNIFORM 
"""
df_bkgd0 = pd.concat([df_vj_flat]).reset_index(drop = True)
df_bkgd0['weights'] = calculate_adjustment(df = df_bkgd0)
""" 
df_bkgd1 VJMATCH  
"""
df_bkgd1 = pd.concat([df_vj_bkgd]).reset_index(drop = True)
df_bkgd1['weights'] = calculate_adjustment(df = df_bkgd1)
""" 
df_bkgd2 VJMATCH + VJUNIFORM 
"""
df_bkgd2 = pd.concat([df_vj_bkgd, df_vj_flat]).reset_index(drop = True)
df_bkgd2['weights'] = calculate_adjustment(df = df_bkgd2)

""" 
df_bkgd3 VJMATCH + VJUNIFORM + REAL
""" 
df_bkgd3 = pd.concat([df_bkgd2[['v_b_gene','j_b_gene', 'cdr3_b_aa', 'weights']], df_real[['v_b_gene','j_b_gene', 'cdr3_b_aa', 'weights']]]).reset_index(drop = True)




"""
FIRST COMPUTE THE TRUTH VALUES BASED RARE TCRS vs. 4.6 Milion
"""
tr_rare = TCRrep(cell_df = df_rare, organism = 'human', chains = ['beta'], store_all_cdr= False, compute_distances = False)
auto_pgen(tr_rare)
# Let's see the proportion of TCRs in a real massive dataset < 30 TCRDIST
df_check = ts.ref_df.rename(columns = {'v_reps':'v_b_gene','j_reps' :'j_b_gene', 'cdr3' : 'cdr3_b_aa'})[['v_b_gene','j_b_gene', 'cdr3_b_aa']].copy()
tr_check  = TCRrep(cell_df = df_check, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
tr_rare = TCRrep(cell_df = df_rare, organism = 'human', chains = ['beta'], store_all_cdr= False)
auto_pgen(tr_rare)

tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_check.clone_df )

tr_rare.clone_df['actual_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
tr_rare.clone_df['actual_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
df_save = tr_rare.clone_df.copy()

tr_rare.rw_beta > 50 = 0

import copy 
tr_rare_copy = copy.deepcopy(tr_rare)
tr_rare.clone_df.to_csv("check.df", index = False)

"""
Compute against backgrounds
"""
def compare_to_background(df_bkgd):

	tr_bkgd = TCRrep(cell_df = df_bkgd, organism = 'human', chains =['beta'], compute_distances = False, store_all_cdr = False)
	tr_rare.compute_rect_distances(df = tr_rare.clone_df, df2 = tr_bkgd.clone_df)
	tr_rare.clone_df['raw_hit_freq'] = np.sum(tr_rare.rw_beta <=30, axis = 1) / tr_rare.rw_beta.shape[1]
	tr_rare.clone_df['raw_hit_count']= np.sum(tr_rare.rw_beta <=30, axis = 1)
	tr_rare.clone_df['raw_bias'] = tr_rare.clone_df['raw_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 
	weights = tr_bkgd.clone_df.weights
	tr_rare.clone_df['adj_hit_count'] =  np.sum((1*(tr_rare.rw_beta <=30) * np.array(weights)[None,:]), axis = 1)
	tr_rare.clone_df['adj_hit_freq'] =  np.sum((1*(tr_rare.rw_beta <=30) *np.array(weights)[None,:]), axis = 1) / tr_rare.rw_beta.shape[1]
	tr_rare.clone_df['adj_bias'] = tr_rare.clone_df['adj_hit_freq'] / tr_rare.clone_df['actual_hit_freq'] 

	df_save2 = tr_rare.clone_df.copy()
	df_save2['actual_hit_freq'] = df_save2['actual_hit_freq'].apply(np.format_float_scientific,precision =1)
	df_save2['raw_hit_freq'] = df_save2['raw_hit_freq'].apply(np.format_float_scientific,precision =1)
	df_save2['adj_hit_freq'] = df_save2['adj_hit_freq'].apply(np.format_float_scientific,precision =1)
	df_save2[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen_cdr3_b_aa', 'raw_hit_count', 'actual_hit_freq','raw_hit_freq','adj_hit_freq','raw_bias', 'adj_bias' ]]

	return df_save2


df0 = compare_to_background(df_bkgd0) # df_bkgd0 VJUNIFORM 
df1 = compare_to_background(df_bkgd1) # df_bkgd1 VMATCH
df2 = compare_to_background(df_bkgd2) # df_bkgd2 VJMATCH + VUNIFORM  
df3 = compare_to_background(df_bkgd3) # df_bkgd3 VJMATCH + VJUNIFORM + REAL

dx = [ 'pgen_cdr3_b_aa','cdr3_b_aa', 'v_b_gene', 'j_b_gene','actual_hit_freq']
cx = ['raw_hit_count', 'adj_hit_freq', 'raw_bias','adj_bias']
pd.concat([df0[dx],df0[cx], df1[cx], df2[cx], df3[cx]],axis = 1, keys = ['DESCRIPTION','VJUNIFORM','VJMATCH',' VJMATCH+VUNIFORM','JMATCH+VUNIFORM+SAMPLE']).to_csv('TABLE1.csv')

```
