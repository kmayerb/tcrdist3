"""
background.py 

Sept 17, 2020
"""

import itertools
import math
import numpy as np
import os
import pandas as pd
import sys

from collections import Counter
from tcrdist.pgen import OlgaModel
from tcrsampler.sampler import TCRsampler
from tcrdist.repertoire import TCRrep
from tcrdist.automate import auto_pgen
from tcrdist.neighbors import compute_ecdf, bkgd_cntl_nn2

flatten = lambda l: [item for sublist in l for item in sublist]

def make_gene_usage_counter(df, cols = ['v_b_gene', 'j_b_gene']):
    """
    >>>
    """
    gene_usage_list= df[cols].to_dict('split')['data']
    gene_usage_list = [(x[0],x[1]) for x in gene_usage_list]
    gene_usage_counter = Counter(gene_usage_list )
    return(gene_usage_counter)



######################
#### CONSTRUCTION ####
######################

def sim_all_cdr3_gen(n = 100, recomb_type="VDJ", chain_folder = "human_T_beta", cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa']):

    def expand_grid(dct):
        rows = itertools.product(*dct.values())
        return pd.DataFrame.from_records(rows, columns=dct.keys())
    
    omb = OlgaModel(recomb_type=recomb_type, chain_folder = chain_folder)
    
    all_vs = [x for x in omb.pgen_model.V_allele_names if x.endswith('*01')]
    all_js = [x for x in omb.pgen_model.J_allele_names if x.endswith('*01')]
    all_possible_beta = expand_grid(
        {'V': all_vs ,
         'J': all_js }
    )


    find_nones = list()
    results = list()
    
    for i,r in all_possible_beta.iterrows():
        
        e = omb.gen_cdr3s(V = r['V'], J = r['J'], n = n)
        results.append(pd.DataFrame({cols[2] : e,  cols[0]: r['V'], cols[1]:r['J']}))
    
        if e[0] is None:
            find_nones.append( [r['V'],r['J']])
 
    return results, find_nones

# Using 4.6M TCRs; Determine frequencies or subsitute minimum frequences

def get_gene_frequencies(ts, df, cols = ['v_b_gene','j_b_gene']):
    min_pV = np.min(list(ts.v_occur_freq.values()))
    min_pJ = np.min(list(ts.j_occur_freq.values()))
    min_pVJ = np.min(list(ts.vj_occur_freq.values()))
    df['pV'] = df[cols[0]].apply(lambda x : ts.v_occur_freq.get(x, min_pV))
    df['pJ'] = df[cols[1]].apply(lambda x : ts.j_occur_freq.get(x, min_pJ))
    df['pVJ'] = [ts.vj_occur_freq.get((r[cols[0]], r[cols[1]]), min_pVJ) for i,r in df[cols].iterrows()] 
    return df


def calculate_adjustment(df, cols = ['v_b_gene','j_b_gene'], adjcol = 'pVJ'):
    """
    
    Returns
    -------
    pd.Series
    """
    xdf = df[cols].groupby(cols).size().reset_index()
    xdf.columns = ['v_b_gene','j_b_gene', 'n']
    # For each V,J pairing compute frequency in this reference
    xdf['ref_freq'] = xdf['n'] / xdf['n'].sum()
    df = df.merge(xdf, how = "left", on = cols).reset_index()
    return  pd.Series(df[adjcol] / df['ref_freq'])




def make_flat_vj_background(ts = None, n=200 , size =100000  ,  cols = ['v_b_gene','j_b_gene']):
    """
    Parameters
    ----------
    ts : TCRsampler
    
    n : int
        Default 200, number of TCRs to generate using Olga
    cols : list
        Default ['v_b_gene','j_b_gene']

    Returns 
    -------
    df : DataFrame

    Makes a flat background where every V,J pair is equally represeented. 

    Enrichment factors for each pair are based on frequency distribution
    of VJ pairing in a TCRsamler
    """
    if size/n < 135:
        raise ValueError(f"Based on size = {size}, increase to alteast {size/1000} to have sufficient TCRs per VJ pairing") 
    if ts is None:
        from tcrsampler.sampler import TCRsampler
        ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')

    results, find_nones = sim_all_cdr3_gen(n = n)
    dfopt = pd.concat(results)
    dfopt = dfopt[dfopt.cdr3_b_aa.notna()]
    # import numpy as np
    # min_pV = np.min(list(ts.v_occur_freq.values()))
    # min_pJ = np.min(list(ts.j_occur_freq.values()))
    # min_pVJ = np.min(list(ts.vj_occur_freq.values()))
    # dfopt['pV'] = dfopt.v_b_gene.apply(lambda x : ts.v_occur_freq.get(x, min_pV))
    # dfopt['pJ'] = dfopt.j_b_gene.apply(lambda x : ts.j_occur_freq.get(x, min_pJ))
    # dfopt['pVJ'] = [ts.vj_occur_freq.get((r[cols[0]], r[cols[1]]), min_pVJ) for i,r in dfopt[cols].iterrows()] 

    min_n = dfopt.groupby(cols).size().min()
    import math 
    n = math.ceil(size / dfopt.groupby(cols).size().shape[0])
    min_n = min(min_n, n)
    parts = list()
    for i,g in dfopt.groupby(cols):
        parts.append(g.sample(min_n))
    df = pd.concat(parts).reset_index(drop = True)
    #df.to_csv("olga_optimized_human_T_beta.csv", index = False)
    df = get_gene_frequencies(ts = ts, df = df, cols = cols)
    return df



def make_vj_matched_background(
    gene_usage_counter, 
    ts = None, 
    size = 100000, 
    recomb_type="VDJ", 
    chain_folder = "human_T_beta",
    cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa']):
    """
    gene_usage_counter : collections.Counter
    size : int
    recomb_type : str
        Default "VDJ", 
    chain_folder : str
        Default is for human beta "human_T_beta",
    cols : list 
        Default is for beta ['v_b_gene', 'j_b_gene', 'cdr3_b_aa']
    
    Example
    -------
    >>> ix =[['TRBV19*01', 'TRBJ2-5*01', 3],['TRBV24-1*01', 'TRBJ2-4*01', 3]]
    >>> df_rare= pd.concat([pd.DataFrame({'cdr3_b_aa' : flatten(ts.sample([[x[0], x[1], x[2]]])) , 'v_b_gene':x[0], 'j_b_gene':x[1]}) for x in ix]).reset_index(drop = True)
    >>> gene_usage_counter = make_gene_usage_counter(df_rare)
    >>> make_vj_matched_background(gene_usage_counter, size = 10)
          v_b_gene    j_b_gene            cdr3_b_aa        pV        pJ       pVJ
    0  TRBV24-1*01  TRBJ2-4*01      CATPVAGVAKNIQYF  0.011942  0.042163  0.000420
    1  TRBV24-1*01  TRBJ2-4*01       CATSPRGSLSIQYF  0.011942  0.042163  0.000420
    2  TRBV24-1*01  TRBJ2-4*01    CATSDLGGGGIHNIQYF  0.011942  0.042163  0.000420
    3    TRBV19*01  TRBJ2-5*01    CASSISDRGKFSETQYF  0.006788  0.089505  0.000394
    4  TRBV24-1*01  TRBJ2-4*01    CATSDLPARTRENIQYF  0.011942  0.042163  0.000420
    5  TRBV24-1*01  TRBJ2-4*01      CATSDPQGAKNIQYF  0.011942  0.042163  0.000420
    6    TRBV19*01  TRBJ2-5*01  CASSISCGRNLGGQETQYF  0.006788  0.089505  0.000394
    7    TRBV19*01  TRBJ2-5*01    CASSCKPSGGYQETQYF  0.006788  0.089505  0.000394
    8    TRBV19*01  TRBJ2-5*01     CASSSGTSHKLETQYF  0.006788  0.089505  0.000394
    9    TRBV19*01  TRBJ2-5*01          CASSDRETQYF  0.006788  0.089505  0.000394
    """
    
    olga_model_beta = OlgaModel(recomb_type=recomb_type, chain_folder = chain_folder)
    total_seqs = np.sum(list(gene_usage_counter.values()))
    adjust_factor = size / total_seqs
    
    dfs = list()
    adjust_depth = 1
    for k,v in gene_usage_counter.items():
        try:
            cdr3s = olga_model_beta.gen_cdr3s(V = k[0], J = k[1], n = v * math.ceil(adjust_factor))
            df = pd.DataFrame({cols[2]: cdr3s})
            df[cols[0]] = k[0]
            df[cols[1]] = k[1]
            dfs.append(df)
        except AttributeError:
            pass
    
    df = pd.concat(dfs).reset_index(drop = True)
    df = df[df.cdr3_b_aa.notna()][cols]
    if ts is None:
        from tcrsampler.sampler import TCRsampler
        ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
    df = get_gene_frequencies(ts = ts, df = df, cols = cols)
    df = df.reset_index(drop = True)
    return(df)


from tcrsampler.sampler import TCRsampler


def sample_britanova(size = 100000, random_state =24082020):
    df = _get_britanova_human_beta_chord_blood_subject_stratified_background(size = size , random_state =random_state)
    return df

def _get_britanova_human_beta_chord_blood_subject_stratified_background(size = 100000, random_state =24082020):
    """
    Produce a background, stratfied by 8 subjects up to 960,000 TCR clones. 

    Unique TCRs are returned without consideration of their clonal frequency.

    Parameters
    ----------
    size : int 
        Size of background
    random_state : int
        Seed for random. sample
    """
    
    """Check for background file. If not present, download"""
    if not 'britanova_human_beta_t_cb.tsv.sampler.tsv' in TCRsampler.currently_available_backgrounds():
        TCRsampler.download_background_file('britanova_human_beta_t_cb.tsv.sampler.tsv.zip')
    else:
        pass 
        # print("CONGRATS 'britanova_human_beta_t_cb.tsv.sampler.tsv' ALREADY INSTALLED")

    ts = TCRsampler(default_background='britanova_human_beta_t_cb.tsv.sampler.tsv')
    # In [10]: ts.ref_df.subject.value_counts()
    # Out[10]:
    # A5-S18.txt    1073416
    # A5-S17.txt     825507
    # A5-S13.txt     692050
    # A5-S12.txt     573373
    # A5-S16.txt     559980
    # A5-S11.txt     519582
    # A5-S14.txt     302288
    # A5-S15.txt     120302 (NOTE THIS IS THE SMALLED STAMPLE)

    total = size  #100K
    nsubject = 8
    import math
    per_sample = math.ceil( total / nsubject)
    if per_sample > 120000:
        raise ValueError("Size: {size} exceed max size (960000) for valid stratification based on smallest sample")

    samples = []
    for subject_name, subject_df in ts.ref_df.groupby('subject'):
        if subject_name == 'A5-S15.txt':
            samples.append(subject_df.sample(per_sample, 
                replace = False,
                random_state = random_state).copy().reset_index(drop = True))
        else:
            samples.append(subject_df.sample(per_sample, 
                replace = False, 
                random_state = random_state).copy().reset_index(drop = True))

    bitanova_unique_clones_sampled = pd.concat(samples).reset_index(drop = True)
    bitanova_unique_clones_sampled = bitanova_unique_clones_sampled[['v_reps', 'j_reps', 'cdr3']].rename(columns = {'v_reps':'v_b_name', 'j_reps':'j_b_name','cdr3':'cdr3_b_aa'})
    return bitanova_unique_clones_sampled 

