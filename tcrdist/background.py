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
import warnings

from collections import Counter
from tcrdist.pgen import OlgaModel
from tcrsampler.sampler import TCRsampler

from tcrdist.automate import auto_pgen


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
    xdf.columns = cols + ['n']
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
    df = df[df[cols[2]].notna()][cols]
    
    if ts is None:
        from tcrsampler.sampler import TCRsampler
        ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
        ts = get_stratified_gene_usage_frequency(ts, replace = True)
    df = get_gene_frequencies(ts = ts, df = df, cols = cols)
    df = df.reset_index(drop = True)
    return(df)


def get_stratified_gene_usage_frequency(ts = None, replace = True):
    """
    MODIFIES A TCRsampler instance with esitmates vj_occur_freq_stratified by subject
    
    Parameters
    ----------
    ts : tcrsampler.sampler.TCRsampler

    replace : bool 
        if True, ts.v_occur_freq is set to ts.v_occur_freq_stratified 
        so other functions will work as befor.

    Returns
    -------
    ts : tcrsampler.sampler.TCRsampler

    """

    if ts is None:
        ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')

    # (1/uniqueTCR_sample_depth) / nsubject
    nsubjects = len(ts.ref_df.subject.value_counts())
    inverse_tcrs_per_subject = (1/ts.ref_df.subject.value_counts()) / nsubjects
    # <weights df> 
    ws_df = pd.DataFrame({'subject': inverse_tcrs_per_subject.index, 'sweight': inverse_tcrs_per_subject}).reset_index(drop = True)
    # left join <ws_df> to provide a subject specific weight 
    df = ts.ref_df.merge(ws_df, how = 'left', on = 'subject').copy()
    # All sweights should sum to 1.0, up to rounding error
    assert np.isclose(df.sweight.sum() ,1.0)

    # SUBJECT STRATIFIED V,J FREQUENCIES
    # For each V,J combo take the weighted sum across all samples
    df_vj_occur_freq = df[['sweight','v_reps','j_reps']].groupby(['v_reps','j_reps']).sum().reset_index().rename(columns = {'sweight': 'pVJ'})
    assert np.isclose(df_vj_occur_freq.pVJ.sum() ,1.0)
    df_vj_occur_freq
    # Covert to a dictionary keyed on (V,J)
    ts.vj_occur_freq_stratified = { (x[0],x[1]): x[2] for x in df_vj_occur_freq.to_dict('split')['data']}

    # SUBJECT STRATIFIED VFREQUENCIES
    df_v_occur_freq = df[['sweight','v_reps']].groupby(['v_reps']).sum().reset_index().rename(columns = {'sweight': 'pV'})
    assert np.isclose(df_v_occur_freq.pV.sum() ,1.0)
    df_v_occur_freq
    # Covert to a dictionary keyed on (V,J)
    ts.v_occur_freq_stratified = { x[0]:x[1] for x in df_v_occur_freq.to_dict('split')['data']}

    # SUBJECT STRATIFIED JFREQUENCIES
    df_j_occur_freq = df[['sweight','j_reps']].groupby(['j_reps']).sum().reset_index().rename(columns = {'sweight': 'pJ'})
    assert np.isclose(df_j_occur_freq.pJ.sum() ,1.0)
    df_j_occur_freq
    # Covert to a dictionary keyed on (V,J)
    ts.j_occur_freq_stratified = { x[0]:x[1] for x in df_j_occur_freq.to_dict('split')['data']}
    
    if replace:
        warnings.warn("REPLACING ts.vj_occur_freq WITH ts.vj_occur_freq_stratified", stacklevel=2)
        warnings.warn("REPLACING ts.v_occur_freq  WITH ts.v_occur_freq_stratified", stacklevel=2)
        warnings.warn("REPLACING ts.j_occur_freq  WITH ts.j_occur_freq_stratified", stacklevel=2)
        ts.vj_occur_freq = ts.vj_occur_freq_stratified 
        ts.v_occur_freq  = ts.v_occur_freq_stratified 
        ts.j_occur_freq  = ts.j_occur_freq_stratified 

    return ts 


def sample_britanova(size = 100000, random_state =24082020):
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
    ts = get_stratified_gene_usage_frequency(ts = ts, replace = True)
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
    bitanova_unique_clones_sampled = bitanova_unique_clones_sampled[['v_reps', 'j_reps', 'cdr3']].rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene','cdr3':'cdr3_b_aa'})
    return bitanova_unique_clones_sampled 


def _synthesize_human_beta_vj_background(ts,fn = None, df = None):
    """
    _build_vj_background

    Parameters
    ----------
    ts: tcrsampler.TCRsampler() 
        a TCRsampler instance, with gene usage frequencies (ideally computed get_stratified_gene_usage_frequency()
    fn: str
        file path to MIRA set of TCRs
    df : pandas DataFrame
        MIRA set of TCRs
    Returns
    -------
    df_vj_bkgd : Pandas DataFrame
        A set of background TCRs with the same V and J gene usage as the input set. 
        These are generated using OLGA (Sethna et al.)
    """
    if fn is not None and df is not None:
        raise ValueError("_build_vj_background can accept <df> or <fn> arguments but not both")
    if fn is not None:
        # Load a set set of TCRs.
        df_target = pd.read_csv(fn)
    if df is not None:
        df_target = df.copy()
    # Subset columns
    df_target = df_target[['v_b_gene','j_b_gene','cdr3_b_aa']]
    # Make a gene usage counter
    gene_usage_counter = make_gene_usage_counter(df_target)                             # 2
    print("MAKING A V-GENE, J-GENE MATCHED BACKGROUND.")
    # Check that sampler is using sample stratified frequencies. 
    assert ts.v_occur_freq is ts.v_occur_freq_stratified
    print("USING STRATIFIED FREQUENCIES.")
    # Make a vj matched background. 
    #   Note: <size> aregument should be greater than desired, because Olga can return none due to non-productive CDR3s.
    df_vj_bkgd = make_vj_matched_background(ts = ts,
        gene_usage_counter = gene_usage_counter,    
        size = 150000, 
        recomb_type="VDJ", 
        chain_folder = "human_T_beta",
        cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
    # Sample to get the desired number of TCRs from teh v,j matched set
    df_vj_bkgd = df_vj_bkgd.sample(100000, random_state = 1).reset_index(drop = True)
    print("CALCULATE INVERSE PROBABILITY WEIGHT ADJUSTMENT.")
    # Calculate the invese weighting adjustmetn
    df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ")
    df_vj_bkgd['source'] = "vj_matched"
    # Combine
    return df_vj_bkgd

def _synthesize_human_alpha_vj_background(ts,fn = None, df = None):
    """
    _build_vj_background

    Parameters
    ----------
    ts: tcrsampler.TCRsampler() 
        a TCRsampler instance, with gene usage frequencies (ideally computed get_stratified_gene_usage_frequency()
    fn: str
        file path to MIRA set of TCRs
    df : pandas DataFrame
        MIRA set of TCRs
    Returns
    -------
    df_vj_bkgd : Pandas DataFrame
        A set of background TCRs with the same V and J gene usage as the input set. 
        These are generated using OLGA (Sethna et al.)
    """
    if fn is not None and df is not None:
        raise ValueError("_build_vj_background can accept <df> or <fn> arguments but not both")
    if fn is not None:
        # Load a set set of TCRs.
        df_target = pd.read_csv(fn)
    if df is not None:
        df_target = df.copy()
    # Subset columns
    df_target = df_target[['v_a_gene','j_a_gene','cdr3_a_aa']]
    # Make a gene usage counter
    gene_usage_counter = make_gene_usage_counter(df_target, cols = ['v_a_gene', 'j_a_gene'])                             # 2
    print("MAKING A V-GENE, J-GENE MATCHED BACKGROUND.")
    # Check that sampler is using sample stratified frequencies. 
    assert ts.v_occur_freq is ts.v_occur_freq_stratified
    print("USING STRATIFIED FREQUENCIES.")
    # Make a vj matched background. 
    #   Note: <size> aregument should be greater than desired, because Olga can return none due to non-productive CDR3s.
    df_vj_bkgd = make_vj_matched_background(ts = ts,
        gene_usage_counter = gene_usage_counter,    
        size = 150000, 
        recomb_type="VJ", 
        chain_folder = "human_T_alpha",
        cols = ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'])
    # Sample to get the desired number of TCRs from teh v,j matched set
    df_vj_bkgd = df_vj_bkgd.sample(100000, random_state = 1).reset_index(drop = True)
    print("CALCULATE INVERSE PROBABILITY WEIGHT ADJUSTMENT.")
    # Calculate the invese weighting adjustmetn
    df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ", cols = ['v_a_gene','j_a_gene'])
    df_vj_bkgd['source'] = "vj_matched"
    # Combine
    return df_vj_bkgd


def _synthesize_mouse_beta_vj_background(ts, fn = None, df = None):
    """
    _build_vj_background

    Parameters
    ----------
    ts: tcrsampler.TCRsampler() 
        a TCRsampler instance, with gene usage frequencies (ideally computed get_stratified_gene_usage_frequency()
    fn: str
        file path to MIRA set of TCRs
    df : pandas DataFrame
        MIRA set of TCRs
    Returns
    -------
    df_vj_bkgd : Pandas DataFrame
        A set of background TCRs with the same V and J gene usage as the input set. 
        These are generated using OLGA (Sethna et al.)
    """
    if fn is not None and df is not None:
        raise ValueError("_build_vj_background can accept <df> or <fn> arguments but not both")
    if fn is not None:
        # Load a set set of TCRs.
        df_target = pd.read_csv(fn)
    if df is not None:
        df_target = df.copy()
    # Subset columns
    df_target = df_target[['v_b_gene','j_b_gene','cdr3_b_aa']]
    # Make a gene usage counter
    gene_usage_counter = make_gene_usage_counter(df_target, cols = ['v_b_gene', 'j_b_gene'])                             # 2
    print("MAKING A V-GENE, J-GENE MATCHED BACKGROUND.")
    # Check that sampler is using sample stratified frequencies. 
    assert ts.v_occur_freq is ts.v_occur_freq_stratified
    print("USING STRATIFIED FREQUENCIES.")
    # Make a vj matched background. 
    #   Note: <size> aregument should be greater than desired, because Olga can return none due to non-productive CDR3s.
    df_vj_bkgd = make_vj_matched_background(ts = ts,
        gene_usage_counter = gene_usage_counter,    
        size = 150000, 
        recomb_type="VDJ", 
        chain_folder = "mouse_T_beta",
        cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
    # Sample to get the desired number of TCRs from teh v,j matched set
    df_vj_bkgd = df_vj_bkgd.sample(100000, random_state = 1).reset_index(drop = True)
    print("CALCULATE INVERSE PROBABILITY WEIGHT ADJUSTMENT.")
    # Calculate the invese weighting adjustmetn
    df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ")
    df_vj_bkgd['source'] = "vj_matched"
    # Combine
    return df_vj_bkgd

def _synthesize_mouse_alpha_vj_background(ts, fn = None, df = None):
    """
    _build_vj_background

    Parameters
    ----------
    ts: tcrsampler.TCRsampler() 
        a TCRsampler instance, with gene usage frequencies (ideally computed get_stratified_gene_usage_frequency()
    fn: str
        file path to MIRA set of TCRs
    df : pandas DataFrame
        MIRA set of TCRs
    Returns
    -------
    df_vj_bkgd : Pandas DataFrame
        A set of background TCRs with the same V and J gene usage as the input set. 
        These are generated using OLGA (Sethna et al.)
    """
    if fn is not None and df is not None:
        raise ValueError("_build_vj_background can accept <df> or <fn> arguments but not both")
    if fn is not None:
        # Load a set set of TCRs.
        df_target = pd.read_csv(fn)
    if df is not None:
        df_target = df.copy()
    # Subset columns
    df_target = df_target[['v_a_gene','j_a_gene','cdr3_a_aa']]
    # Make a gene usage counter
    gene_usage_counter = make_gene_usage_counter(df_target, cols = ['v_a_gene', 'j_a_gene'] )                             # 2
    print("MAKING A V-GENE, J-GENE MATCHED BACKGROUND.")
    # Check that sampler is using sample stratified frequencies. 
    assert ts.v_occur_freq is ts.v_occur_freq_stratified
    print("USING STRATIFIED FREQUENCIES.")
    # Make a vj matched background. 
    #   Note: <size> aregument should be greater than desired, because Olga can return none due to non-productive CDR3s.
    df_vj_bkgd = make_vj_matched_background(ts = ts,
        gene_usage_counter = gene_usage_counter,    
        size = 150000, 
        recomb_type="VJ", 
        chain_folder = "mouse_T_alpha",
        cols = ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'])
    # Sample to get the desired number of TCRs from teh v,j matched set
    df_vj_bkgd = df_vj_bkgd.sample(100000, random_state = 1, replace = True).reset_index(drop = True)
    print("CALCULATE INVERSE PROBABILITY WEIGHT ADJUSTMENT.")
    print(df_vj_bkgd)
    # Calculate the invese weighting adjustmetn
    df_vj_bkgd['weights'] = calculate_adjustment(df = df_vj_bkgd, adjcol = "pVJ", cols = ['v_a_gene','j_a_gene'])
    df_vj_bkgd['source'] = "vj_matched"
    # Combine
    return df_vj_bkgd
