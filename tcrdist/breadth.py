"""
breadth.py 

This module concerns the estimation of clonal breadth, whether it be at the 
pathogen, protein, or epitope level. Once meta-clonotype have been defined they 
can be used to search for biochemically similar TCRs in bulk repertoires 
that are likely to share antigen recognition. It is possible that a
single TCR clonotype may be conformant with multiple TCR meta-clonotypes, 
so an accurate estimate of clonal breadth must avoid double counting clonotypes.  

To estimate clonal breadth of antigen-associated TCRs within 
a bulk repertoires with (N) productive clonotypes and (M) total 
productive templates, we use a set of (X) previously defined
antigen-associated meta-clonotypes (defined as a (i) Centroid TRV,CDR3, 
(ii) TCR-specific RADIUS, (iii) MOTIF. 

1. Compute the TCRdist between each centroid TCRij for i {1...i...X} 
and all bulk clones {1...j..M} using rectangular search with the 
tcrdist.repertoires.TCRrep.compute_sparse_rect_distances(), producing 
a sparse distance matrix. 

2. Next perform a long-form tabulation that records the network formed between 
all meta-clonotype centroids and bulk sequences within the specified radius. 
This is performed with the function tcrdist.breadth.long_form_tabulation().

The network is represented as a Pandas DataFrame. Where centroid sequences
are recorded as "cdr3_b_aa", "v_b_gene", "j_b_gene" and the conformant sequence
in the bulk repertoire is "cdr3_b_aa_hit", 'v_b_gene_hit', 'j_b_gene_hit'. 
Crucially there is a column "MOTIF" which indicates whether the CDR3 of 
the hit sequence is conformant with the regular expression in the column 
"regex". 

3. The long-form Pandas DataFrame can then be used as input to the function 
tcrdist.breadth.estimate_breadth_and_depth(). The unit of analysis -- that is, 
whether breadth refers to pathogen, protein, or epitope specific breadth --
can be specified in with the argument <breadth_cols>. Crucially, when 
running tcrdist.breadth.long_form_tabulation() the argument 
<search_cols> must include a column indicating the association between 
a metac-clonotype and a particular 'protein' or 'epitope'  
e.g., ['tag', 'protein', 'epitope', cdr3_b_aa', 'v_b_gene', 'j_b_gene', 
'pgen','regex', 'radius']

An example is provided in the documentation. 
"""
import numpy as np
import os
import pandas as pd
import math
from tcrdist.public import _neighbors_sparse_variable_radius

def long_form_tabulation(
	clone_df1,
	clone_df2,
	csrmat, 
	search_cols = ['tag','cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen','regex', 'radius']):
    '''
    Tabulate metaclonotype-discovered hits in 
    bulk data as a long-form Pandas DataFrame
 
    Parameters
    ----------
    clone_df1 : DataFrame
        Typically TCRrep.clone_df of metaclonotypes 
    clone_df2 : DataFame
        Typically TCRrep.clone_df of bulk sequences
    csrmat : scipy.sparse.csrmat
        Typically TCRrep.rw_beta
    search_cols: list
        Default ['tag','feature','cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen','regex', 'radius']
 
    Returns 
    -------
    long_form_df : DataFrame
        Columns: 
        'tag'                  - Metaclonotype Target Name
        'feature_key'          - Metaclonotype Feature (CDR3X+TRVX+TRJX+22)
        'cdr3_b_aa'            - Metaclonotype CDR3
        'v_b_gene'             - Metaclonotype V
        'j_b_gene'             - Metaclonotype J
        'pgen'.                - Metaclonotype pgen
        'regex'                - Metaclonotype regex
        'radius'               - Metaclonotype radius
        'cdr3_b_aa_hit'        - Discovered match CDR3
        'v_b_gene_hit'         - Discovered match CDR3
        'j_b_gene_hit'         - Discovered match CDR3
        'productive_frequency' - Discovered match sample Frequency
        'count'                - Discovered match sample template count
        'dist'                 - TCRdist to match
        'MOTIF'                - regex match (boolean)
    '''
    from tcrdist.public import _neighbors_sparse_variable_radius
    import re

    clone_df1['RADIUS'] = \
        _neighbors_sparse_variable_radius(
            csrmat      = csrmat, 
            radius_list = clone_df1.radius.to_list())
 
    # Tidy Long Form DataFrame to Record All Hits
    long_form = list()
    for i,r in clone_df1.iterrows():
        if len(r['RADIUS']) > 0:
            hit_features = clone_df2.iloc[r['RADIUS'],:]
            hit_features = hit_features[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'productive_frequency', 'count']]
            
            dists = csrmat[i,r['RADIUS']].toarray()[0].tolist() 
            hit_features['dist'] = dists
 
            search_features = r[search_cols]
            feature_key = search_features['cdr3_b_aa'] + '+' + search_features['v_b_gene'] + '+' + search_features['j_b_gene'] + '+' + str(search_features['radius'])
            
            search_features = search_features[search_cols]
            search_features['feature_key'] = feature_key
 
            hit_features['feature_key']    = feature_key
            
            dx = pd.DataFrame(search_features).\
            transpose().\
            merge( hit_features, how = 'left', left_on = 'feature_key', right_on = 'feature_key', suffixes = ['','_hit'])
            
            long_form.append(dx)
 
    long_form_df = pd.concat(long_form)
 
    def _re(regex, cdr3_b_aa_hit):
        return re.search(pattern = regex, string = cdr3_b_aa_hit) is not None

    long_form_df['MOTIF'] = long_form_df.apply(lambda x: _re(x['regex'], x['cdr3_b_aa_hit']), axis=1)
    
    return(long_form_df)


def estimate_breadth_and_depth(df,
    N,
    M, 
    breadth_cols = ['protein'], 
    clonal_cols = ['cdr3_b_aa_hit', 'v_b_gene_hit', 'j_b_gene_hit'],
    motif = True, 
    exact = False):
    """
    Computes Breadth and Depth of Pathogen, Protein or Epitope specific 
    clones.

	Parameters 
	----------
    df : DataFrame
        The result of long_from_tabulation()
    N : int
        number of unique productive TCR clonotypes 
    M : int 
        number of total sequenced productive TCRs 
    breadth_cols : list 
        ['protein'], the columns over which to compute breadth. 
        The Default ['protein'] will compute breadth for each 
        protein. If None, the breadth estimate will be based on 
        all clones specific to all proteins and epitopes. That 
        is it likely be a pathogen specific breadth.
    clonal_cols : list
        A list of columns that define what is a clone.
        to avoid
    motif : bool 
        If True, df.query("MOTIF == True") is applied, such that only 
        clones that matched a MOTIF are considered in an estimate.
    exact : bool 
        If true df.query("dist <= 0") is applied, such that only clones that 
        exactly matched a search centroid are considered. You probably don't 
        want ot use this except for comparisons sake. (Note that perfect matches
        TCRdist = 0, are designated as 0 in sparse format) 

	Returns 
	-------


    Notes
    -----
    This function is called after the function long_form_tabulation()
    and takes its output as input <df>, N and M can be extracted from 
    N =tr_bulk.clone_df.shape[0] , 
    M= tr_bulk.clone_df['count'].sum()))

    For protein specific breadth, input long-form DataFrame must contain 
    column specifying protein and argument <breadth_cols> should be 
    ['protein'. Similarly for epitope specific breadth <breadth_cols> should be 
    ['protein', 'epitope']. Setting the argument  <breadth_cols> to None
    will result in breadth and depth estimates based on the full set 
    of hits.

    Examples 
    --------

    lf_df = long_form_tabulation(
        clone_df1 = tr_search.clone_df,
        clone_df2 = tr_bulk.clone_df, 
        csrmat = tr_search.rw_beta, 
        search_cols = ['tag', 'protein','protein_coordinate','cdr3_b_aa',
                       'v_b_gene', 'j_b_gene', 'pgen','regex', 'radius'])

    # PROTEIN-specific breadth, RADIUS+MOTIF conformant clones 
    estimate_breadth_and_depth(df=lf_df, 
        breadth_cols = ['protein'], 
        N =tr_bulk.clone_df.shape[0] , 
        M= tr_bulk.clone_df['count'].sum(), 
        motif = True,
        exact = False)

    # PROTEIN-specific breadth, RADIUS conformant clones 
    estimate_breadth_and_depth(df=lf_df, 
        breadth_cols = ['protein'], 
        N =tr_bulk.clone_df.shape[0] , 
        M= tr_bulk.clone_df['count'].sum(), 
        motif = False,
        exact = False))\

    # PROTEIN-specific breadth, EXACT MATCHING clones 
    print(estimate_breadth_and_depth(df=lf_df, 
        breadth_cols = ['protein'], 
        N =tr_bulk.clone_df.shape[0] , 
        M= tr_bulk.clone_df['count'].sum(), 
        motif = True, 
        exact = True)
    
    # PATHOGEN-specific breadth, RADIUS+MOTIF conformant clones 
    print(estimate_breadth_and_depth(df=lf_df, 
        breadth_cols = None, 
        N =tr_bulk.clone_df.shape[0] , 
        M= tr_bulk.clone_df['count'].sum(), 
        motif = True, 
    exact = False))
    """
    if motif:
        df = df.query("MOTIF == True")
    if exact: 
        df = df.query("dist <= 0")

    if breadth_cols is None:
        print("No columns specified, breadth estimates are based on all results")

    # To avoid double counting we only take the first instance of 
    # each unique CDR3,V,G clonotype.
    df_unique = df.\
        sort_values(['tag','dist'], ascending=True).\
        groupby(by = clonal_cols).\
        first()
    # Next we compute the log2(1 + count) for each unique hit, note zeros evaluate to log2(1) = 0
    # and need not be calculated in the sum
    df_unique['log2_count'] = df_unique['count'].apply(lambda x : np.log2(1+x))
    df_unique['clone_count'] = 1 
    # Next we groupby at the level of interest
    if breadth_cols is None:
        df_unique['scope'] = "all eptopes"
        df_report = df_unique.groupby('scope')[['productive_frequency', 'clone_count', 'count','log2_count']].sum()
    else:
        df_report = df_unique.groupby(breadth_cols)[['productive_frequency', 'clone_count', 'count','log2_count']].sum()
    
    df_report['breadth'] = df_report['clone_count'].apply(lambda x : x / N)
    df_report['depth'] = df_report['log2_count'].apply(lambda x : x - np.log2(M))
    df_report['1_over_N'] = 1/N
    df_report['1_over_clone_count'] = 1/df_report['clone_count']
    df_report['SE_breadth'] = df_report['breadth'] * np.sqrt(df_report['1_over_N']+df_report['1_over_clone_count'])
    
    return df_report.sort_values('breadth', ascending= False)



def get_safe_chunk(search_clones, bulk_clones,target = 10**7):
	"""
	This function help pick a chunk size that prevents excessive memory use,
	With two CPU, 10*7 should keep total overall memory demand below 1GB
	"""
	ideal_divisor = (search_clones * bulk_clones) / target
	if ideal_divisor < 1:
		ideal_chunk_size = search_clones
		print(ideal_chunk_size)
	else:
		ideal_chunk_size = math.ceil((search_clones)/ ideal_divisor)
		print(ideal_chunk_size)
	return ideal_chunk_size

