import logging
from collections import Counter
import os
import pandas as pd
import numpy as np
import pwseqdist as pw
from tcrdist.rep_funcs import _pw, _pws

import scipy
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

from tcrdist.repertoire import TCRrep
from tcrdist.swap_gene_name import adaptive_to_imgt

__all__ = [ 'import_adaptive',
            'ispublic',
            'simple_cluster_index',
            'cluster_index_to_df',
            'get_centroid_seq',
            'default_tcrdist_and_image']


def import_adaptive_file(   adaptive_filename, 
                            organism = "human", 
                            chain = "beta",
                            return_valid_cdr3_only = True, 
                            count = 'productive_frequency',
                            version_year= 2020,
                            sep = "\t",
                            subject = None,
                            epitope = None,
                            log = True, 
                            swap_imgt_dictionary = None,
                            additional_cols = None):
    """
    Prepare tcrdist3 input from 2020 Adaptive File containing 'bio_identity', 'productive_frequency', 'templates', and 'rearrangement'.

    Parameters 
    ----------
    adaptive_filename : str
        path to the Adaptive filename 
    version : int
        version_year
    epitope : str or None
        name of epitope if known
    subject : str or None
        If none the filename will be used as the subject
    use_as_count : str
        name of column to be used as count (could be 'productive_frequency' or 'templates')
    sep : str
        seperatore in Adaptive file
    organism : str
        'human' or 'mouse'
    chain : str
        'beta' or 'alpha'
    log : bool
        If True, write a log.
    swap_imgt_dictionary : dict or None
        If None, the default dictionary adaptive_to_imgt is used
    additional_cols : None or List
        list of any additional columns you want to keep
    
    Returns 
    -------
    bulk_df : pd.DataFrame
    """
    try:
        bulk_df = pd.read_csv(adaptive_filename, sep= sep, usecols = ['bio_identity', 'productive_frequency', 'templates', 'rearrangement'])
    except ValueError as e:
        raise Exception('Bulk Adpative TCR file was missing required columns') from e


    if swap_imgt_dictionary is None:
        swap_imgt_dictionary = adaptive_to_imgt

    if log: logging.basicConfig(filename='prepbulk.log',level=logging.DEBUG, format='tcrdist3:%(asctime)s\n\t %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    item_names = {'alpha': ["cdr3_b_aa","v_a_gene", "j_a_gene","cdr3_b_nucseq"], 
                  'beta' : ["cdr3_b_aa","v_b_gene","j_b_gene","cdr3_b_nucseq"]}[chain]

    # Coerce strings to numeric
    bulk_df['productive_frequency'] = pd.to_numeric(bulk_df['productive_frequency'],errors='coerce')

    # Parse bio-identity
    ns= {0:"cdr3_aa", 1:"v_gene", 2:"j_gene"}
    cdr_v_j = bulk_df['bio_identity'].str.split("+", expand = True).\
        rename(columns = lambda x: ns[x])
    bulk_df[[item_names[0], 'v_gene', 'j_gene']] = cdr_v_j

    # Convert Names from Adapative to IMGT
    bulk_df[item_names[1]] = bulk_df['v_gene'].apply(lambda x : swap_imgt_dictionary[organism].get(x))
    bulk_df[item_names[2]] = bulk_df['j_gene'].apply(lambda x : swap_imgt_dictionary[organism].get(x))
    # Record Invalid Names
    invalid_v_names =  Counter(bulk_df['v_gene'][  bulk_df[item_names[1]].isna() ].to_list())
    invalid_j_names =  Counter(bulk_df['j_gene'][  bulk_df[item_names[2]].isna() ].to_list())

     # Validate CDR sequences
    bulk_df['valid_cdr3'] = bulk_df[item_names[0]].apply(lambda cdr3: _valid_cdr3(cdr3)) 
    # Count number of valid seqs
    valid = np.sum(bulk_df['valid_cdr3'])
    
    # Assign subject baesd on the < subject > argument 
    if subject is None:
        bulk_df['subject'] = adaptive_filename
    else: 
        bulk_df['subject'] = subject

    # Assign a user supplied or blank epitope baesd on the < epitope > argument 
    if epitope is None:
        bulk_df['epitope'] = 'X' 
    else: 
        bulk_df['epitope'] = epitope

    if additional_cols is None:
        bulk_df = bulk_df[['subject','productive_frequency', 'templates','epitope',item_names[0],item_names[1],item_names[2],'valid_cdr3','rearrangement']].copy()
        bulk_df = bulk_df.rename(columns = {'rearrangement':item_names[3]})
    else:
        selections = ['subject','productive_frequency', 'templates','epitope',item_names[0],item_names[1],item_names[2],'valid_cdr3','rearrangement'] + additional_cols
        bulk_df = bulk_df[selections].copy()
        bulk_df = bulk_df.rename(columns = {'rearrangement':item_names[3]})

    # Logging
    if return_valid_cdr3_only:
        bulk_df = bulk_df[bulk_df['valid_cdr3']]
        if log: logging.info(f"VALID CDR3 ({valid }) / ({bulk_df.shape[0]}) CLONES")
        if log: logging.info(f"OMITTING INVALID CDR3s FROM FINAL RESULTS")
        if log: logging.info(bulk_df[bulk_df.valid_cdr3 == False][['subject', "cdr3_b_aa"]])
        if log: logging.info(f"Final Reults has Column Names {bulk_df.columns}")
        if log: logging.info(f"Invalid ADAPTIVE V-gene names {invalid_v_names} not included\n")
        if log: logging.info(f"Invalid ADAPTIVE J-gene names {invalid_j_names} not included\n")
    
    # Assign count the productive_frequency based on the < use_as_count > argument
    bulk_df['count'] = bulk_df[count].copy()

    return bulk_df


def _valid_cdr3(cdr3):
    """ Return True iff all amino acids are part of standard amino acid list"""
    if not isinstance(cdr3, str):
        return False
    else:
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        valid = np.all([aa in amino_acids for aa in cdr3])
        return valid


def ispublic(gr, var = "subject", n = 1):
    """
    Return True if a cluster public, defined as comprised of members from multiple individuals 
    or cell subsets (e.g., CD4/CD8) 
    
    Parameters
    ----------
    gr : group
        within pandas Data.Frame.groupby
    var : str
        variable name of class that group most transcend to be considered public
    m : int
        number of unique values of selected variable to be considered public
    
    Returns
    -------
    r : bool
        True if a cluster public
    
    """
    r = len(gr[var].value_counts()) > n
    if r:
        return 'public',len(gr[var].value_counts())
    else:
        return 'private',len(gr[var].value_counts())


def simple_cluster_index(   
                            pw_distances = None,
                            method = 'ward',
                            criterion = "distance",
                            t = 75):
    """
    Get 'cluster_index' 

    Parameters
    ----------

    t : int
        scipy.cluster.hierarchy.fcluster param t
    criterion : str 
        scipy.cluster.hierarchy.fcluster param criterion 
    method : str 
        scipy.cluster.linkage param method 
    t : int
        scipy.cluster.hierarcy.fcluster param t

    Notes
    -----
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html
    """

    pw_distances
    
    compressed_dmat = scipy.spatial.distance.squareform(pw_distances, force = "vector")
    Z = linkage(compressed_dmat, method = method)
    cluster_index = fcluster(Z, t = t, criterion = criterion)
    return cluster_index


def cluster_index_to_df(cluster_index):
    """
    Parameters
    ----------
    cluster_index : np.ndarray 

    Returns
    -------
    cluster_df : pd.DataFrame

    Notes
    -----
    
    cluster_df format:
    
    cluster_id                                          neighbors  K_neighbors
        4  [16, 25, 26, 29, 32, 50, 61, 68, 69, 94, 103, ...           24
        92  [35, 38, 41, 105, 131, 146, 181, 186, 189, 206...           18

    """
    dl = dict()
    for k,v in enumerate(cluster_index):
        dl.setdefault(v, [])
        dl[v].append(k)

    cluster_df = pd.DataFrame({'neighbors' : pd.Series(dl)}).sort_index().reset_index().rename(columns = {'index':'cluster_id'})
    cluster_df['K_neighbors'] = cluster_df.neighbors.str.len()
    cluster_df = cluster_df.sort_values(by = 'K_neighbors', ascending = False)
    
    return cluster_df


def get_centroid_seq(df, cdr3_name = 'cdr3_b_aa'):
    """
    Given a list of sequences, returns the sequence with the minimum 
    sum of distances to all other seqs in the list.

    Parameters
    ----------
    seqs : list
        list of strings (amino acid rep)
    metric : func
        defaults to pw.metrics.nw_hamming_metric

    Returns
    -------
    centroid_seq : str

    Example 
    -------
    >>> seqs = ['CASSEILAALGTQYF', 'CASSWTSRETQYF', 'CASSLAQETQYF', 'CASSLAPGDVSQYF', 'CASSWDQETQYF', 'CASSLWWDSGANVLTF', 'CASSLARTLSSGANVLTF', 'CASIPGTLFTFSGANVLTF', 'CASSFASSGANVLTF', 'CASSYRLLSGANVLTF']	
    >>> get_centroid_seq(seqs)
    'CASSFASSGANVLTF'

    Notes 
    -----
    In case of multiple occurrences of the minimum values, the indices 
    corresponding to the first occurrence are returned.
    """
    #import pwseqdist as pw
    #from scipy.spatial.distance import squareform
    seqs = df[cdr3_name]
    if len(seqs) < 3:
        return df.head(1)[cdr3_name], None, None, None

    metrics = {
            "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
            "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
            "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
            "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist}

    # Define weights
    weights = { "cdr3_b_aa" : 3,
                "pmhc_b_aa" : 1,
                "cdr2_b_aa" : 1,
                "cdr1_b_aa" : 1}

    kargs = {
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
    
    dmat = _pws(df = df,
        metrics = metrics, 
        weights = weights, 
        store = False,
        uniquify=False,
        kargs = kargs)
        
    dmat = dmat['tcrdist']

    dmat = dmat.astype(int)
    iloc_idx= dmat.sum(axis = 0).argmin()
    centroid_seq  = df[cdr3_name].to_list()[iloc_idx]
    loc_idx = df.index.to_list()[iloc_idx]
    return centroid_seq, dmat, iloc_idx, loc_idx

def get_centroid_seq_alpha(df, cdr3_name = 'cdr3_a_aa'):
    """
    Given a list of sequences, returns the sequence with the minimum 
    sum of distances to all other seqs in the list.

    Parameters
    ----------
    seqs : list
        list of strings (amino acid rep)
    metric : func
        defaults to pw.metrics.nw_hamming_metric

    Returns
    -------
    centroid_seq : str

    Example 
    -------
    >>> seqs = ['CASSEILAALGTQYF', 'CASSWTSRETQYF', 'CASSLAQETQYF', 'CASSLAPGDVSQYF', 'CASSWDQETQYF', 'CASSLWWDSGANVLTF', 'CASSLARTLSSGANVLTF', 'CASIPGTLFTFSGANVLTF', 'CASSFASSGANVLTF', 'CASSYRLLSGANVLTF']	
    >>> get_centroid_seq(seqs)
    'CASSFASSGANVLTF'

    Notes 
    -----
    In case of multiple occurrences of the minimum values, the indices 
    corresponding to the first occurrence are returned.
    """
    #import pwseqdist as pw
    #from scipy.spatial.distance import squareform
    seqs = df[cdr3_name]
    if len(seqs) < 3:
        return df.head(1)[cdr3_name], None, None, None

    metrics = {
            "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
            "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
            "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
            "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist}

    # Define weights
    weights = { "cdr3_a_aa" : 3,
                "pmhc_a_aa" : 1,
                "cdr2_a_aa" : 1,
                "cdr1_a_aa" : 1}

    kargs = {
                "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}
    
    dmat = _pws(df = df,
        metrics = metrics, 
        weights = weights, 
        store = False,
        uniquify=False,
        kargs = kargs)
        
    dmat = dmat['tcrdist']

    dmat = dmat.astype(int)
    iloc_idx= dmat.sum(axis = 0).argmin()
    centroid_seq  = df[cdr3_name].to_list()[iloc_idx]
    loc_idx = df.index.to_list()[iloc_idx]
    return centroid_seq, dmat, iloc_idx, loc_idx



def bulk_adaptive_dataset_to_tcrdist3_clone_df( bulk_filename = None,
                                                minimum_file = None,
                                                maximum_file = None,
                                                organism = 'human',
                                                chains = ['beta'],
                                                epitope = None):
    """
    bulk_filename : str 
        input adaptive file e.g., "KHBR20-00150_TCRB.tsv"
    minimum_file : str
        minimum (CDR3,V,J,freq  only) output path e.g.,"KHBR20-00150_TCRB.tsv.tcrdist3.v_min.csv",
    maximum_file : str
        maximum (All CDRs, V,J, freq, subject,) output path "KHBR20-00150_TCRB.tsv.tcrdist3.v_max.csv",
    organism = 'human',
    chains = ['beta'],
    epitope = None)
    Bulk adaptive (Beta, Human) data to a tcrdist3 clone DataFrame. 

    2020 adaptive data has column bio_identity, productive_frequency and template. 

    * We will use productive frequency for counts.
    * We will only consider CDR3 with all valid amino acids
    """

    log = True
    logging.basicConfig(filename='prepbulk.log',level=logging.DEBUG, format='KMB:%(asctime)s\n\t %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Load only three columns 
    try:
        bulk_df = pd.read_csv(bulk_filename, sep= "\t", usecols = ['bio_identity', 'productive_frequency', 'templates'])
    except ValueError as e:
        raise Exception('Bulk Adpative TCR file was missing required columns') from e

    # Coerce strings to numeric
    bulk_df['productive_frequency'] = pd.to_numeric(bulk_df['productive_frequency'],errors='coerce')

    # Parse bio-identity
    ns= {0:"cdr3_b_aa", 1:"v_gene", 2:"j_gene"}
    cdr_v_j = bulk_df['bio_identity'].str.split("+", expand = True).\
        rename(columns = lambda x: ns[x])
    bulk_df[["cdr3_b_aa","v_gene","j_gene"]] = cdr_v_j

    # Convert Names from Adapative to IMGT
    bulk_df['v_b_gene'] = bulk_df['v_gene'].apply(lambda x : adaptive_to_imgt['human'].get(x))
    bulk_df['j_b_gene'] = bulk_df['j_gene'].apply(lambda x : adaptive_to_imgt['human'].get(x))

     # Validate CDR sequences
    bulk_df['valid_cdr3'] = bulk_df['cdr3_b_aa'].apply(lambda cdr3: _valid_cdr3(cdr3)) 
    # Count number of valid seqs
    valid = np.sum(bulk_df['valid_cdr3'])
    
    # Assign subject the filename
    bulk_df['subject'] = bulk_filename
    # Assign a blank epitope
    if epitope is None:
        bulk_df['epitope'] = 'X' 
    else: 
        bulk_df['epitope'] = epitope

    # Report % Valid
    print(f"VALID CDR3 ({valid }) / ({bulk_df.shape[0]}) CLONES")
    if log: logging.info(f"VALID CDR3 ({valid }) / ({bulk_df.shape[0]}) CLONES")
    print(f"OMITTING INVALID CDR3s FROM FINAL RESULTS")
    if log: logging.info(f"OMITTING INVALID CDR3s FROM FINAL RESULTS")
    # Log what was dropped
    print(bulk_df[bulk_df.valid_cdr3 == False][['subject', "cdr3_b_aa"]])
    if log: logging.info(bulk_df[bulk_df.valid_cdr3 == False][['subject', "cdr3_b_aa"]])
    print(f"Final Reults has Column Names {bulk_df.columns}")
    if log: logging.info(f"Final Reults has Column Names {bulk_df.columns}")

    bulk_df = bulk_df[['subject','productive_frequency', 'templates','epitope','cdr3_b_aa','v_b_gene','j_b_gene','valid_cdr3']].copy()
    bulk_df = bulk_df[bulk_df['valid_cdr3']]
    
    # Asign Count the productive_frequency
    bulk_df['count'] = bulk_df['productive_frequency'].copy()

    tr = TCRrep(cell_df = bulk_df, 
                organism = 'human', 
                chains = ['beta'], 
                compute_distances = False,
                infer_index_cols = True,
                deduplicate=False,
                cpus = 1,
                db_file = 'alphabeta_gammadelta_db.tsv')

    tr.clone_df[['cdr3_b_aa','v_b_gene','j_b_gene','productive_frequency']].to_csv(minimum_file, sep = "\t", index = False)
    tr.clone_df.to_csv(maximum_file, sep = "\t", index = False)
    
    return tr.clone_df.copy()

def default_dist_clust_centroids(infile, cpus = 1, cdr3_b_aa_weight = 5, max_dist = 200):
    
    from tcrdist.repertoire import TCRrep
    
    print(infile)
    
    df = pd.read_csv(infile)
    df['count'] = 1

    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['beta'], 
                compute_distances = False,
                infer_index_cols = True,
                deduplicate=False,
                cpus = cpus,
                store_all_cdr=False,
                db_file = 'alphabeta_gammadelta_db.tsv')
    
    # Overweight CDR3B
    tr.weights_b['cdr3_b_aa'] = cdr3_b_aa_weight 

    icols = ['cell_type','subject','v_b_gene','j_b_gene','cdr3_b_aa','cdr3_b_nucseq','cdr1_b_aa','cdr2_b_aa','pmhc_b_aa']
    # Manually assign index cols    
    tr.index_cols = icols
    # Deduplicate
    tr.deduplicate()
    # Compute Distances
    tr.compute_distances()
    
    
    # Cluster based on the max_dist
    ci = simple_cluster_index(tr.pw_beta, t = max_dist)

    # Get a DataFrame for clusters (It is returned largest to smallest)
    ci_df = cluster_index_to_df(cluster_index = ci)  

    # Determine degree of public sharing in each cluster
    publicities = list()
    for i,r in ci_df.iterrows():
        clone_cluster_df = tr.clone_df.iloc[r['neighbors'],]
        publicity,n_subjects = ispublic(clone_cluster_df)
        publicities.append((publicity,n_subjects))
    # Store public status in a DataFrame
    public_df = pd.DataFrame(publicities).rename(columns = {0:'public',1:'n_subjects'})

    # Iterate through the ci_df DataFrame
    counter = 0
    centroids = list()
    cluster_df_list = list()
    for i,r in ci_df.iterrows():
        counter = counter + 1
        #print(r['neighbors'])
        clone_cluster_df = tr.clone_df.iloc[r['neighbors'],]
        cluster_df_list.append(clone_cluster_df )
        publicity = ispublic(clone_cluster_df)
        try:
            centroid, dmatrix, iloc_ind, loc_ind, = get_centroid_seq(df = tr.clone_df.iloc[r['neighbors'],])
            #print(centroid)
            #print(dmatrix)
            #print(np.max(dmatrix))
            #print(clone_cluster_df.v_b_gene.value_counts())
            #print(clone_cluster_df.j_b_gene.value_counts())
            #print(clone_cluster_df.subject.value_counts())
            #print(clone_cluster_df)
            #print("----")
            #print(clone_cluster_df.iloc[iloc_ind,])
            #print(centroid)
            assert clone_cluster_df.iloc[iloc_ind,]['cdr3_b_aa'] == centroid
            #print(loc_ind)
            #print(tr.clone_df.iloc[loc_ind,])
            assert tr.clone_df.iloc[loc_ind,]['cdr3_b_aa'] == centroid
            # store the centroid
            centroids.append( clone_cluster_df.iloc[iloc_ind,].reset_index(drop=True).copy())
        except Exception as e:
            print(e)
            print("fail")
            centroids.append( clone_cluster_df.iloc[0,].reset_index(drop=True).copy())

    # Save the Information in a centroids_df

    # renames --- {0: 'cell_type', 1: 'subject', 2: 'v_b_gene', 3: 'j_b_gene', 4: 'cdr3_b_aa', 5: 'cdr3_b_nucs, 6: 'cdr1_b_aa', 7: 'cdr2_b_aa', 8: 'pmhc_b_aa', 9: 'count', 10: 'clone_id'}
    renames = {i:v for i,v in enumerate(tr.clone_df.columns.to_list())}
    centroids_df = pd.DataFrame(centroids).rename(columns = renames)
    centroids_df['neighbors']     = ci_df['neighbors'].to_list()
    centroids_df['K_neighbors']   = ci_df['K_neighbors'].to_list()
    centroids_df['cluster_id']    = ci_df['cluster_id'].to_list()
    centroids_df['public']        = public_df['public'].to_list()
    centroids_df['n_subjects']    = public_df['n_subjects'].to_list()
    centroids_df['size_order']    = list(range(0,centroids_df.shape[0]))
    
    tr.centroids_df = centroids_df.copy()
    return tr


def get_basic_centroids(tr, max_dist = 200, look = False, cdr3_name = 'cdr3_b_aa'):

    # Cluster based on the max_dist
    ci = simple_cluster_index(tr.pw_beta, t = max_dist)

    # Get a DataFrame for clusters (It is returned largest to smallest)
    ci_df = cluster_index_to_df(cluster_index = ci)  

    # Determine degree of public sharing in each cluster
    publicities = list()
    for i,r in ci_df.iterrows():
        clone_cluster_df = tr.clone_df.iloc[r['neighbors'],]
        publicity,n_subjects = ispublic(clone_cluster_df)
        publicities.append((publicity,n_subjects))
    # Store public status in a DataFrame
    public_df = pd.DataFrame(publicities).rename(columns = {0:'public',1:'n_subjects'})

    # Iterate through the ci_df DataFrame
    counter = 0
    centroids = list()
    cluster_df_list = list()
    for i,r in ci_df.iterrows():
        counter = counter + 1
        
        clone_cluster_df = tr.clone_df.iloc[r['neighbors'],]
        cluster_df_list.append(clone_cluster_df )
        publicity = ispublic(clone_cluster_df)
        try:
            centroid, dmatrix, iloc_ind, loc_ind, = get_centroid_seq(df = tr.clone_df.iloc[r['neighbors'],])
            assert clone_cluster_df.iloc[iloc_ind,][ cdr3_name] == centroid
            assert tr.clone_df.iloc[loc_ind,][ cdr3_name] == centroid
            centroids.append( clone_cluster_df.iloc[iloc_ind,].reset_index(drop=True).copy())
            if look: 
                #print(dmatrix)
                print(r['neighbors'])
                print(f"MAX DISTANCE IN CLUSTER: {np.max(dmatrix)}" )
                print("----")
        except Exception as e:
            print("Using first sequence when not enough available to find centroid.\n")
            #print('FAIL')
            centroids.append( clone_cluster_df.iloc[0,].reset_index(drop=True).copy())


    # Save the Information in a centroids_df

    # renames --- {0: 'cell_type', 1: 'subject', 2: 'v_b_gene', 3: 'j_b_gene', 4: 'cdr3_b_aa', 5: 'cdr3_b_nucs, 6: 'cdr1_b_aa', 7: 'cdr2_b_aa', 8: 'pmhc_b_aa', 9: 'count', 10: 'clone_id'}
    renames = {i:v for i,v in enumerate(tr.clone_df.columns.to_list())}
    centroids_df = pd.DataFrame(centroids).rename(columns = renames)
    centroids_df['neighbors']     = ci_df['neighbors'].to_list()
    centroids_df['K_neighbors']   = ci_df['K_neighbors'].to_list()
    centroids_df['cluster_id']    = ci_df['cluster_id'].to_list()
    centroids_df['public']        = public_df['public'].to_list()
    centroids_df['n_subjects']    = public_df['n_subjects'].to_list()
    centroids_df['size_order']    = list(range(0,centroids_df.shape[0]))
    
    tr.centroids_df = centroids_df.copy()
    return tr