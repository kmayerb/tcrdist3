"""
tcrdist3: Functional Programming Approach For Higher Memory Demand Use Cases
"""

import pwseqdist as pw
import pandas as pd
import numpy as np
from scipy import sparse

import os
from tcrdist import memory
import secrets
import parmap
import shutil
import multiprocessing

__all__ = ['_pw',
           '_pws',
           'compute_pw_sparse_out_of_memory',
           'compute_pws_sparse',
           'pw2dense']

def pw2dense(pw, maxd):
    """Make a pairwise distance matrix dense
    assuming -1 is used to encode D = 0"""
    pw = np.asarray(pw.todense())
    pw[pw == 0] = maxd + 1
    # pw[np.diag_indices_from(pw)] = 0
    pw[pw == -1] = 0
    return pw


def _pws(df, metrics, weights, kargs, df2 = None, cpu = 1, uniquify = True, store = False):
    """    
    _pws performs pairwise distance calculation across a multiple 
    columns of a Pandas DataFrame. This naturally permits calculation of
    a CDR-weighted tcrdistance that incorporates dissimilarity across 
    multiple complimentarity determining regions (see example below):

    Parameters 
    ----------
    df : pd.DataFrame
        Clones DataFrame containing, at a minimum, columns with CDR sequences
    df2 : pd.DataFrame or None
        Second clones DataFrame containing, at a minimum, columns with CDR sequences
    metrics : dict
        Dictionary of functions, specifying the distance metrics to apply to each CDR
    weights : dict
        Weights determining the contributions of each CDR distance to the aggregate distance
    kargs : dict
        Dictionary of Dictionaries
    cpu : int 
        Number of available cpus
    use_numba : bool
        If True, use must use a numba compatible metric 
    store : bool
        If False, only full tcrdist is returned. If True, 
        all component distance matrices are returned.
        
    Returns
    -------
    s : dictionary with tcr_distance.

    Example
    -------
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.rep_funcs import _pw, _pw2
    
    # Define metrics for each region
    metrics = { "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist}

    # Define weights
    weights = { 
                "cdr3_a_aa" : 3,
                "pmhc_a_aa" : 1,
                "cdr2_a_aa" : 1,
                "cdr1_a_aa" : 1,
                "cdr3_b_aa" : 3,
                "pmhc_b_aa" : 1,
                "cdr2_b_aa" : 1,
                "cdr1_b_aa" : 1}

    kargs = {   "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}

    df = pd.DataFrame("dash2.csv")
    _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    """
    metric_keys = list(metrics.keys())
    weight_keys = list(weights.keys())
    assert metric_keys == weight_keys, "metrics and weights keys must be identical"
    
    if kargs is not None:
        kargs_keys  = list(kargs.keys())
        assert metric_keys == kargs_keys,  "metrics and kargs keys must be identical"
    
    tcrdist = None
    s = dict()
    for k in metric_keys:
        if df2 is None:
            pw_mat = _pw(seqs1 = df[k].values, metric = metrics[k], ncpus = cpu, uniqify= uniquify, **kargs[k])
        else:
            pw_mat = _pw(seqs1 = df[k].values, seqs2 = df2[k].values, metric = metrics[k], ncpus = cpu, uniqify= uniquify, **kargs[k])
            
        if store:
           s[k] = pw_mat 
        if tcrdist is None:
            tcrdist = np.zeros(pw_mat.shape, dtype=np.int16)
        tcrdist = tcrdist + (weights[k] * pw_mat)
    
    s['tcrdist'] = tcrdist
    return s


def _pw(metric, seqs1, seqs2=None, ncpus=1, uniqify= True, use_numba = False, **kwargs):
    """
    This is a wrapper for accessing pwseqdist version > 0.2.
    No matter what, it returns squareform results
    """
    pw_mat = pw.apply_pairwise_rect(metric = metric, 
                                    seqs1  = seqs1, 
                                    seqs2  = seqs2, 
                                    ncpus  = ncpus, 
                                    uniqify = uniqify, 
                                    use_numba = use_numba,
                                    **kwargs)
    # if len(pw_mat.shape) == 1:
    #     from scipy.spatial.distance import squareform
    #     pw_mat = squareform(pw_mat)

    return pw_mat


def compute_pws_sparse(df, metrics, weights, kargs, radius=50, df2=None, cpu=1, chunk_size=500, store=False, pm_pbar=True):
    """   
    compute_pw_sparse performs pairwise distance calculation across multiple 
    columns of a Pandas DataFrame. This naturally permits calculation of
    a CDR-weighted tcrdistance that incorporates dissimilarity across 
    multiple complimentarity determining regions (see example below):

    Radius is applied per chain.
    TODO: allow for different radius per chain

    Parameters 
    ----------
    df : pd.DataFrame
        Clones DataFrame containing, at a minimum, columns with CDR sequences
    df2 : pd.DataFrame or None
        Second clones DataFrame containing, at a minimum, columns with CDR sequences
    metrics : dict
        Dictionary of functions, specifying the distance metrics to apply to each CDR
    weights : dict
        Weights determining the contributions of each CDR distance to the aggregate distance
    kargs : dict
        Dictionary of Dictionaries
    cpu : int 
        Number of available cpus
    radius : int
        Only distances <= radius will be retained.
    chunk_size : int
        Number of rows that will make up each row chunk of [df x df2]
        which will be computed densely.
    store : bool
        IGNORED
        
    Returns
    -------
    s : dictionary with tcr_distance.

    Example
    -------
    import pwseqdist as pw
    import pandas as pd
    from tcrdist.rep_funcs import _pw, _pw2
    
    # Define metrics for each region
    metrics = { "cdr3_a_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_a_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr3_b_aa" : pw.metrics.nb_vector_tcrdist,
                "pmhc_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr2_b_aa" : pw.metrics.nb_vector_tcrdist,
                "cdr1_b_aa" : pw.metrics.nb_vector_tcrdist}

    # Define weights
    weights = { 
                "cdr3_a_aa" : 3,
                "pmhc_a_aa" : 1,
                "cdr2_a_aa" : 1,
                "cdr1_a_aa" : 1,
                "cdr3_b_aa" : 3,
                "pmhc_b_aa" : 1,
                "cdr2_b_aa" : 1,
                "cdr1_b_aa" : 1}

    kargs = {   "cdr3_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_a_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr3_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':3, 'ctrim':2, 'fixed_gappos':False},
                "pmhc_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr2_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True},
                "cdr1_b_aa" : {'use_numba': True, 'distance_matrix': pw.matrices.tcr_nb_distance_matrix, 'dist_weight': 1, 'gap_penalty':4, 'ntrim':0, 'ctrim':0, 'fixed_gappos':True}}

    df = pd.DataFrame("dash2.csv")
    _pws(df = df, metrics = metrics, weights= weights, kargs=kargs, cpu = 1, store = False)
    """
    metric_keys = [k for k in metrics.keys() if not 'cdr3' in k]
    weight_keys = [k for k in weights.keys() if not 'cdr3' in k]
    assert metric_keys == weight_keys, "metrics and weights keys must be identical"
    
    if kargs is not None:
        kargs_keys  = [k for k in kargs.keys() if not 'cdr3' in k]
        assert metric_keys == kargs_keys,  "metrics and kargs keys must be identical"
    
    n1 = df.shape[0]
    
    """Compute all but CDR3 as normal, but do not reexpand.
    Computing unique distances should not be memory or CPU intensive"""
    tcrdist = None
    components = dict()
    for k in metric_keys:
        if df2 is None:
            seqs2 = None
        else:
            seqs2 = df2[k].values

        """With reexapnd = False, returns: pw_mat, uind_i1, uind_i2"""
        pwmat, ind1, ind2 = pw.apply_pairwise_rect(metric=metrics[k], 
                                                    seqs1=df[k].values, 
                                                    seqs2=seqs2, 
                                                    ncpus=min(cpu, 2),
                                                    uniqify=True, 
                                                    reexpand=False,
                                                    **kargs[k])

        components[k] = (pwmat * weights[k], ind1, ind2) 

    """Can't do this because it will be huge. Also can't compute list of
    potential D < radius because that also could be too large.
    But need to pre-compute these non-CDR3 because otherwise the CDR3
    matrix won't be sparse enough with the radius.

    Solution: chunk the computation here. Only compute subsets of potential
    seqs pairs sparsely and spread across processors."""

    if cpu > 1 and n1 > chunk_size:
        """Chunk along df (rows) only"""
        chunk_func = lambda l, n: [np.array(l[i:i + n], dtype=np.int64) for i in range(0, len(l), n)]
        # chunksz = max(len(pw_indices) // cpu, 1)
        """Chunked indices is a list of arrays of indices"""
        """List of the chunked [chunk_size,] arrays"""
        chunked_indices = chunk_func(np.arange(n1, dtype=np.int64), chunk_size)
        
        with multiprocessing.Pool(cpu) as pool:
            dists = parmap.map(memory._sparse_cdr3_tcrdist_shard,
                               chunked_indices,
                               components,
                               df,
                               metrics,
                               weights,
                               kargs,
                               radius,
                               df2,
                               pm_parallel=True,
                               pm_pool=pool,
                               pm_pbar=pm_pbar)

        full_S = sparse.vstack(dists)
    else:
        full_S = memory._sparse_cdr3_tcrdist_shard(np.arange(n1, dtype=np.int64),
                                                components,
                                                df,
                                                metrics,
                                                weights,
                                                kargs,
                                                radius,
                                                df2)  
    
    return {'tcrdist': full_S}



def compute_n_tally_out_of_memory(fragments, 
                                 matrix_name = "rw_beta", 
                                 to_file = False, 
                                 to_memory = True,
                                 pm_processes = 2,
                                 **kwargs):
    """
    Parameters 
    ---------
    fragments : tuple 
        3-part tuple
            0 : TCRrep instance
            1 : list of rows in each fragment
            2 : filename holding the distances as a .npz
    matrix_name : str
        For 'beta' chain dists use 'rw_beta'
    to_file : bool
        If True, than a file is written to disk containing all the n_tally
    to_memory : bool
        If True, the n_tally is loaded directly to memory 
    **kwargs
        Keyword arguments are passed directly to tcridist.memory.gen_n_tally_on_fragment,
        including:  
            x_cols : list
                categorical variable 'epitope']  
            count_col : str
                column for counts  
            knn_neighbors : bool
                option to choose a fixed number of neighbors (NOT RECOMMENDED) 
            knn_radius : int
                maximum radius for finding neighbors

    Notes
    -----
    We envision that user may want to run nndif 
    multiple times with different categorical variables
    and may not want to go through the computational 
    intensive steps of computing TCRdistance each time. 

    This can be accomodated via setting cleanup to False:

    compute_pw_sparse_out_of_memory(cleanup=False) 

    Example 
    -------
    from scipy import sparse
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory
    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df,               
                organism = 'mouse',
                chains = ['beta'],
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = True,
                store_all_cdr = False)

    S, chunks = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 1000, cleanup = False)
    # dest contains all the shards
    """

    # [(<tcrdist.repertoire.TCRrep at 0x14035b6d0>,
    #   range(0, 500),
    #   'd3be945e8956/0.rw_beta.npz'),

    # rearrange fragments in order (tr, ind, .npz, .csv)
    fragments =  [(x[0], x[1], x[2], f"{x[2]}.nndif.csv") for x in fragments ] 
    #fragments = [(tr, f"{dest}/{i}.{matrix_name}.npz", ind, f"{dest}/{i}.nndif.cvs") for i,ind in enumerate(row_chunks)] 
    csvfragments = parmap.starmap(memory.gen_n_tally_on_fragment, fragments, **kwargs, pm_pbar=True, pm_processes = pm_processes)
    if to_file:
        dest =os.path.dirname(csvfragments[0])
        nndiff_file = memory._concat_to_file(dest =dest, fragments =csvfragments)
        return nndiff_file

    if to_memory:
        nndiff = memory._concat_to_memory(fragments = csvfragments)
        return nndiff



def compute_pw_sparse_out_of_memory(tr,
                                    row_size      = 500,
                                    pm_processes  = 2,
                                    pm_pbar       = True,
                                    max_distance  = 50,
                                    matrix_name   = 'rw_beta',
                                    reassemble    = True,
                                    cleanup       = True):
    """
    Instead of calling TCRrep.compute_distances(), this 
    function permits a parallelizable approach that does 
    not require holding a large matrix in memory. 

    Default behavior is to reassemble a scipy
    sparse matrix from a set of sub matrices written to disk fragment. 
    With <reassemble = True>  function returns a scipy sparse matrix. 
    Space savings are achieved because any value above <max_distance> is set to zero. 
    True zero distances are set to 1. 

    Can be used to form a network of TCRs with tcrdistances < max_distance,

    Parameters
    ----------
    tr : TCRrep
        TCRrep instance with clone_df
    row_size : int
        How many rows to process in memory at once
    pm_processes : int
        Numbe of concurrent parallel processes to run at once
    pm_bar : bool 
        If True, show progress bar.
    max_distance : int
        Max distance
    matrix_name : str
        Name of matrix to return (i.e, 'rw_beta' or 'rw_alpha')
    reassemble: True
        If true, makes one matrix from all the sparse sub matrices. 
    cleanup: bool,
        if True, deletes temporary files. 

    Returns
    -------
    csr_full : sparse scipy matrix
            
    dest : str
        name of the folder that holds fragments

    
    Examples
    --------
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory

    df = pd.read_csv("dash.csv")
                          #(1)
    tr = TCRrep(cell_df = df,               #(2)
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = True,
                store_all_cdr = False)

    S = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 1000)
    # S is a <1920x1920 sparse matrix of type '<class 'numpy.int16'>'
    M = S.todense()
    M[M==1] = 0 
    np.all(M == tr.pw_beta)
    S, chunks = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 50)
    print(S)
    # S is a <1920x1920 sparse matrix of type '<class 'numpy.int16'>'
    """ 
    dest = secrets.token_hex(6)
    os.mkdir(dest)
    print(f"CREATED /{dest}/ FOR HOLDING DISTANCE OUT OF MEMORY")
    row_chunks = memory._partition(range(tr.clone_df.shape[0]), row_size)

    smatrix_chunks = [(tr, ind,  f"{dest}/{i}.{matrix_name}.npz") for i,ind in enumerate(row_chunks)]
    csrfragments = parmap.starmap(memory.gen_sparse_rw_on_fragment, 
            smatrix_chunks, 
            matrix_name = matrix_name, 
            max_distance=max_distance, 
            pm_pbar=pm_pbar, 
            pm_processes = pm_processes)
    if reassemble:
        csr_full = memory.collapse_csrs([x[2] for x in smatrix_chunks])
        print(f"RETURNING scipy.sparse csr_matrix w/dims {csr_full.shape}")
    else: 
        csr_full = None
        
    if cleanup: 
        assert os.path.isdir(dest)
        print(f"CLEANING UP {dest}")
        shutil.rmtree(dest)
    
    
    return csr_full, smatrix_chunks



def compute_pw_sparse_out_of_memory2(tr,
                                    row_size      = 500,
                                    pm_processes  = 2,
                                    pm_pbar       = True,
                                    max_distance  = 50,
                                    reassemble    = True,
                                    cleanup       = True,
                                    assign        = True):
    """
    Instead of calling TCRrep.compute_distances(), this 
    function permits a parallelizable approach that does 
    not require holding a large matrix in memory. 

    Default behavior is to reassemble a scipy
    sparse matrix from a set of sub matrices written to disk fragment. 
    With <reassemble = True>  function returns a scipy sparse matrix. 
    Space savings are achieved because any value above <max_distance> is set to zero. 
    True zero distances are set to -1. 

    Can be used to form a network of TCRs with tcrdistances < max_distance,

    Parameters
    ----------
    tr : TCRrep
        TCRrep instance with clone_df
    row_size : int
        How many rows to process in memory at once
    pm_processes : int
        Numbe of concurrent parallel processes to run at once
    pm_bar : bool 
        If True, show progress bar.
    max_distance : int
        Max distance
    matrix_name : str
        Name of matrix to return (i.e, 'rw_beta' or 'rw_alpha')
    reassemble: True
        If true, makes one matrix from all the sparse sub matrices. 
    cleanup: bool,
        if True, deletes temporary files. 
    assign : bool 
        if True, assigns pw sparse matrices to TCRrep object.
        That is TCRrep.pw_beta, TCRrep.pw_alpha will be assigned 
        the reassembled spare matrces.

    Returns
    -------
    csr_full : sparse scipy matrix
            
    dest : str
        name of the folder that holds fragments

    
    Examples
    --------
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory

    df = pd.read_csv("dash.csv")
                          #(1)
    tr = TCRrep(cell_df = df,               #(2)
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = True,
                store_all_cdr = False)

    S = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 1000)
    # S is a <1920x1920 sparse matrix of type '<class 'numpy.int16'>'
    M = S.todense()
    M[M==1] = 0 
    np.all(M == tr.pw_beta)
    S, chunks = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 50)
    print(S)
    # S is a <1920x1920 sparse matrix of type '<class 'numpy.int16'>'
    """ 

    # Early warning to save heartache
    if assign is True and reassemble is False:
        raise ValueError("If you want to assign results to a TCRrep instance, you must set reassemble to True")

    dest = secrets.token_hex(6)
    os.mkdir(dest)
    print(f"CREATED /{dest}/ FOR HOLDING DISTANCE OUT OF MEMORY")
    row_chunks = memory._partition(range(tr.clone_df.shape[0]), row_size)


    smatrix_chunks = [(tr, ind,  f"{dest}/{i}") for i,ind in enumerate(row_chunks)]
    csrfragments = parmap.starmap(memory.gen_sparse_rw_on_fragment2, 
            smatrix_chunks,  
            max_distance=max_distance, 
            pm_pbar=pm_pbar, 
            pm_processes = pm_processes)

    if reassemble:
        csr_full_dict = dict()
        for chain in tr.chains:
            chain_str = f"rw_{chain}"
            csr_full = memory.collapse_csrs([f"{x[2]}.{chain_str}.npz" for x in smatrix_chunks])
            print(f"RETURNING scipy.sparse csr_matrix w/dims {csr_full.shape}")
            csr_full_dict[chain] = csr_full
    else: 
        csr_full_dict= None
       
    if assign:
        for chain in tr.chains:
            setattr(tr, f"pw_{chain}", csr_full_dict[chain])

    if cleanup: 
        assert os.path.isdir(dest)
        print(f"CLEANING UP {dest}")
        shutil.rmtree(dest)
    
    return csr_full_dict, smatrix_chunks




def compute_n_tally_out_of_memory2(fragments, 
                                 to_file = False, 
                                 to_memory = True,
                                 pm_processes = 2,
                                 **kwargs):
    """
    Parameters 
    ---------
    fragments : tuple 
        3-part tuple
            0 : TCRrep instance
            1 : list of rows in each fragment
            2 : filename holding the distances as a .npz
    matrix_name : str
        For 'beta' chain dists use 'rw_beta'
    to_file : bool
        If True, than a file is written to disk containing all the n_tally
    to_memory : bool
        If True, the n_tally is loaded directly to memory 
    **kwargs
        Keyword arguments are passed directly to tcridist.memory.gen_n_tally_on_fragment,
        including:  
            x_cols : list
                categorical variable 'epitope']  
            count_col : str
                column for counts  
            knn_neighbors : bool
                option to choose a fixed number of neighbors (NOT RECOMMENDED) 
            knn_radius : int
                maximum radius for finding neighbors

    Notes
    -----
    We envision that user may want to run nndif 
    multiple times with different categorical variables
    and may not want to go through the computational 
    intensive steps of computing TCRdistance each time. 

    This can be accomodated via setting cleanup to False:

    compute_pw_sparse_out_of_memory(cleanup=False) 

    Example 
    -------
    from scipy import sparse
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory
    df = pd.read_csv("dash.csv")
    tr = TCRrep(cell_df = df,               
                organism = 'mouse',
                chains = ['beta'],
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = True,
                store_all_cdr = False)

    S, chunks = compute_pw_sparse_out_of_memory(tr, matrix_name = "rw_beta", max_distance = 1000, cleanup = False)
    # dest contains all the shards
    """

    # [(<tcrdist.repertoire.TCRrep at 0x14035b6d0>,
    #   range(0, 500),
    #   'd3be945e8956/0.rw_beta.npz'),

    # rearrange fragments in order (tr, ind, .npz, .csv)
    fragments =  [(x[0], x[1], x[2], f"{x[2]}.nndif.csv") for x in fragments ] 
    #fragments = [(tr, f"{dest}/{i}.{matrix_name}.npz", ind, f"{dest}/{i}.nndif.cvs") for i,ind in enumerate(row_chunks)] 
    csvfragments = parmap.starmap(memory.gen_n_tally_on_fragment2, fragments, **kwargs, pm_pbar=True, pm_processes = pm_processes)
    if to_file:
        dest =os.path.dirname(csvfragments[0])
        nndiff_file = memory._concat_to_file(dest =dest, fragments =csvfragments)
        return nndiff_file

    if to_memory:
        nndiff = memory._concat_to_memory(fragments = csvfragments)
        return nndiff



