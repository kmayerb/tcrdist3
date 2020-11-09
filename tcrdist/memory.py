import numpy as np
import pandas as pd
import copy
import os
import parmap 
import secrets
from progress.bar import IncrementalBar
from hierdiff.tally import neighborhood_tally
import pwseqdist


from scipy import sparse
"""
Large datasets can consume more memory than is typically available. 

The memory module solves this problem with functions
that allow for TCRdistance computation to 
be broken down into fragments that can be stored on disk
and reassembled in a sparser format. 
"""

def _partition(l, n):
    """
    _partition breaks up a list in parts of size n

    Parameters
    ----------
    l : list
    
    n : int

    >>> _partition([1,2,3,4,5,6], 2)
    [[1, 2], [3, 4], [5, 6]]
    >>> _partition([1,2,3,4,5,6], 3)
    [[1, 2, 3], [4, 5, 6]]
        """
    n = int(n)
    r = [l[i:i + n] for i in range(0, len(l), n)]
    return r


def gen_sparse_rw_on_fragment(tcrrep, ind, outfile, matrix_name = 'rw_beta', max_distance=50):
    """
    gen_sparse_rw_on_fragment generates a sparse matrix of distances
    on a fragment of overall comparisons. 

    Suppose a clone_df has m unique clones. A full matrix is m x m. 
    this computes distances of the rectangular fragment 
    m[ind,] x m. The fragment is converted to a sparse compressed 
    format. All true zeros are set to 1, and all values above <max_distance>
    are set to one, for a massive space savings.

    tcrrep : TCRrep
        TCRrep instance 
    ind : list
        line of row index position specifying the rows of the clone_df
    outfile : str
        string name of the fragment
    matrix_name  : str
        the matrix attribute to be stored ('rw_beta')

    Example
    -------
    See usage in tcrdist.rep_funcs.compute_pw_sparse_out_of_memory

    """
    tr = copy.deepcopy(tcrrep)
    tr.compute_rect_distances(df = tr.clone_df.iloc[ind,], df2 = tr.clone_df, store = False)
    M = getattr(tr, matrix_name)
    # Set all true 0s to -1
    M[M == 0] = -1
    # Set all values greater than max_distance to zero to create massive scarcity
    M[M > max_distance] = 0 
    # convert to scr matrix spa
    M = sparse.csr_matrix(M)
    # Write to outfile
    sparse.save_npz(outfile, M)
    del tr
    del M
    return True

def gen_sparse_rw_on_fragment2(tcrrep, ind, outfile, max_distance=50):
    """
    gen_sparse_rw_on_fragment generates a sparse matrix of distances
    on a fragment of overall comparisons. 

    Suppose a clone_df has m unique clones. A full matrix is m x m. 
    this computes distances of the rectangular fragment 
    m[ind,] x m. The fragment is converted to a sparse compressed 
    format. All true zeros are set to 1, and all values above <max_distance>
    are set to one, for a massive space savings.

    tcrrep : TCRrep
        TCRrep instance 
    ind : list
        line of row index position specifying the rows of the clone_df
    outfile : str
        string name of the fragment
    matrix_name  : str
        the matrix attribute to be stored ('rw_beta')

    Example
    -------
    See usage in tcrdist.rep_funcs.compute_pw_sparse_out_of_memory

    """
    tr = copy.deepcopy(tcrrep)
    tr.compute_rect_distances(df = tr.clone_df.iloc[ind,], df2 = tr.clone_df, store = False)
    for chain in tr.chains:
        M = getattr(tr, f"rw_{chain}")
        # Set all true 0s to 1s
        M[M == 0] = 1
        # Set all values greater than max_distance to zero to create massive scarcity
        M[M > max_distance] = 0 
        # convert to scr matrix spa
        M = sparse.csr_matrix(M)
        # Write to outfile
        x = f"rw_{chain}"
        sparse.save_npz(f"{outfile}.{x}.npz", M)
        del M
    
    del tr
    return True


def collapse_csrs(list_of_csr_filenames, axis = 0):
    """
    Given a list of filenames referring to the sparse
    representation of ordered row-wise (axis = 0) chunks
    of a matrix, collapse chunks into a single 
    matrix in sparse format.

    Parameters
    ----------
    list_of_csr_filenames : list of strings
        list of string pointing to filenames of .npz sparse matrices
    axis : int
        if 0, combine with vstack, else if 1

    Example
    -------
    >>> import numpy as np
    >>> from scipy import sparse
    >>> A = np.array([[1, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 1], [0, 0, 0, 2, 0, 0]])
    >>> B = np.array([[3, 3, 3, 1, 0, 0], [5, 5, 2, 0, 0, 1], [6, 7, 0, 2, 0, 0]])
    >>> A_s = sparse.csr_matrix(A)
    >>> B_s = sparse.csr_matrix(B)
    >>> sparse.save_npz("A.npz", A_s)
    >>> sparse.save_npz("B.npz", B_s)
    >>> AB_s = collapse_csrs(["A.npz","B.npz"], axis = 0)
    >>> AB = AB_s.todense()
    >>> assert np.all(AB == np.concatenate([A,B], axis = 0))
    """
    assert axis in [0,1], "axis must be 0 (for hstack on rows) or 1 (for vstack on columns)."
    if axis == 0:
        full_spase_matrix = sparse.vstack([sparse.load_npz(f_npz) for f_npz in list_of_csr_filenames])
    if axis == 1:
        full_spase_matrix = sparse.hstack([sparse.load_npz(f_npz) for f_npz in list_of_csr_filenames])
    return full_spase_matrix



def gen_n_tally_on_fragment(tcrrep, 
                            ind, 
                            infile, 
                            outfile,
                            x_cols = ['epitope'], 
                            count_col='count', 
                            knn_neighbors= None, 
                            knn_radius =50):
    import copy
    tr = copy.deepcopy(tcrrep)
    #with open(infile ,"rb") as frag:
    #   rwmat = np.load(frag)
    rwmat = sparse.load_npz(infile)
    rwmat = np.asarray(rwmat.todense())
    rwmat[rwmat == 0] = 500
    ndif = neighborhood_tally(  df_pop = tr.clone_df, 
                                pwmat = rwmat,#tr.pw_beta[ind,], 
                                x_cols = x_cols, 
                                df_centroids=tr.clone_df.iloc[ind,],
                                count_col=count_col, 
                                knn_neighbors= knn_neighbors, 
                                knn_radius =knn_radius)

    ndif.to_csv(outfile, index = False)
    
    del ndif
    
    return outfile


def gen_n_tally_on_fragment2(tcrrep, 
                            ind, 
                            infile, 
                            outfile,
                            x_cols = ['epitope'], 
                            count_col='count', 
                            knn_neighbors= None, 
                            knn_radius =50):
    import copy
    tr = copy.deepcopy(tcrrep)
    
    rwmat_dict = dict()
    for chain in tr.chains:
        rwmat = sparse.load_npz(f"{infile}.rw_{chain}.npz")
        rwmat = np.asarray(rwmat.todense())
        rwmat[rwmat == 0] = 500
        rwmat_dict[chain] = rwmat

    if len(tr.chains) == 1:
        rwmat = rwmat_dict[tr.chains[0]]
    if len(tr.chains) == 2:
        rwmat = np.add(rwmat_dict[tr.chains[0]], rwmat_dict[tr.chains[1]])
        
    ndif = neighborhood_tally(  df_pop = tr.clone_df, 
                                pwmat = rwmat,#tr.pw_beta[ind,], 
                                x_cols = x_cols, 
                                df_centroids=tr.clone_df.iloc[ind,],
                                count_col=count_col, 
                                knn_neighbors= knn_neighbors, 
                                knn_radius =knn_radius)

    ndif.to_csv(outfile, index = False)
    
    del ndif
    
    return outfile


def _concat_to_file(dest, fragments, filename = "nndif.csv" ):
    """
    _concat_to_file take a list of .csv filenames, <fragments>. 
    and combines them into one concatenated file 
    written written to the <dest> directory on disk.

    Parameters
    ----------
    dest : str
        location where file is to be written 
    fragments : list
        list of .csv files
    Returns
    -------
    concatenated_filename : str
        filename of the concatenated file
    """
    from progress.bar import IncrementalBar
    bar = IncrementalBar(f'Processing Files', max = len(fragments), suffix='%(percent)d%%')
    mysep = ","
    concatenated_filename = os.path.join(dest, filename)
    counter = 0
    for f in fragments:
        df = pd.read_csv(os.path.join(f), sep = ',')
        if counter == 0:
            df.to_csv(concatenated_filename, header=True, index = False, sep = mysep)
        if counter > 0:
            df.to_csv(concatenated_filename, header=False, index = False, sep = mysep, mode = "a")
        del df
        counter =+1 
        bar.next()
    bar.next(); bar.finish()
    print(f"WROTE TO {concatenated_filename}")
    return concatenated_filename
    

def _concat_to_memory(fragments):
    """
    This
    Parameters
    ----------
    fragments : list
        list of .csv files
    Returns
    -------
    concatenated_df : DataFrame
        
    """
    from progress.bar import IncrementalBar
    bar = IncrementalBar(f'Processing Files', max = len(fragments), suffix='%(percent)d%%')
    mysep = ","
    dfs = list()
    for f in fragments:
        dfs.append(pd.read_csv(os.path.join(f), sep = ','))
        bar.next()
    concatenated_df = pd.concat(dfs)
    bar.finish()
    return concatenated_df


def _sparse_cdr3_tcrdist_shard(chunk1_indices, components, df, metrics, weights, kargs, radius, df2):
    """Returns a sparse matrix with 0/NA meaning D > radius and -1 meaning D = 0"""
    
    """For now not trying to re-expand less as I expand components and find D > radius"""
    # keep_ind = np.ones(chunk1_indices.shape[0], dtype=bool)

    if df2 is None:
        n2 = df.shape[0]
    else:
        n2 = df2.shape[0]

    metric_keys = [k for k in metrics.keys() if not 'cdr3' in k]
    
    metric_cdr3_keys = [k for k in metrics.keys() if 'cdr3' in k]
    weight_cdr3_keys = [k for k in weights.keys() if 'cdr3' in k]
    assert metric_cdr3_keys == weight_cdr3_keys, "metrics and weights cdr3 keys must be identical"
    
    if kargs is not None:
        kargs_cdr3_keys  = [k for k in kargs.keys() if 'cdr3' in k]
        assert metric_cdr3_keys == kargs_cdr3_keys,  "metrics and kargs cdr3 keys must be identical"

    tcrdist = np.zeros((chunk1_indices.shape[0], n2), dtype=np.int16)
    for k in metric_keys:
        """Re-expand chunk of rows and all of columns: weights already applied"""
        tcrdist = tcrdist + components[k][0][components[k][1][chunk1_indices], :][:, components[k][2]]

    """Make a list of pairs to compute CDR3 distance for this chunk.
    Then create sparse CSR matrix for the chunk and return"""
    for k in metric_cdr3_keys:
        """Compute CDR3s sparsely"""
        """Only compute CDR3 distance for if the running total <= radius:
        Not currently possible with running rect and sparse doesn't avoid redundant distances. Opting
        for computing the whole chunk with pwrect"""

        """Since this is just a chunk i could use regular pwrect and not worry about running.
        Running won't be very helpful since it won't incorporate the radius from the other components"""
        # apply_running_rect(metric, seqs1, seqs2, radius=radius, density_est=density_est, *args, ncpus=1, uniqify=True, alphabet=parasail_aa_alphabet, **kwargs)
        if df2 is None:
            seqs2 = df[k].values
        else:
            seqs2 = df2[k].values
        pwrect = pwseqdist.apply_pairwise_rect(metric=metrics[k], 
                                        seqs1=df[k].values[chunk1_indices], 
                                        seqs2=seqs2, 
                                        ncpus=1, 
                                        uniqify=True, 
                                        reexpand=True,
                                        **kargs[k])
        tcrdist = tcrdist + weights[k] * pwrect
    zind = np.nonzero(tcrdist == 0)
    tcrdist[zind] = -1
    keep_ind = np.nonzero(tcrdist <= radius)
    S = sparse.csr_matrix((tcrdist[keep_ind], keep_ind), shape=tcrdist.shape)
    return S