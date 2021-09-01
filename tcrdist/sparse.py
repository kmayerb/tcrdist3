from scipy import sparse
import numpy as np

def add_sparse_pwd(a, b, neg1_fmt=True):
    """
    Add sparse pairwise distance matrices a and b.
    Can work with two "formats": true zero distances as -1
    or true zero distances as 0. Returns same format as received.
    Does not attempt to add sparse zeros (ie NA),
    so returns an intersection of data in a and b (summed).
    Returns a scipy.sparse.csr matrix (can take csr or any
    format as inputs).
    Parameters
    ----------
    a, b : scipy.sparse.csr or coo
        Sparse pairwise distance matrices. D=0 can be 
        represented as -1 or 0 but should be specified
        using neg1_fmt
    neg1_fmt : bool
        Indicator whether to assume D=0 is represented as -1
    Returns
    -------
    tot : scipy.sparse.csr
        Sparse matrix representing sum of a and b
        where a and b both have values.
    """
    
    a_coo = a.tocoo()
    b_coo = b.tocoo()
    nr, nc = a_coo.shape
    
    ind_a = a_coo.row.astype(np.int64) * nc + a_coo.col
    ind_b = b_coo.row.astype(np.int64) * nc + b_coo.col
    # ind_ab = np.array(tuple(set(ind_a).intersection(set(ind_b))))
    ind_ab = np.intersect1d(ind_a, ind_b, assume_unique=True)
    keep_a = np.in1d(ind_a, ind_ab, assume_unique=True)
    keep_b = np.in1d(ind_b, ind_ab, assume_unique=True)
    if neg1_fmt:
        a_data = a_coo.data[keep_a]
        a_data[a_data == -1] = 0
        b_data = b_coo.data[keep_b]
        b_data[b_data == -1] = 0
        tot = a_data + b_data
    else:
        tot = a_coo.data[keep_a] + b_coo.data[keep_b]
    """Need to add -1 as placeholder upon creation of new sparse matrix"""
    tot[tot == 0] = -1
    aplusb = sparse.coo_matrix((tot,
                                (a_coo.row[keep_a], a_coo.col[keep_a])), shape=(nr, nc))
    aplusb = aplusb.tocsr()
    if not neg1_fmt:
        """Convert it back to 0 = true zero distance"""
        aplusb.data[aplusb.data == -1] = 0
    return aplusb
