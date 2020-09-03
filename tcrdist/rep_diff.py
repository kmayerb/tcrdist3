import pandas as pd
import numpy as np

import hierdiff as hd

__all__ = ['neighborhood_diff',
           'hcluster_diff',
           'member_summ']

def neighborhood_diff(clone_df, pwmat, x_cols, count_col='count', knn_neighbors=None, knn_radius=None, subset_ind=None, cluster_ind=None, test_method='fishers'):
    """Tests for association of categorical variables in x_cols with the neighborhood
    around each TCR in clone_df. The neighborhood is defined by the K closest neighbors
    using pairwise distances in pwmat, or defined by a distance radius, knn_radius.

    Uses hierdiff package (available on PyPI) for tallying counts in each cluster
    and running tests.

    The statistical tests made available by this function are limited and meant only
    as a way to scan for signals. More sophisticated testing/modeling frameworks
    should be considered for real-world problems.

    Use test_method = None to return a table of counts for all neighborhoods that can be saved as
    a CSV and used to run other, more sophisticated tests (e.g. edgeR or other regressions).

    Use Fisher's exact test (test_method='fishers') to detect enrichment/association of the neighborhood
    with one binary variable. For example, test the 2 x 2 table for each clone:

    +----+----+-------+--------+
    |         |  Neighborhood  |
    |         +-------+--------+
    |         | MEM+  |   MEM- |
    +----+----+-------+--------+
    |VAR |  0 | a     |    b   |
    |    +----+-------+--------+
    |    |  1 | c     |    d   |
    +----+----+-------+--------+

    Use the chi-squared test (test='chi2') to detect association across multiple variables.
    Note that with sparse neighborhoods Chi-squared tests are unreliable.

    Use the Cochran-Mantel-Haenszel test (test='chm') to test stratified 2 x 2 tables:
    one VAR vs. neighborhood, over several strata defined in other variables.
    Use x_cols[0] as the primary (binary) variable and other x_cols for the categorical
    strata-defining variables. This tests the overall null that OR = 1 for x_cols[0].
    A test is also performed for homogeneity of the ORs among the strata (Breslow-Day test).

    Params
    ------
    clone_df : pd.DataFrame [nclones x metadata]
        Contains metadata for each clone.
    pwmat : np.ndarray [nclones x nclones]
        Square distance matrix for defining neighborhoods
    x_cols : list
        List of columns to be tested for association with the neighborhood. 
        (If test_method ='fishers' than x_cols must be length 1 and contain a binary variable)
    count_col : str
        Column in clone_df that specifies counts.
        Default none assumes count of 1 cell for each row.
    knn_neighbors : int
        Number of neighbors to include in the neighborhood.
    knn_radius : float
        Radius for inclusion of neighbors within the neighborhood.
        Specify K or R but not both.
    cluster_ind : None or np.ndarray
        Indices into df specifying the neighborhoods for testing.
    test_method : str or None
        Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
        Cochran-Mantel-Haenszel test ("chm") for testing.

    Returns
    -------
    res : pd.DataFrame [nclones x results]
        Results from testing the neighborhood around each clone."""
    res = hd.neighborhood_tally(df_pop=clone_df,
                                  pwmat=pwmat,
                                  x_cols=x_cols,
                                  count_col=count_col,
                                  knn_neighbors=knn_neighbors,
                                  knn_radius=knn_radius)
    if not test_method is None:
        res = hd.cluster_association_test(res, y_col='cmember', method=test_method)
    return res

def hcluster_diff(clone_df, pwmat, x_cols, Z=None, count_col='count', subset_ind=None, hclust_method='complete', optimal_ordering=True, test_method='fishers'):
    """Tests for association of categorical variables in x_cols with each cluster/node
    in a hierarchical clustering of TCR clones with distances in pwmat.

    Uses hierdiff package (available on PyPI) for tallying counts in each cluster
    and running tests.

    The statistical tests made available by this function are limited and meant only
    as a way to scan for signals. More sophisticated testing/modeling frameworks
    should be considered for real-world problems.

    Use test_method = None to return a table of counts for all neighborhoods that can be saved as
    a CSV and used to run other, more sophisticated tests (e.g. edgeR or other regressions).

    Use Fisher's exact test (test='fishers') to detect enrichment/association of the neighborhood/cluster
    with one binary variable.

    Tests the 2 x 2 table for each clone:

    +----+----+-------+--------+
    |         |    Cluster     |
    |         +-------+--------+
    |         |  MEM+ |  MEM-  |
    +----+----+-------+--------+
    |VAR |  0 |   a   |    b   |
    |    +----+-------+--------+
    |    |  1 |   c   |    d   |
    +----+----+-------+--------+

    Use the chi-squared test (test='chi2') to detect association across multiple categorical variables.
    Note that with small clusters Chi-squared tests are unreliable.

    Use the Cochran-Mantel-Haenszel test (test='chm') to test stratified 2 x 2 tables:
    one VAR vs. neighborhood, over several strata defined in other variables.
    Use x_cols[0] as the primary (binary) variable and other x_cols for the categorical
    strata-defining variables. This tests the overall null that OR = 1 for x_cols[0].
    A test is also performed for homogeneity of the ORs among the strata (Breslow-Day test).

    Params
    ------
    clone_df : pd.DataFrame [nclones x metadata]
        Contains metadata for each clone.
    pwmat : np.ndarray [nclones x nclones]
        Square distance matrix for defining neighborhoods
    x_cols : list
        List of columns to be tested for association with the neighborhood
    count_col : str
        Column in clone_df that specifies counts.
        Default none assumes count of 1 cell for each row.
    subset_ind : None or np.ndarray with partial index of df, optional
        Provides option to tally counts only within a subset of df, but to maintain the clustering
        of all individuals. Allows for one clustering of pooled TCRs,
        but tallying/testing within a subset (e.g. participants or conditions)
    hclust_method : str
        Method for hierarchical clustering, passed to the scipy.clustering.hierarchy
        linkage function.
    optimal_ordering : bool
        Flag passed to the scipy.clustering.hierarchy linkage function to improve
        visual tree layout. Can be slow for large trees.
    test_method : str or None
        Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
        Cochran-Mantel-Haenszel test ("chm") for testing.

    Returns
    -------
    res : pd.DataFrame [nclusters x results]
        Results from testing each cluster.
    Z : linkage matrix [clusters, 4]
        Clustering result returned from scipy.cluster.hierarchy.linkage"""
    res, Z = hd.hcluster_tally(df=clone_df,
                                  pwmat=pwmat,
                                  x_cols=x_cols,
                                  Z=Z,
                                  count_col=count_col,
                                  subset_ind=subset_ind,
                                  method='complete',
                                  optimal_ordering=optimal_ordering)
    if not test_method is None:
        res = hd.cluster_association_test(res, y_col='cmember', method=test_method)
    return res, Z


def member_summ(res_df, clone_df, key_col = 'neighbors_i', count_col='count', addl_cols=[], addl_n=1):
    """Return additional summary info about each result (row)) based on the members of the cluster.

    This is helpful for preparing strings to add to the tooltip in hierdiff.plot_hclust_props.
    
    Parameters
    ----------
    res_df : pd.DataFrame [nclusters x result cols]
        Returned from neighborhood_diff or hcluster_diff
    clone_df : pd.DataFrame [nclones x metadata]
        Contains metadata for each clone.
    key_col : str
        Column in res_df that specifies the iloc of members in the clone_df
    count_col : str
        Column in clone_df that specifies counts.
        Default none assumes count of 1 cell for each row.
    addl_cols : list
        Columns to summarize
    addl_n : int
        Number of top N clones to include in the summary of
        each cluster.

    Returns
    -------
    summ : pd.DataFrame [nclusters x summary columns]
        Columns that can be joined with res_df


    Example
    -------
    summ_df = member_summ(res_df, clone_df)
    res_df = res_df.join(summ_df, how='left')"""
    def _top_N_str(m, col, count_col, N):
        gby = m.groupby(col)[count_col].agg(np.sum)
        gby = 100 * gby / gby.sum()
        gby = gby.sort_values(ascending=False)
        out = ', '.join(['%s (%2.1f%%)' % (idx, v) for idx,v in gby.iteritems()][:N])
        return out
    
    split = []
    for resi, res_row in res_df.iterrows():
        m = clone_df.iloc[res_row[key_col]]

        mode_i = m[count_col].idxmax()
        summ = {}
        for c in [c for c in clone_df.columns if 'cdr3' in c]:
            summ[c] = _top_N_str(m, c, count_col, 1)
        for c in [c for c in clone_df.columns if 'gene' in c]:
            summ[c] = _top_N_str(m, c, count_col, 3)

        x_val_cols = [c for c in res_df.columns if 'x_val_' in c]
        x_freq_cols = [c for c in res_df.columns if 'x_freq_' in c]
        
        for label_col, freq_col in zip(x_val_cols, x_freq_cols):
            summ[res_row[label_col]] = np.round(res_row[freq_col], 3)

        for c in [c for c in addl_cols]:
            summ[c] = _top_N_str(m, c, count_col, addl_n)
        summ = pd.Series(summ, name=resi)
        split.append(summ)
    summ = pd.DataFrame(split)
    return summ
