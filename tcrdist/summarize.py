"""
Simple set of functions for summarizing over a group
"""

import pandas as pd
import numpy as np
import re


def filter_is(df, col, val):
    return df[df[col] == val]

def filter_in(df, col, vals):
    ind = [np.any(x in vals) for x in df[col]]
    return df[ind]

def filter_gt(df, col, val):
    return df[df[col] > val]

def filter_lt(df, col, val):
    return df[df[col] < val]


def is_subset(set1, set2):
    """
    Return True, if all members of set2 are contained in set1. 
    i.e, set2 is a subset of set1. See example for clarity:
    Example
    -------
    >>> is_subset(set(["A","B"]), set(["A"]))
    True
    >>> is_subset(set(["A","B"]), set(["A","C"]))
    False
    """
    return set(set2)-set(set1) == set()


def is_almost_subset(set1, set2, min_set_diff = 2):
    """
    Return True, if no more than min_set_diff members of set2 
    are contained in set1. 
    i.e, set2 is almost a subset of set1. See example for clarity:
    
    Example
    -------
    >>>  is_almost_subset(set(["A","B","C","D"]), set(["A", "K"]), 2)
    True
    >>> is_almost_subset(set(["A","B","C","D"]), set(["A", "K"]), 1)
    False

    """
    return len(set(set2)-set(set1)) < min_set_diff


def test_for_subsets(list_of_sets):
    """

    test_for_subsets (formerly known as turtles_all_the_way_down)

    For a ranked list of sets, return a vector where 1 
    indicated the set is not a subset of any of sets
    that come before it in the list. 

    This is useful for eliminating clusters (sets) of TCRs
    which are smaller than a higher ranked and larger set 
    that contains all its members. See example for clarity:

    Example 
    -------
    >>> test_for_subsets([["A","B","C"], ["A","C","D"], ["A","D"], ["B","E"],["B","C"]])
    [1, 1, 0, 1, 0] 
    >>> test_for_subsets([ [1,2,3], [1,3,4], [1,4], [2,5],[2,3]])
    [1, 1, 0, 1, 0]
    """
    tracker = [1]
    if isinstance(list_of_sets, pd.Series):
        list_of_sets = list_of_sets.to_list()
    checked_sets = [list_of_sets[0]]
    for s in list_of_sets[1:]:
        if np.any([is_subset(cs, s) for cs in checked_sets]):
            tracker.append(0)
        else: 
            tracker.append(1)
            checked_sets.append(s)
    assert len(tracker) == len(list_of_sets)
    return tracker


def test_for_almost_subsets(list_of_sets, thr = 3):
    """

    test_for_subsets (formerly known as turtles_all_the_way_down)

    For a ranked list of sets, return a vector where 1 
    indicated the set is not a subset of any of sets
    that come before it in the list. 

    This is useful for eliminating clusters (sets) of TCRs
    which are smaller than a higher ranked and larger set 
    that contains all its members. See example for clarity:

    Example 
    -------
    >>> test_for_almost_subsets([["A","B","C"], ["A","C","D"], ["A","D"], ["B","E"],["B","C"]], 1)
    [1, 1, 0, 1, 0] 
    >>> test_for_almost_subsets([ [1,2,3], [1,3,4], [1,4], [2,5],[2,3]], 1)
    [1, 1, 0, 1, 0]
    """
    tracker = [1]
    if isinstance(list_of_sets, pd.Series):
        list_of_sets = list_of_sets.to_list()
    checked_sets = [list_of_sets[0]]
    for s in list_of_sets[1:]:
        if np.any([is_almost_subset(cs, s, thr) for cs in checked_sets]):
            tracker.append(0)
        else: 
            tracker.append(1)
            checked_sets.append(s)
    assert len(tracker) == len(list_of_sets)
    return tracker



def _dist_summ(data, precision = 1, scientific = True):
    """
    Summarise distribution [as min,q1,median,q3, max]
    
    Parameters
    ----------
    data : list
        List of numeric data
    precision : int
        How many integers precision in scientific notation = 1, 
    scientific : bool
        Default is True, to return result in scientific notation 

    Examples
    --------
    >>> _dist_summ([1,2,3,4,5])
    ['1.e+00', '2.e+00', '3.e+00', '4.e+00', '5.e+00']
    
    _dist_summ([1,2,3,4,5], scientific=False)
    [1, 2.0, 3.0, 4.0, 5]

    """
    dmin = np.min(data)
    dQ1 = np.percentile(data, q = 25, interpolation = 'midpoint') 
    dmedian = np.median(data)
    dQ3 = np.percentile(data, q = 75, interpolation = 'midpoint') 
    dmax = np.max(data)
    r = [dmin, dQ1, dmedian, dQ3, dmax]
    if scientific:
        return [np.format_float_scientific(s, precision = precision) for s in r]
    else:
        return r 



def _select(df, iloc_rows, col = 'cdr3_b_aa'):
    return df.iloc[iloc_rows,][col].to_list()

def _summ(df, indices, column = None , f=None, fdf = None, **kwargs):
    """
    _summ implements a split, apply some function, combine result routine. 
    
    Parameters
    ----------
    f : callable
        a function callable on a list of series
    fdf : callable
        a function callable on a dataframe
    df : pd.DataFrame
        DataFrame
    indices : list
        list of lists containing integers corresponding to the iloc rows of a the < df >
    column : str or None
        column name, should be None if using a fdf
        
    Returns
    -------
    summary : list of identical lenght to indices

    Examples
    --------
    >>> from tcrdist.summarize import _summ, _occurs_N_str, _top_N_str
    >>> df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    >>> _summ(df, indices = [[0,1], [2,3]], column = 'numvar', f = np.median)
    [5.5, 51.5]
    >>> _summ(df, indices = [[0,1], [2,3]], column = 'catvar', f = _occurs_N_str, N = 2)
    ['b (50.0%), a (50.0%)', 'c (50.0%), b (50.0%)']
    >>> _summ(df, indices = [[0,1], [2,3]], column = 'catvar', fdf = _top_N_str, **{'col': 'catvar', 'count_col': 'numvar','N':2})
    ['a (90.9%), b (9.1%)', 'b (97.1%), c (2.9%)']
    """
    summary = list()
    for ind in indices:
        if f is not None:
            if isinstance(df.iloc[ind, ][column], pd.Series):
                selection = df.iloc[ind, ][column].to_list()
            else:
                selection = df.iloc[ind, ][column]
            summary.append(f(selection, **kwargs))
        elif fdf is not None:
            selection = df.iloc[ind, ]
            summary.append(fdf(selection, **kwargs))
        else:
            raise(ValueError("No function (f) or function on a DataFrame (fdf) were supplied\n"))
    assert len(summary) == len(indices)
    return summary

def _occurs_N_str(m, N):
    """
    Return occurances in a pd.Series as a string 

    Example
    -------
    >>> _occurs_N_str(["a","b","b","c"], 1)
    'b (50.0%)'    
    >>> _occurs_N_str(["a","b","b","c"], 2)
    'b (50.0%), c (25.0%)' 
    >>> _occurs_N_str(["a","b","b","c"], 3)
    'b (50.0%), c (25.0%), a (25.0%)'
    """
    if isinstance(m, pd.Series):
        gby = m.value_counts()
    else:
        m = pd.Series(m)
        gby = m.value_counts()
    gby = 100 * gby / gby.sum()
    gby = gby.sort_values(ascending=False)
    out = ', '.join(['%s (%2.1f%%)' % (idx, v) for idx,v in gby.iteritems()][:N])
    return out


def _top_N_str(m, col, count_col, N):
    """
    Example
    -------
    >>> df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    >>> _top_N_str(df, col = 'catvar', count_col ='numvar', N=2)
    'b (88.6%), a (8.8%)'
    """
    gby = m.groupby(col)[count_col].agg(np.sum)
    gby = 100 * gby / gby.sum()
    gby = gby.sort_values(ascending=False)
    out = ', '.join(['%s (%2.1f%%)' % (idx, v) for idx,v in gby.iteritems()][:N])
    return out


def _extract_percentage(s, key):
    """
    extractor for pattern '%s (%2.1f%%)', see examples for clarity

    Parameter
    ---------
    s : str
        string pattern '%s (%2.1f%%)','%s (%2.1f%%)','%s (%2.1f%%)'
    k : str
        key for the percentage you want to extract

    Returns 
    -------
    tuple (str, float)
    
    Examples 
    --------
    >>> _extract_percentage('naive_CD8 (100.0%)', 'naive_CD8')
    ('naive_CD8', '100.0')

    >>> _extract_percentage('naive_CD8 (100.0%)', 'PBMC')
    ('PMBC', 0.0)

    >>> _extract_percentage('naive_CD8 (94.1%), PBMC (5.9%)', 'PBMC')
    ('PBMC', '5.9')
    """
    ls = s.split(",")
    try: 
        rs = [re.search(pattern = '([A-Za-z0-9_]+) [(]([0-9]+[\.][0-9])%[)]', string = s) for s in ls]
        rgs =  [reg.groups() for reg in rs]
        return (key, {k:v for k,v in rgs}[key])
    except:
        return key, 0.0



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