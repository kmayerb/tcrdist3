"""
Simple set of functions for summarizing over a group
"""

import pandas as pd
import numpy as np
import re

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
            selection = df.iloc[ind, ][column].to_list()
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