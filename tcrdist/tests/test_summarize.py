"""
August 14, 2020
Unit tests for tcrdist.summarize
"""
import pytest
import numpy as np
import pandas as pd
from tcrdist.summarize import _summ, _occurs_N_str, _top_N_str

def test__top_N_str():
    df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    assert _top_N_str(m = df, col = 'catvar', count_col ='numvar', N=2) == 'b (88.6%), a (8.8%)'

def test__occurs_N_str():
    assert _occurs_N_str(["a","b","b","c"], 1) == 'b (50.0%)'    
    assert _occurs_N_str(["a","b","b","c"], 2) == 'b (50.0%), c (25.0%)' or _occurs_N_str(["a","b","b","c"], 2) == 'b (50.0%), a (25.0%)' 
    assert _occurs_N_str(["a","b","b","c"], 3) == 'b (50.0%), c (25.0%), a (25.0%)' or 'b (50.0%), a (25.0%), c (25.0%)'
 
def test__summ_0():
    """
    ValueError: No function (f) or function on a DataFrame (fdf) were supplied
    """
    df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    with pytest.raises(ValueError):
        r =  _summ(df, indices = [[0,1], [2,3]], column = 'numvar')

def test__summ_1():
    df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    r =  _summ(df, indices = [[0,1], [2,3]], column = 'numvar', f = np.median)
    e = [5.5, 51.5]
    assert np.all(r == e)

def test__summ_2():
    df = pd.DataFrame({'catvar':["a","b","b","b","b","c"], "numvar":[10,1,1,1,100,3]})
    r = _summ(df, indices = [[0,1,3], [3,4,5]], column = 'catvar', f = _occurs_N_str, N = 2)
    e = ['b (66.7%), a (33.3%)', 'b (66.7%), c (33.3%)']
    assert np.all(r == e)

def test__summ_3():
    df = pd.DataFrame({'catvar':["a","b","b","c"], "numvar":[10,1,100,3]})
    r = _summ(df, indices = [[0,1], [2,3]], column = 'catvar', fdf = _top_N_str, **{'col': 'catvar', 'count_col': 'numvar','N':2})
    e = ['a (90.9%), b (9.1%)', 'b (97.1%), c (2.9%)']
    assert np.all(r == e)

