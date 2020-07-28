import pytest
import Levenshtein
def my_own_metric(s1,s2):   
    return Levenshtein.distance(s1,s2)
def test_example_6():
    """
    If you hate object oriented programming, just show me the functions. 
    No problem. 
    
    Maybe you only care about the CDR3 on the beta chain.

    def my_own_metric(s1,s2):   
        return Levenshtein.distance(s1,s2)
    """  
    import multiprocessing
    import pandas as pd
    from tcrdist.rep_funcs import _pws, _pw

    df = pd.read_csv("dash2.csv")

    # 
    dmat = _pw( metric = my_own_metric,
                seqs1 = df['cdr3_b_aa'].values,
                ncpus=2,
                uniqify=True,
                use_numba=False)
