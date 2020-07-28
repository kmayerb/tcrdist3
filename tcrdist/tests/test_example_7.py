import pytest
import Levenshtein
def my_own_metric(s1,s2):   
    return Levenshtein.distance(s1,s2)
def test_example_7():
    """
    If you don't want to use OOP, but you I still want a multi-CDR 
    tcrdistances on a single chain, using you own metric 

    def my_own_metric(s1,s2):   
        return Levenshtein.distance(s1,s2)    
    """
    import multiprocessing
    import pandas as pd
    from tcrdist.rep_funcs import _pws, _pw
    
    df = pd.read_csv("dash2.csv")
    
    metrics_b = {
        "cdr3_b_aa" : my_own_metric ,
        "pmhc_b_aa" : my_own_metric ,
        "cdr2_b_aa" : my_own_metric ,
        "cdr1_b_aa" : my_own_metric }

    weights_b = { 
        "cdr3_b_aa" : 1,
        "pmhc_b_aa" : 1,
        "cdr2_b_aa" : 1,
        "cdr1_b_aa" : 1}

    kargs_b = {  
        'cdr3_b_aa' : 
            {'use_numba': False},
        'pmhc_b_aa' : {
            'use_numba': False},
        'cdr2_b_aa' : {
            'use_numba': False},
        'cdr1_b_aa' : {
            'use_numba': False}
        }

    dmats =  _pws(df = df , 
                metrics = metrics_b, 
                weights = weights_b, 
                kargs   = kargs_b , 
                cpu     = 1, 
                uniquify= True, 
                store   = True)

    print(dmats.keys())