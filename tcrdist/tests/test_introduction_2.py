import pytest 



def test_introduction_2():
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.adpt_funcs import import_adaptive_file, adaptive_to_imgt
    
    df = import_adaptive_file(adaptive_filename = 'Adaptive2020.tsv')

    tr = TCRrep(cell_df = df, 
                organism = 'human', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    """Lookup *01 IMGT allele corresponding with an Adaptive gene name"""
    assert adaptive_to_imgt['human']['TCRBV30'] == 'TRBV30*01'    