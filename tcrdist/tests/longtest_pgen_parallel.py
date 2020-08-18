import pytest



def longtest_pgen_with_parmap():
    """
    Test speed up of computation of many pgens using parmap to make use of more than one cpu

    For 1000 CDR3 

    Finished 'olga_in_series' in 26.3842 secs with 1 core
    Finished 'olga_in_parmap' in 6.2384 secs with 6 cores
    """
    import numpy as np
    import pandas as pd
    import parmap
    from tcrdist.pgen import OlgaModel
    from tcrdist.speed import timer
    from tcrdist.adpt_funcs import _valid_cdr3

    from tcrdist.setup_tests import download_and_extract_zip_file
    download_and_extract_zip_file( 'cdr3_beta_500K.zip', source = "dropbox", dest = ".")

    olga_beta  = OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ")
    
    n = 1000
    df = pd.read_csv('cdr3_beta_500K.csv')
    inputlist = df.iloc[:,0].to_list()[0:n]
    inputlist  = [x for x in inputlist if _valid_cdr3(x)]

    @timer
    def olga_in_series(f = olga_beta.compute_aa_cdr3_pgen, input = inputlist):
        return [f(x) for x in input]

    @timer
    def olga_in_parmap(f = olga_beta.compute_aa_cdr3_pgen, input = inputlist, **kwargs):
        return parmap.map(f, input, pm_pbar=True, **kwargs)

    r1 = olga_in_series(f = olga_beta.compute_aa_cdr3_pgen, input = inputlist)
    r2 = olga_in_parmap(f = olga_beta.compute_aa_cdr3_pgen, input = inputlist)
    
    assert np.all(r1 == r2)


def test_pgen_with_parmap():
    """
    Really simple example of using multiple cpus to 
    speed up computation of pgens with olga.
    """
    import parmap
    from tcrdist.pgen import OlgaModel
    olga_beta  = OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ")
    parmap.map(olga_beta.compute_aa_cdr3_pgen, ['CASSYRVGTDTQYF', 'CATSTNRGGTPADTQYF', 'CASQGDSFNSPLHF', 'CASSPWTGSMALHF'])