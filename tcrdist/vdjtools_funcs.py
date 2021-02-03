import pandas as pd
from tcrdist.adpt_funcs import _valid_cdr3
from tcrdist.paths import path_to_base
from tcrdist.repertoire_db import RefGeneSet
import warnings
import os

def import_vdjtools(   vdj_tools_file,
                        chain = 'beta',
                        organism = 'human',
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        validate = True):
    """
    Import VDJtools formated input .tsv or .tsv.gz file. 
    Input file must have columns ['count', 'freq', 'cdr3aa', 'v', 'j','cdr3nt'], 
    see (https://vdjtools-doc.readthedocs.io/en/master/input.html#vdjtools-format)
    for examples of how you can use VDJtool to directly import 
    and convert from the following: 
    MiTCR, MiGEC, IgBlast (MIGMAP), ImmunoSEQ, VDJdb, Vidjil, MiXCR

    Parameters
    -----------
    vdj_tools_file : str
        e.g., os.path.join(path_to_base, 'tcrdist','data','formats','vdj.M_15_CD8_beta.clonotypes.TRB.txt.gz')
    chain : str
        'alpha','beta','gamma', or 'delta'
    organism : str
        'human', 'mouse'
    db_file : str
        'alphabeta_gammadelta_db.tsv',
    validate : bool 
        If True, only clones with valid CDR3AA, V-Gene and J-Gene are returned)
        CDR3 length must be >= 5 for tcrdist3 
    include_nucseq : True
        If True, retain nucletoide sequence

    Returns
    -------
    df or valid_df : DataFrame
        with columns : ['count', 'freq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_nucseq','valid_v', 'valid_j', 'valid_cdr3']
    
    """

    assert chain in ['alpha','beta','gamma','delta']
    assert organism in ['human','mouse']
    
    _chain  = { 'alpha' : {'cdr3aa':'cdr3_a_aa','v':'v_a_gene', 'j':'j_a_gene', 'cdr3nt':'cdr3_a_nucseq', 'AB' : 'A'},
                'beta' :  {'cdr3aa':'cdr3_b_aa','v':'v_b_gene', 'j':'j_b_gene', 'cdr3nt':'cdr3_b_nucseq', 'AB' : 'B'},
                'gamma' : {'cdr3aa':'cdr3_g_aa','v':'v_g_gene', 'j':'j_g_gene', 'cdr3nt':'cdr3_g_nucseq', 'AB' : 'A'},
                'delta' : {'cdr3aa':'cdr3_d_aa','v':'v_d_gene', 'j':'j_d_gene', 'cdr3nt':'cdr3_d_nucseq', 'AB' : 'B'} }.\
                get(chain)
    
    
    df = pd.read_csv(vdj_tools_file, sep = "\t")
    
    df = df[['count', 'freq', 'cdr3aa', 'v', 'j','cdr3nt']].\
        rename(columns = {
            'cdr3aa': _chain['cdr3aa'],
            'v'     : _chain['v'],
            'j'     : _chain['j'],
            'cdr3nt': _chain['cdr3nt']
            })

    df[_chain['v']] = df[_chain['v']].apply(lambda x: f"{x}*01")
    df[_chain['j']] = df[_chain['j']].apply(lambda x: f"{x}*01")
        # CHECK AND ADD ALELE
    all_genes =  RefGeneSet(db_file = db_file).all_genes[organism]
    all_valid_genes = [x for x in all_genes.keys() if all_genes[x].chain == _chain['AB']]

    warnings.warn(f"ADDING *01 allele to each V- and J-gene") 
    df['valid_v'] = df[_chain['v']].apply(lambda x : all_genes.get(x) is not None )
    df['valid_j'] = df[_chain['j']].apply(lambda x : all_genes.get(x) is not None )
    df['valid_cdr3'] = df[_chain['cdr3aa']].apply(lambda x : _valid_cdr3(x) and len(x) >=5 )

    x1 = df.shape[0]
    invalid_df = df[~(df['valid_v'] & df['valid_j'] & df['valid_cdr3'])].copy()
    valid_df   = df[(df['valid_v'] & df['valid_j'] & df['valid_cdr3'])].reset_index(drop = True).copy()
    xiv = invalid_df.shape[0]
    warnings.warn(f"{xiv} of {x1} had invalid {_chain['v']}, {_chain['j']}, or {_chain['cdr3aa']}")
    
    if validate:
        warnings.warn(f"REMOVED {xiv} of {x1} with invalid {_chain['v']}, {_chain['j']}, or {_chain['cdr3aa']}")
        return valid_df
    else:
        warnings.warn(f"Invalid clones were not removed, do so manually before proceeding")
        return df
        





