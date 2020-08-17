import pytest 



def test_pgen_1():
    """
    How to add pgen estimates to human alpha/beta CDR3s
    """
    import pandas as pd
    from tcrdist.pgen import OlgaModel
    from tcrdist import mappers 
    from tcrdist.repertoire import TCRrep
    from tcrdist.setup_tests import download_and_extract_zip_file

    df = pd.read_csv("dash_human.csv")

    tr = TCRrep(cell_df = df.sample(5, random_state = 3), 
                organism = 'human', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv', 
                store_all_cdr = False)

    olga_beta  = OlgaModel(chain_folder = "human_T_beta", recomb_type="VDJ")
    olga_alpha = OlgaModel(chain_folder = "human_T_alpha", recomb_type="VJ")

    tr.clone_df['pgen_cdr3_b_aa'] = olga_beta.compute_aa_cdr3_pgens(
        CDR3_seq = tr.clone_df.cdr3_b_aa)
    
    tr.clone_df['pgen_cdr3_a_aa'] = olga_alpha.compute_aa_cdr3_pgens(
        CDR3_seq = tr.clone_df.cdr3_a_aa)

    tr.clone_df[['cdr3_b_aa', 'pgen_cdr3_b_aa', 'cdr3_a_aa','pgen_cdr3_a_aa']]
    """
                cdr3_b_aa  pgen_cdr3_b_aa          cdr3_a_aa  pgen_cdr3_a_aa
    0    CASSETSGRSPYEQYF    1.922623e-09    CAVRPGYSSASKIIF    1.028741e-06
    1  CASSPGLASPYSYNEQFF    1.554569e-11    CAVRVLMEYGNKLVF    1.128865e-08
    2       CASSSRSTDTQYF    5.184760e-07     CAGAGGGSQGNLIF    2.512580e-06
    3        CASSIGVYGYTF    8.576919e-08  CAFMSNAGGTSYGKLTF    1.832054e-07
    4  CASSLLVSGVSSTDTQYF    2.143508e-11       CAVIGEGGKLTF    5.698448e-12
    """