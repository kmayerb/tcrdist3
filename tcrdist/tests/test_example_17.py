"""
example 17
"""

def test_example_report_public():
    """
    tcrdist3 is particularly useful for finding 
    what we term quasi-public meta-clonotypes, 
    collections of biochemically similar TCRs 
    recognizing the same peptide-MHC. 

    The easist way to define meta-clonotypes
    is to compute pairwise distances between 
    TCRs found in an antigen-enriched 
    subrepertoire, abbreviated below as 
    <aesr>
    """
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic
    
        # <aesr_fn> antigen-enriched subrepertoire
    fn = os.path.join('tcrdist', 'data','covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv')
        # <aesr_df> antigen-enriched subrepertoire
    df = pd.read_csv(fn)
        # <tr> TCR repertoire
    tr = TCRrep(
        cell_df = df[['cohort','subject','v_b_gene', 'j_b_gene','cdr3_b_aa']].copy(), 
        organism = 'human', 
        chains = ['beta'], 
        db_file = 'alphabeta_gammadelta_db.tsv', 
        compute_distances = True)

        # <tp> TCRpublic class for reporting publicities, fixed radius 18, 'nsubject > 3'
    tp = TCRpublic(
        tcrrep = tr, 
        output_html_name = "quasi_public_clones.html")

        # by calling, .report() an html report is made
    public = tp.report()
        
        # Also, the following datafarme are available
        # <clone_df> pd.DataFrame clone_df from tr.clone_df 
        # with neighbors and summary measures appended
    public['clone_df']
        # <nn_summary> pd.DataFrame with just summary measures
    public['nn_summary']
        # <quasi_public_df> Non-redundant groups of quasipublic clones
    public['quasi_public_df']
