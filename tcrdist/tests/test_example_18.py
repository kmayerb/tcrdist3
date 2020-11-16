"""
example 18
"""

def test_example_report_public_adjustin_attributes ():
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic
    
    # <aesr_fn> antigen-enriched subrepertoire
    aesr_fn = os.path.join(
        'tcrdist',
        'data',
        'covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv')
    
    # <aesr_df> antigen-enriched subrepertoire
    aesr_df = pd.read_csv(aesr_fn)
    # <tr> TCR repertoire
    tr = TCRrep(
        cell_df = aesr_df[[
            'cohort',
            'subject',
            'v_b_gene', 
            'j_b_gene',
            'cdr3_b_aa']].copy(), 
        organism = 'human', 
        chains = ['beta'], 
        db_file = 'alphabeta_gammadelta_db.tsv', 
        compute_distances = True)
    # <tp> TCRpublic class for reporting publicities 
    tp = TCRpublic(
        tcrrep = tr, 
        output_html_name = "quasi_public_clones2.html")
    # set to True, if we want a universal radius
    tp.fixed_radius = True
    # must then specify maximum distance for finding similar TCRs
    tp.radius = 18
    # set criteria for being quasi-public
    tp.query_str = 'nsubject > 6'
    # Add additional columns to be summarized in the report
    tp.kargs_member_summ['addl_cols'] = ['subject', 'cohort']
    # Add cohort.summary to the labels column so it shows up in the report
    tp.labels.append("cohort.summary")
    # by calling, .report() an html report is made
    public = tp.report()
