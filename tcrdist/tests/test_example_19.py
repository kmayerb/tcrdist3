



def test_example_variable_radius():
    """
    Instead of enforcing a fixed radius, 
    use a radius specific to each
    centroid, specified in an additional 
    column.
    """
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic
    
    fn = os.path.join(
        'tcrdist',
        'data',
        'covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.radius.csv')

    df = pd.read_csv(fn)

    tr = TCRrep(cell_df = df[['cohort','subject','v_b_gene', 'j_b_gene','cdr3_b_aa', 'radius']], 
                organism = "human", 
                chains = ["beta"])

    tp = TCRpublic(
        tcrrep = tr, 
        output_html_name = "quasi_public_clones3.html")

    # set to True, if we want a universal radius
    tp.fixed_radius = False
    # must then specify maximum distance for finding similar TCRs
    tp.radius = None
    # set criteria for being quasi-public
    tp.query_str = 'nsubject > 5'
    # Add additional columns to be summarized in the report
    tp.kargs_member_summ['addl_cols'] = ['subject', 'cohort']
    # Add cohort.summary to the labels column so it shows up in the report
    tp.labels.append("cohort.summary")
    tp.labels.append("cdr3s")
    # Change number of subjects to display
    tp.kargs_member_summ['addl_n'] = 10
    # by calling, .report() an html report is made
    public = tp.report()













    # # <aesr_fn> antigen-enriched subrepertoire
    # fn = os.path.join(
    #     'tcrdist',
    #     'data',
    #     'covid19',
    #     'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv')

    # df = pd.read_csv(fn)
    
    # df = df[['subject',
    #     'cohort',
    #     'v_b_gene',
    #     'j_b_gene',
    #     'cdr3_b_aa']]

    # fn2 = os.path.join(
    #     'tcrdist',
    #     'data',
    #     'covid19',
    #     'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.bE5ctrl.centers.csv')
    # df2 = pd.read_csv(fn2)
    # df2 = df2[[
    #     'cdr3_b_aa',
    #     'v_b_gene',
    #     'j_b_gene',
    #     'pgen', 
    #     'max_radi',
    #     'chi2joint']].\
    #     rename(columns = {'max_radi':'radius'}).\
    #     groupby(['v_b_gene','j_b_gene','cdr3_b_aa']).\
    #     first().reset_index(drop = False)
    # outfn = os.path.join(
    #     'tcrdist',
    #     'data',
    #     'covid19',
    #     'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.radius.csv')
    # df.merge(df2, how = "left", on = [ 'cdr3_b_aa','v_b_gene','j_b_gene']).to_csv(outfn, index = False)

