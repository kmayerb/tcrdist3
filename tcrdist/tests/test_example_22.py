



def test_ranked_centers_example():
    """
    Report of meta-clonotypes using two dataframes.
    <df>  has all TCRS
    <df2> has a subset of TCRS in <df>, specifiyint which 
    are to be used as centroids.
    """
    import os
    import pandas as pd
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import TCRpublic  
    from tcrdist.tree import _default_sampler_olga
    from progress.bar import IncrementalBar
    from palmotif import compute_pal_motif, svg_logo
    from tcrdist.public import make_motif_logo

    output_html_name = "custom_report.html"
    # <fn> filename for all TCRs in an antigen-enriched repertoire
    fn = os.path.join('tcrdist','data','covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.bE5ctrl.centers.csv')
    df = pd.read_csv(fn, sep = ",")
    df = df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen', 'max_radi']].\
        rename(columns= {'max_radi':'radius'}).copy()

    # <fn>2 filename for priority TCRs
    fn2 = os.path.join('tcrdist','data','covid19',
        'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv.bE5ctrl.centers.csv.ranked_centers.tsv')
    df2 = pd.read_csv(fn2, sep = "\t").\
        rename(columns= {'max_radi':'radius'}).copy()

    # Compute distances between all TCRs
    tr = TCRrep(cell_df = df, 
        organism = 'human',
        chains = ['beta'], 
        compute_distances = True)

    # Initialize a tcrsampler, this will be used to make background motifs
    tcrsampler = _default_sampler_olga(organism = "human", chain = "beta")()

    # Iterate through each row of the df2, making a logo for each.
    svgs = list()
    svgs_raw = list()
    bar = IncrementalBar("Making Logos", max = df2.shape[0])
    for i,r in df2.iterrows():
        bar.next()
        svg,svg_raw=make_motif_logo(tcrsampler = tcrsampler,
                            clone_df = tr.clone_df,
                            pwmat = tr.pw_beta,
                            centroid = r['cdr3_b_aa'],
                            v_gene = r['v_b_gene'],
                            radius = r['radius'],
                            pwmat_str = 'pw_beta',
                            cdr3_name = 'cdr3_b_aa',
                            v_name = 'v_b_gene',
                            gene_names = ['v_b_gene','j_b_gene'])
        svgs.append(svg)
        svgs_raw .append(svg_raw)
    bar.next(); bar.finish()
    df2['svg'] = svgs
    df2['svg_raw'] = svgs_raw

    def shrink(s):
        """reduce size of svg graphic"""
        s = s.replace('height="100%"', 'height="20%"')
        s = s.replace('width="100%"', 'width="20%"')
        return s

    # Choose columns to include in the report
    labels = [  'cdr3_b_aa', 
                'v_b_gene',
                'j_b_gene',
                'radius',
                'regex', 
                'target_hits',
                'nsubject',
                'chi2joint']

    with open(output_html_name, 'w') as output_handle:
        for i,r in df2.iterrows():
            #import pdb; pdb.set_trace()
            svg, svg_raw = r['svg'],r['svg_raw']
            output_handle.write("<br></br>")
            output_handle.write(shrink(svg))
            output_handle.write(shrink(svg_raw))
            output_handle.write("<br></br>")
            output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
            output_handle.write("<br></br>")


