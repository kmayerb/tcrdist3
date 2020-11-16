



def test_sparse_example_report_qpublic():
    """
    Making a meta-clonotype report from a 
    scipy.sparse TCRdist matrix.
    """
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.public import _neighbors_sparse_fixed_radius, _neighbors_sparse_variable_radius
    from tcrdist.summarize import test_for_subsets
    from tcrdist.tree import _default_sampler_olga
    from tcrdist.public import make_motif_logo_from_index

    df = pd.read_csv("dash.csv").query('epitope == "PA"')
    tr = TCRrep(cell_df = df,               #(2)
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False)
        # When setting the radius to 50, the sparse matrix 
        # will convert any value > 50 to 0. True zeros are 
        # repressented as -1.
    radius = 50
    tr.cpus = 1
        # Notice that we called .compute_sparse_rect_distances instead of .compute_distances
    tr.compute_sparse_rect_distances(df = tr.clone_df, radius = radius)

        # There are two functions for finding neighbors from a sparse TCRdist matrix. 
        # For 1 fixed radius: _neighbors_sparse_fixed_radius()
        # For a radius per row: _neighbors_sparse_variable_radius()
    tr.clone_df['radius'] = 12 
    tr.clone_df['neighbors'] = \
        _neighbors_sparse_variable_radius(
            csrmat = tr.rw_beta, 
            #radius = 12)
            radius_list = tr.clone_df['radius'].to_list())

        # <K_neighbors>the number of neighbors per TCR
    tr.clone_df['K_neighbors'] = tr.clone_df['neighbors'].apply(lambda x: len(x))
        # <nsubject> the number of subject (nsubject) neighboring the TCR (
    tr.clone_df['nsubject'] = tr.clone_df['neighbors'].apply(lambda x: len(tr.clone_df['subject'][x].unique()))
        # nsubject > 1 implies quasi-publicity)
    tr.clone_df['qpublic'] = tr.clone_df['nsubject'].apply(lambda x: x >1 )
    
        # e.g., For the report, focus on TCRs with more than 5 neighboring subjects 
    quasi_public_df = tr.clone_df.query('nsubject > 5').copy().\
        sort_values('nsubject', ascending = False)
        # test_for_subsets()> allows us to remove TCRs with identical neighborhoods
    quasi_public_df['unique']  = test_for_subsets(quasi_public_df['neighbors'])
    quasi_public_df = quasi_public_df[quasi_public_df['unique'] == 1].copy()
        # declare a sampler for generating a backgrond comparison
    ts = _default_sampler_olga(organism = 'mouse', chain = 'beta')()

        # make a background-subtracted logo <svg> and raw log <svg_raw> for each TCR
    svgs = list()
    svgs_raw = list()
    for i,r in quasi_public_df.iterrows():
        svg, svg_raw  = make_motif_logo_from_index(tcrsampler = ts,
                                                   ind = r['neighbors'],
                                                   centroid = r['cdr3_b_aa'],
                                                   clone_df = tr.clone_df,
                                                   cdr3_name = 'cdr3_b_aa',
                                                   v_name = 'v_b_gene',
                                                   gene_names = ['v_b_gene','j_b_gene'])
        svgs.append(svg)
        svgs_raw.append(svg_raw)

        # Output a html report
    output_html_name = 'quasi_public_df_report.html'
    quasi_public_df['svg'] = svgs
    quasi_public_df['svg_raw'] = svgs_raw
        # Specific columns to include in the report
    labels = [  'cdr3_b_aa', 
                'v_b_gene',
                'j_b_gene',
                'radius', 
                'K_neighbors',
                'nsubject']

    def shrink(s):
        """reduce size of svg graphic"""
        s = s.replace('height="100%"', 'height="20%"')
        s = s.replace('width="100%"', 'width="20%"')
        return s

    with open(output_html_name, 'w') as output_handle:
        for i,r in quasi_public_df.iterrows():
            #import pdb; pdb.set_trace()
            svg, svg_raw = r['svg'],r['svg_raw']
            output_handle.write("<br></br>")
            output_handle.write(shrink(svg))
            output_handle.write(shrink(svg_raw))
            output_handle.write("<br></br>")
            output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
            output_handle.write("<br></br>")