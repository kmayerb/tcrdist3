"""
An example for the gallery.
"""

def test_gallery_hdiff():
    """
    All imports are provided here, and are repeated 
    step-wise below, for clarity, and for
    module cut-and-paste. This example
    performs paired alpha-beta analysis,
    but code blocks can be used for single
    chain analysis as well.
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.rep_diff import hcluster_diff, member_summ
    from tcrsampler.sampler import TCRsampler
    from tcrdist.adpt_funcs import get_centroid_seq
    from tcrdist.summarize import _select
    from palmotif import compute_pal_motif, svg_logo
    from hierdiff import plot_hclust_props
    """
    Load a subset of data that contains paired alpha-beta
    chain mouse TCR receptors that recognized 
    the PA or PB1 epitopes (present in mouse influenza). 
    """
    import pandas as pd
    df = pd.read_csv("dash.csv")
    conditional = df['epitope'].apply( lambda x: x in ['PA','PB1'])
    """
    For illustrative/testing purposes, randomly subset the data to include 
    only 100 clones. Increase for more informative plot.
    """
    df = df[conditional].\
        reset_index(drop = True).\
        sample(100, random_state = 3).\
        reset_index(drop = True).\
        copy()
    """
    Load DataFrame into TCRrep instance, 
    which automatically computes attributes:
    1. .clone_df DataFrame
    2. .pw_beta nd.array 
    3. .pw_alpha nd.array 
    """
    from tcrdist.repertoire import TCRrep
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta','alpha'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    """
    Apply hcluster_diff, which hierarchically clusters.
    
    Note
    ----
    pwmat could easily be tr.pw_beta or tr.pw_alpha if 
    clustering should be done on a single chain.
    """
    from tcrdist.rep_diff import hcluster_diff
    tr.hcluster_df, tr.Z =\
        hcluster_diff(clone_df = tr.clone_df, 
                      pwmat    = tr.pw_beta + tr.pw_alpha,
                      x_cols = ['epitope'], 
                      count_col = 'count')

    """
    Load a custom background, mouse appropriate dataset to sample CDR3s 
    according to the V and J gene usage frequencies observed in each node.
    See the tcrsampler package for more details 
    (https://github.com/kmayerb/tcrsampler/blob/master/docs/getting_default_backgrounds.md)
    """
    from tcrsampler.sampler import TCRsampler

    t = TCRsampler()
    t.download_background_file("ruggiero_mouse_sampler.zip")
    tcrsampler_beta = TCRsampler(default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv')
    tcrsampler_alpha = TCRsampler(default_background = 'ruggiero_mouse_alpha_t.tsv.sampler.tsv')

    """
    Add an SVG graphic to every node of the tree 
    aligned to the cluster centroid.
    """
    from tcrdist.adpt_funcs import get_centroid_seq
    from tcrdist.summarize import _select
    from palmotif import compute_pal_motif, svg_logo

    """Beta Chain"""
    svgs_beta = list()
    for i,r in tr.hcluster_df.iterrows():

        dfnode = tr.clone_df.iloc[r['neighbors_i'],]
        if dfnode.shape[0] > 2:
            centroid, *_ = get_centroid_seq(df = dfnode)
        else:
            centroid = dfnode['cdr3_b_aa'].to_list()[0]
        print(f"BETA-CHAIN: {centroid}")

        gene_usage_beta = dfnode.groupby(['v_b_gene','j_b_gene']).size()
        sampled_rep = tcrsampler_beta.sample( gene_usage_beta.reset_index().to_dict('split')['data'],
                        flatten = True, depth = 10)
        sampled_rep  = [x for x in sampled_rep if x is not None]
        motif, stat = compute_pal_motif(
                        seqs = _select(df = tr.clone_df, 
                                       iloc_rows = r['neighbors_i'], 
                                       col = 'cdr3_b_aa'),
                        refs = sampled_rep, 
                        centroid = centroid)
        
        svgs_beta.append(svg_logo(motif, return_str= True))

    """Add Beta SVG graphics to hcluster_df"""
    tr.hcluster_df['svg_beta'] = svgs_beta


    """Alpha Chain"""
    svgs_alpha = list()
    for i,r in tr.hcluster_df.iterrows():

        dfnode = tr.clone_df.iloc[r['neighbors_i'],]
        if dfnode.shape[0] > 2:
            centroid, *_ = get_centroid_seq(df = dfnode)
        else:
            centroid = dfnode['cdr3_a_aa'].to_list()[0]
        print(f"ALPHA-CHAIN: {centroid}")
        gene_usage_alpha = dfnode.groupby(['v_a_gene','j_a_gene']).size()
        sampled_rep = tcrsampler_alpha.sample( gene_usage_alpha.reset_index().to_dict('split')['data'], 
                        flatten = True, depth = 10)
        
        sampled_rep  = [x for x in sampled_rep if x is not None]
        motif, stat = compute_pal_motif(
                        seqs = _select(df = tr.clone_df, 
                                       iloc_rows = r['neighbors_i'], 
                                       col = 'cdr3_a_aa'),
                        refs = sampled_rep, 
                        centroid = centroid)

        svgs_alpha.append(svg_logo(motif, return_str= True))
    
    """Add Alpha SVG graphics to hcluster_df"""
    tr.hcluster_df['svg_alpha'] = svgs_alpha
    """
    Produce summary information for tooltips. 
    For instance, describe percentage of TCRs with 
    a given epitope at a given node.
    """
    res_summary = member_summ(  res_df = tr.hcluster_df,
                                clone_df = tr.clone_df, 
                                addl_cols=['epitope'])

    tr.hcluster_df_detailed = \
        pd.concat([tr.hcluster_df, res_summary], axis = 1)
    """
    Write D3 html for interactive denogram graphic. 
    Specify desired tooltips.
    """
    from hierdiff import plot_hclust_props
    html = plot_hclust_props(tr.Z,
                title='PA Epitope Example',
                res=tr.hcluster_df_detailed,
                tooltip_cols=['cdr3_b_aa','v_b_gene', 'j_b_gene','svg_alpha','svg_beta'],
                alpha=0.00001, colors = ['blue','gray'],
                alpha_col='pvalue')

    with open('hierdiff_example_PA_v_PB1.html', 'w') as fh:
        fh.write(html)