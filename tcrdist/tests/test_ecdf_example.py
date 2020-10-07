import pytest

def test_dash_ecdf():
    """
    An empirical distribution function (ECDF) can be created
    for a target TCR and a reference set of TCRs to show
    the proportion of reference TCRs that are within a distance
    D of the target TCR, over a range of distances.

    A plot of the ECDF as a function of increasing D shows the
    density of TCR space in the reference set in the neighborhood
    around the target TCR. This can be very helpful for 
    identifying dense antigen-specific clusters in an antigen
    enriched TCR repertoire, where the "reference" set is 
    actually an experimentally enriched repertoire (e.g. 
    pMHC:tetramer or AIM sorting). Or the ECDF can be helpful
    for identifying a radius around a TCR that retains high
    antigen specificity, by showing that the neighborhood
    is extremely sparse in an large unsorted/bulk TCR repertoire.
    
    """
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrsampler.sampler import TCRsampler
    from tcrdist.ecdf import distance_ecdf, make_ecdf_step

    df = pd.read_csv('dash.csv')
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    TCRsampler.download_background_file(download_file = 'ruggiero_mouse_sampler.zip')
    ts = TCRsampler(default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv')
    ts.build_background(stratify_by_subject = True, use_frequency = False)
    
    tmp = df[['v_b_gene', 'j_b_gene']].applymap(lambda s: s.split('*')[0] + '*01')

    freqs = tmp.groupby(['v_b_gene', 'j_b_gene']).agg(lambda v: v.shape[0])
    freqs = list(freqs.to_frame().to_records())
    ref = ts.sample(freqs, depth=100, seed=110820)
    ref_df = pd.concat([pd.DataFrame({'cdr3_b_aa':ref[i]}).assign(v_b_gene=v, j_b_gene=j) for i,(v,j,_) in enumerate(freqs)])
    ref_df.loc[:, 'count'] = 1
        
    ref_tr = TCRrep(cell_df=ref_df, 
                    organism='mouse', 
                    chains=['beta'],
                    compute_distances=False,
                    store_all_cdr=False)

    tr.compute_rect_distances(df=tr.clone_df, df2=ref_tr.clone_df, store=False)
    
    """TODO: Add weights to correct for sampling. Also add randomly sampled sequences"""
    thresholds, ecdf = distance_ecdf(tr.rw_beta, thresholds=None, weights=None)

    figh = plt.figure(figsize=(5, 5))
    axh = figh.add_axes([0.15, 0.15, 0.6, 0.7], yscale='log')
    plt.ylabel(f'Proportion of reference TCRs')
    plt.xlabel(f'Distance from target TCR clone')

    for tari in range(ecdf.shape[0]):
        x, y = make_ecdf_step(thresholds, ecdf[tari, :])
        axh.plot(x, y, color=k, alpha=0.2)
    
    x, y = make_ecdf_step(thresholds, np.mean(ecdf, axis=0))
    axh.plot(x, y, color='r', alpha=1)