import pytest
import pandas

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
    import numpy as np
    from tcrdist.repertoire import TCRrep
    from tcrsampler.sampler import TCRsampler
    from tcrdist.ecdf import distance_ecdf, make_ecdf_step
    from tcrdist.background import make_gene_usage_counter, make_vj_matched_background, \
                                    make_flat_vj_background, get_gene_frequencies, calculate_adjustment

    import matplotlib.pyplot as plt

    df = pd.read_csv('dash.csv')
    df = df.loc[df['epitope'] == 'PB1']
    tr = TCRrep(cell_df = df, 
                organism = 'mouse', 
                chains = ['beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    TCRsampler.download_background_file(download_file='wiraninha_sampler.zip')
    cols = ['v_b_gene','j_b_gene']

    refs = []
    for ts_fn in [f'wirasinha_mouse_beta_s_{i}.tsv.sampler.tsv' for i in '48']:
        ts = TCRsampler(default_background=ts_fn)
        ts.build_background(stratify_by_subject=True, use_frequency=False)
        
        """Sanitize the alleles to *01 for TCRSampler"""
        tmp = df[cols].applymap(lambda s: s.split('*')[0] + '*01')
        freqs = tmp.groupby(cols).size()
        freq_records = list(freqs.to_frame().to_records())
        ref = ts.sample(freq_records, depth=10, seed=110820)
        ref_df = pd.concat([pd.DataFrame({'cdr3_b_aa':ref[i]}).assign(v_b_gene=v, j_b_gene=j) for i,(v,j,_) in enumerate(freq_records)])        

        """Assigns pV, pJ and pVJ to ref_df"""
        ref_df = get_gene_frequencies(ts=ts, df=ref_df) 
        
        xdf = freqs.reset_index()
        xdf.columns = ['v_b_gene','j_b_gene', 'n']
        
        """For each V,J pairing compute frequency in this reference"""
        xdf = xdf.assign(ref_freq=xdf['n'] / xdf['n'].sum())
        ref_df = ref_df.merge(xdf, how='left', on=cols).reset_index()

        """ Assign weights to ref sequences: Pr_actual / Pr_sampling"""
        ref_df = ref_df.assign(weights=ref_df['pVJ'] / ref_df['ref_freq'])
        refs.append(ref_df)

        """Add uniformly sampled sequences"""
        ref_df = ts.ref_df.sample(100, random_state=1)
        refs.append(ref_df)

    ref_df = pd.concat(refs, axis=0)
    ref_tr = TCRrep(cell_df=ref_df[cols + ['cdr3_b_aa', 'weights']], 
                    organism='mouse', 
                    chains=['beta'],
                    compute_distances=False,
                    store_all_cdr=False)

    tr.compute_rect_distances(df=tr.clone_df, df2=ref_tr.clone_df, store=False)
    
    thresholds = np.arange(1, 50)
    thresholds, ref_ecdf = distance_ecdf(tr.rw_beta,
                                     thresholds=thresholds,
                                     weights=ref_tr.clone_df['weights'] * ref_tr.clone_df['count'])

    thresholds, target_ecdf = distance_ecdf(tr.pw_beta,
                                     thresholds=thresholds,
                                     weights=None)

    figh = plt.figure(figsize=(5, 5))
    axh = figh.add_axes([0.15, 0.15, 0.6, 0.7], yscale='log')
    plt.ylabel(f'Proportion of reference TCRs')
    plt.xlabel(f'Distance from target TCR clone')
    for tari in range(ref_ecdf.shape[0]):
        x, y = make_ecdf_step(thresholds, ref_ecdf[tari, :])
        axh.plot(x, y, color='k', alpha=0.2)
    x, y = make_ecdf_step(thresholds, np.mean(ref_ecdf, axis=0))
    axh.plot(x, y, color='r', alpha=1)

    figh = plt.figure(figsize=(5, 5))
    axh = figh.add_axes([0.15, 0.15, 0.6, 0.7], yscale='log')
    plt.ylabel(f'Proportion of target TCRs')
    plt.xlabel(f'Distance from target TCR clone')
    for tari in range(target_ecdf.shape[0]):
        x, y = make_ecdf_step(thresholds, target_ecdf[tari, :])
        axh.plot(x, y, color='k', alpha=0.2)
    x, y = make_ecdf_step(thresholds, np.mean(target_ecdf, axis=0))
    axh.plot(x, y, color='r', alpha=1)

    """Make an "ROC" plot combining the ECDF against the target (sensitivity)
    vs. ECDF against the reference (specificity)"""
    figh = plt.figure(figsize=(7, 5))
    axh = figh.add_axes([0.15, 0.15, 0.6, 0.7], yscale='log', xscale='log')
    plt.ylabel(f'Proportion of target TCRs')
    plt.xlabel(f'Proportion of reference TCRs')
    for tari in range(target_ecdf.shape[0]):
        x, y = make_ecdf_step(ref_ecdf[tari, :], target_ecdf[tari, :])
        axh.plot(x, y, color='k', alpha=0.2)
    x, y = make_ecdf_step(np.mean(ref_ecdf, axis=0), np.mean(target_ecdf, axis=0))
    axh.plot(x, y, color='r', alpha=1)
    yl = plt.ylim()
    xl = plt.xlim()
    #yl = (1e-6, 0.3)
    plt.plot(yl, yl, '--', color='gray')
    plt.xlim(xl)
    plt.ylim(yl)

