import pytest 


def test_fuzzy_div():
    import numpy as np
    from scipy.spatial.distance import squareform
    from tcrdist.diversity import fuzzy_diversity

    n = 100
    counts = np.random.randint(1, 20, size=n)
    dvec = np.round(np.random.rand((n**2 - n) // 2))
    pwmat = squareform(dvec)

    two = fuzzy_diversity(counts, pwmat, order=2, threshold=1)
    two_s = fuzzy_diversity(counts, pwmat, order=2, threshold=1, nsamples=10000, force_sampling=True)
    three = fuzzy_diversity(counts, pwmat, order=3, threshold=1, nsamples=10000)
    four = fuzzy_diversity(counts, pwmat, order=4, threshold=1, nsamples=10000)

    print(two, two_s, three, four)


def test_diversity_example_1():
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist.diversity import generalized_simpsons_entropy, simpsons_difference, fuzzy_diversity

    df = pd.read_csv("dash.csv")

    """Compute the Simpson's diversity index (SDI) with 95% CI
    and the "effective number" of TCRs for a TCR repertoire."""
    pb1_div = generalized_simpsons_entropy(df.loc[df['epitope'] == 'PB1', 'count'])

    """The most common SDI is order=2. But see Grabchak et al.
    who generalized the index/entropy for all integer orders and show it is helpful
    to look at a range of orders:

    Grabchak M, Marcon E, Lang G, Zhang Z (2017) The generalized Simpson’s entropy
        is a measure of biodiversity. PLoS ONE 12(3): e0173305.
        https://doi.org/10.1371/journal.pone.0173305"""
    pb1_div = generalized_simpsons_entropy(df.loc[df['epitope'] == 'PB1', 'count'],
                                           orders=np.arange(2, 20),
                                           alpha=0.05)

    """Compare the SDI for two repertoires"""
    pa_div = generalized_simpsons_entropy(df.loc[df['epitope'] == 'PA', 'count'],
                                           orders=np.arange(2, 20),
                                           alpha=0.05)

    div_diff = simpsons_difference(df.loc[df['epitope'] == 'PB1', 'count'],
                                    df.loc[df['epitope'] == 'PA', 'count'],
                                    orders=np.arange(2, 20))

    """Compute a "fuzzy" Simpson's diversity index using the pairwise distance matrix."""
    tr = TCRrep(cell_df = df.loc[df['epitope'] == 'PB1'], 
                organism = 'mouse', 
                chains = ['alpha','beta'], 
                db_file = 'alphabeta_gammadelta_db.tsv')

    fdiv = fuzzy_diversity(tr.clone_df['count'],
                           tr.pw_alpha + tr.pw_beta,
                           order=2,
                           threshold=75)