"""
test_breadth.py 

The breadth and depth of a disease-specific T-cell response.

This module concerns the estimation of clonal breadth, whether it be at the 
pathogen, protein, or epitope level. Once meta-clonotype have been defined they 
can be used to search for biochemically similar TCRs in bulk repertoires 
that are likely to share antigen recognition. It is possible that a
single TCR clonotype may be conformant with multiple TCR meta-clonotypes, 
so an accurate estimate of clonal breadth must avoid double counting such 
clonotypes.  

To estimate clonal breadth of antigen-associated TCRs within 
a bulk repertoires with (N) productive clonotypes and (M) total 
productive templates, we use a set of (X) previously defined
antigen-associated meta-clonotypes (defined as a (i) Centroid TRV,CDR3, 
(ii) TCR-specific RADIUS, (iii) MOTIF. 

1. Compute the TCRdist between each centroid TCRij for i {1...i...X} 
and all bulk clones {1...j..N} using rectangular search with the 
tcrdist.repertoires.TCRrep.compute_sparse_rect_distances(), producing 
a sparse distance matrix. 

2. Next perform a long-form tabulation that records the network formed between 
all meta-clonotype centroids and bulk sequences within the specified radius. 
This is performed with the function tcrdist.breadth.long_form_tabulation().

The network is represented as a Pandas DataFrame. Where centroid sequences
are recorded as "cdr3_b_aa", "v_b_gene", "j_b_gene" and the conformant sequence
in the bulk repertoire is "cdr3_b_aa_hit", 'v_b_gene_hit', 'j_b_gene_hit'. 
Crucially there is a column "MOTIF" which indicates whether the CDR3 of 
the hit sequence is conformant with the regular expression in the column 
"regex". 

3. The long-form Pandas DataFrame can then be used as input to the function 
tcrdist.breadth.estimate_breadth_and_depth(). The unit of analysis -- that is, 
whether breadth refers to pathogen, protein, or epitope specific breadth --
can be specified in with the argument <breadth_cols>. Crucially, when 
running tcrdist.breadth.long_form_tabulation() the argument 
<search_cols> must include a column indicating the association between 
a metac-clonotype and a particular 'protein' or 'epitope'  
e.g., ['tag', 'protein', 'epitope', cdr3_b_aa', 'v_b_gene', 'j_b_gene', 
'pgen','regex', 'radius']

Note that breadth and depth follow the definitions in Synder et al. 2020 
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418734/), 
see section: The breadth and depth of a disease-specific T-cell response. 

"""
import pytest

#Download many bulk repertoires or just 2
def test_breadth_estimation_workflow(testing_only = True):
    import numpy as np
    import os 
    import pandas as pd
    from tcrdist.setup_tests import download_and_extract_zip_file
    from tcrdist.repertoire import TCRrep
    from tcrdist.breadth import get_safe_chunk
    from tcrdist.breadth import long_form_tabulation
    from tcrdist.breadth import estimate_breadth_and_depth
    import multiprocessing
    import time
    import scipy.sparse
    # For testing use only a subset of meta-clonotypes
  
    ncpus = min(multiprocessing.cpu_count(), 6)
    # Download 2 bulk files for testing and demonstration purposes
    files = [
    '1588BW_20200417_PBMC_unsorted_cc1000000_ImmunRACE_050820_008_gDNA_TCRB.tsv.tcrdist3.tsv',
    '1349BW_unsorted_cc1000000_ImmunRACE_043020_003_gDNA_TCRB.tsv.tcrdist3.tsv']
    if not np.all([os.path.isfile(f) for f in files]):
        download_and_extract_zip_file(
            "ImmunoSeq_MIRA_matched_tcrdist3_ready_2_files.zip", 
            source = "dropbox", 
            dest = ".")
        assert np.all([os.path.isfile(f) for f in files])
    # Download a Meta-Clonotypes File, All Meta-Clonotypes from Pre-Print Manuscript 
    # (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7781332/)
    if not os.path.isfile("bioRxiv_v2_metaclonotypes.tsv.zip"): 
        download_and_extract_zip_file('bioRxiv_v2_metaclonotypes.tsv.zip', source = "dropbox", dest = ".")
        assert os.path.isfile('bioRxiv_v2_metaclonotypes.tsv')
    # ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
    # ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
    # ███████╗█████╗  ███████║██████╔╝██║     ███████║
    # ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
    # ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
    # ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
    meta_clonotypes_filename = 'bioRxiv_v2_metaclonotypes.tsv'
    meta_clonotypes          = pd.read_csv(meta_clonotypes_filename , sep = "\t")
    # Assign feature id as a tag for tracking by meta-clonotype definition
    meta_clonotypes['tag'] = meta_clonotypes['feature'].copy()
    cols = ['cdr3_b_aa','v_b_gene','j_b_gene',
            'pgen','radius','regex','protein','protein_coordinate', 'tag']
    df_search = meta_clonotypes[cols]
    df_search = df_search.assign(count= 1)

    if testing_only: 
        # take a subsample to accelerated unit testing
        df_search = df_search.sample(100, random_state=1)
    tr_search = TCRrep(cell_df = df_search,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
    tr_search.cpus = ncpus
    # ██████╗ ██╗   ██╗██╗     ██╗  ██╗
    # ██╔══██╗██║   ██║██║     ██║ ██╔╝
    # ██████╔╝██║   ██║██║     █████╔╝ 
    # ██╔══██╗██║   ██║██║     ██╔═██╗ 
    # ██████╔╝╚██████╔╝███████╗██║  ██╗
    # ╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝
    bulk_filename = files[0]
    df_bulk = pd.read_csv( bulk_filename, sep = "\t")
    valid_TRBV_CDR3_TRBJ = (df_bulk['valid_cdr3'] == True)&(df_bulk['v_b_gene'].notna())&(df_bulk['j_b_gene'].notna())
    df_bulk = df_bulk[valid_TRBV_CDR3_TRBJ].reset_index(drop = True)
    # Convert templates to counts
    if 'templates' in df_bulk.columns:
        df_bulk = df_bulk[['cdr3_b_aa','v_b_gene','j_b_gene','templates',
        'productive_frequency']].\
        rename(columns = {'templates':'count'})
    else:
        df_bulk = df_bulk[['cdr3_b_aa','v_b_gene','j_b_gene','count','productive_frequency']]
    tr_bulk   = TCRrep(cell_df = df_bulk,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
    # Get dimensions
    search_clones = tr_search.clone_df.shape[0] 
    bulk_clones   = tr_bulk.clone_df.shape[0]
    # Get and ideal chunksize that will control memory usage
    # 10**7 will keep memory > 2GB per CPU. 
    ideal_chunk_size = get_safe_chunk(
        search_clones = tr_search.clone_df.shape[0], 
        bulk_clones = tr_bulk.clone_df.shape[0],
        target = 10**7)
    print(f"IDEAL CHUNK SIZE {ideal_chunk_size},{ncpus} CPUS")
    # ██████╗ ███████╗ ██████╗████████╗    ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
    # ██╔══██╗██╔════╝██╔════╝╚══██╔══╝    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
    # ██████╔╝█████╗  ██║        ██║       ███████╗█████╗  ███████║██████╔╝██║     ███████║
    # ██╔══██╗██╔══╝  ██║        ██║       ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
    # ██║  ██║███████╗╚██████╗   ██║       ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
    # ╚═╝  ╚═╝╚══════╝ ╚═════╝   ╚═╝       ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
    tic = time.perf_counter()
    tr_search.compute_sparse_rect_distances(
        df = tr_search.clone_df, 
        df2 = tr_bulk.clone_df, 
        chunk_size = ideal_chunk_size, 
        radius = 37)
    toc = time.perf_counter()
    print(f"SEARCHING {search_clones} META-CLONOTYPES IN {bulk_clones} BULK CLONES IN {toc - tic:0.2f} sec.")
    # ████████╗ █████╗ ██████╗ ██╗   ██╗██╗      █████╗ ████████╗███████╗
    # ╚══██╔══╝██╔══██╗██╔══██╗██║   ██║██║     ██╔══██╗╚══██╔══╝██╔════╝
    #    ██║   ███████║██████╔╝██║   ██║██║     ███████║   ██║   █████╗  
    #    ██║   ██╔══██║██╔══██╗██║   ██║██║     ██╔══██║   ██║   ██╔══╝  
    #    ██║   ██║  ██║██████╔╝╚██████╔╝███████╗██║  ██║   ██║   ███████╗
    #    ╚═╝   ╚═╝  ╚═╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚══════╝
    df_tab = long_form_tabulation(
        clone_df1   = tr_search.clone_df,
        clone_df2   = tr_bulk.clone_df,
        csrmat      = tr_search.rw_beta, 
        search_cols = ['tag', 'protein','protein_coordinate','cdr3_b_aa', 
        'v_b_gene', 'j_b_gene', 'pgen','regex', 'radius'])
    # ██████╗ ██████╗ ███████╗ █████╗ ██████╗ ████████╗██╗  ██╗
    # ██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗╚══██╔══╝██║  ██║
    # ██████╔╝██████╔╝█████╗  ███████║██║  ██║   ██║   ███████║
    # ██╔══██╗██╔══██╗██╔══╝  ██╔══██║██║  ██║   ██║   ██╔══██║
    # ██████╔╝██║  ██║███████╗██║  ██║██████╔╝   ██║   ██║  ██║
    # ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═════╝    ╚═╝   ╚═╝  ╚═╝
    df_breadth = estimate_breadth_and_depth(df=df_tab, 
        breadth_cols = ['protein'], 
        N =tr_bulk.clone_df.shape[0] , 
        M= tr_bulk.clone_df['count'].sum(), 
        motif = True,
        exact = False)
    df_breadth = df_breadth.assign(file = files[0])
    # ███╗   ███╗ █████╗ ███╗   ██╗██╗   ██╗    ███████╗██╗██╗     ███████╗███████╗
    # ████╗ ████║██╔══██╗████╗  ██║╚██╗ ██╔╝    ██╔════╝██║██║     ██╔════╝██╔════╝
    # ██╔████╔██║███████║██╔██╗ ██║ ╚████╔╝     █████╗  ██║██║     █████╗  ███████╗
    # ██║╚██╔╝██║██╔══██║██║╚██╗██║  ╚██╔╝      ██╔══╝  ██║██║     ██╔══╝  ╚════██║
    # ██║ ╚═╝ ██║██║  ██║██║ ╚████║   ██║       ██║     ██║███████╗███████╗███████║
    # ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝   ╚═╝       ╚═╝     ╚═╝╚══════╝╚══════╝╚══════╝
    import numpy as np
    import os 
    import pandas as pd
    from tcrdist.setup_tests import download_and_extract_zip_file
    from tcrdist.repertoire import TCRrep
    from tcrdist.breadth import get_safe_chunk
    from tcrdist.breadth import long_form_tabulation
    from tcrdist.breadth import estimate_breadth_and_depth
    import multiprocessing
    import time
    import scipy.sparse
    ncpus = min(multiprocessing.cpu_count(), 6)
    destination = "breadth_estimates"
    if not os.path.isdir(destination):
        os.mkdir(destination)
    files = [
    '1588BW_20200417_PBMC_unsorted_cc1000000_ImmunRACE_050820_008_gDNA_TCRB.tsv.tcrdist3.tsv',
    '1349BW_unsorted_cc1000000_ImmunRACE_043020_003_gDNA_TCRB.tsv.tcrdist3.tsv']

    # ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
    # ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
    # ███████╗█████╗  ███████║██████╔╝██║     ███████║
    # ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
    # ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
    # ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
    meta_clonotypes_filename = 'bioRxiv_v2_metaclonotypes.tsv'
    meta_clonotypes          = pd.read_csv(meta_clonotypes_filename , sep = "\t")
    # Assign feature id as a tag for tracking by meta-clonotype definition
    meta_clonotypes['tag'] = meta_clonotypes['feature'].copy()
    cols = ['cdr3_b_aa','v_b_gene','j_b_gene',
            'pgen','radius','regex','protein','protein_coordinate', 'tag']
    df_search = meta_clonotypes[cols]
    df_search = df_search.assign(count= 1)
    if testing_only: 
        # take a subsample to accelerated unit testing
        df_search = df_search.sample(100, random_state=1)
    tr_search = TCRrep(cell_df = df_search,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
    tr_search.cpus = ncpus

    for file in files:
        bulk_filename = file
        df_bulk = pd.read_csv( bulk_filename, sep = "\t")
        valid_TRBV_CDR3_TRBJ = (df_bulk['valid_cdr3'] == True)&(df_bulk['v_b_gene'].notna())&(df_bulk['j_b_gene'].notna())
        df_bulk = df_bulk[valid_TRBV_CDR3_TRBJ].reset_index(drop = True)
        # Convert templates to counts
        if 'templates' in df_bulk.columns:
            df_bulk = df_bulk[['cdr3_b_aa','v_b_gene','j_b_gene','templates','productive_frequency']].rename(columns = {'templates':'count'})
        else:
            df_bulk = df_bulk[['cdr3_b_aa','v_b_gene','j_b_gene','count','productive_frequency']]
        tr_bulk   = TCRrep(cell_df = df_bulk,
                        organism = 'human',
                        chains = ['beta'],
                        db_file = 'alphabeta_gammadelta_db.tsv',
                        compute_distances = False)
        # Get dimensions
        search_clones = tr_search.clone_df.shape[0] 
        bulk_clones   = tr_bulk.clone_df.shape[0]
        # Get and ideal chunksize that will control memory usage
        # 10**7 will keep memory > 2GB per CPU. 
        ideal_chunk_size = get_safe_chunk(tr_search.clone_df.shape[0], tr_bulk.clone_df.shape[0], target = 10**7)
        print(f"IDEAL CHUNK SIZE {ideal_chunk_size},{ncpus} CPUS")
        # ██████╗ ███████╗ ██████╗████████╗    ███████╗███████╗ █████╗ ██████╗  ██████╗██╗  ██╗
        # ██╔══██╗██╔════╝██╔════╝╚══██╔══╝    ██╔════╝██╔════╝██╔══██╗██╔══██╗██╔════╝██║  ██║
        # ██████╔╝█████╗  ██║        ██║       ███████╗█████╗  ███████║██████╔╝██║     ███████║
        # ██╔══██╗██╔══╝  ██║        ██║       ╚════██║██╔══╝  ██╔══██║██╔══██╗██║     ██╔══██║
        # ██║  ██║███████╗╚██████╗   ██║       ███████║███████╗██║  ██║██║  ██║╚██████╗██║  ██║
        # ╚═╝  ╚═╝╚══════╝ ╚═════╝   ╚═╝       ╚══════╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
        tic = time.perf_counter()
        tr_search.compute_sparse_rect_distances(
            df = tr_search.clone_df, 
            df2 = tr_bulk.clone_df, 
            chunk_size = ideal_chunk_size, 
            radius = 37)
        toc = time.perf_counter()
        print(f"SEARCHING {search_clones} META-CLONOTYPES IN {bulk_clones} BULK CLONES IN {toc - tic:0.2f} sec.")
        # ████████╗ █████╗ ██████╗ ██╗   ██╗██╗      █████╗ ████████╗███████╗
        # ╚══██╔══╝██╔══██╗██╔══██╗██║   ██║██║     ██╔══██╗╚══██╔══╝██╔════╝
        #    ██║   ███████║██████╔╝██║   ██║██║     ███████║   ██║   █████╗  
        #    ██║   ██╔══██║██╔══██╗██║   ██║██║     ██╔══██║   ██║   ██╔══╝  
        #    ██║   ██║  ██║██████╔╝╚██████╔╝███████╗██║  ██║   ██║   ███████╗
        #    ╚═╝   ╚═╝  ╚═╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚══════╝
        df_tab = long_form_tabulation(
            clone_df1   = tr_search.clone_df,
            clone_df2   = tr_bulk.clone_df,
            csrmat      = tr_search.rw_beta, 
            search_cols = ['tag', 'protein','protein_coordinate','cdr3_b_aa', 
            'v_b_gene', 'j_b_gene', 'pgen','regex', 'radius'])
        
        # ██████╗ ██████╗ ███████╗ █████╗ ██████╗ ████████╗██╗  ██╗
        # ██╔══██╗██╔══██╗██╔════╝██╔══██╗██╔══██╗╚══██╔══╝██║  ██║
        # ██████╔╝██████╔╝█████╗  ███████║██║  ██║   ██║   ███████║
        # ██╔══██╗██╔══██╗██╔══╝  ██╔══██║██║  ██║   ██║   ██╔══██║
        # ██████╔╝██║  ██║███████╗██║  ██║██████╔╝   ██║   ██║  ██║
        # ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═════╝    ╚═╝   ╚═╝  ╚═╝
        df_breadth = estimate_breadth_and_depth(df=df_tab, 
            breadth_cols = ['protein'], 
            N =tr_bulk.clone_df.shape[0] , 
            M= tr_bulk.clone_df['count'].sum(), 
            motif = True,
            exact = False)
        df_breadth = df_breadth.assign(file = files[0])
        # save rectangular sparse matrix, tabulation, and breadth esimates for future reference
        scipy.sparse.save_npz(os.path.join(destination, f"{file}.rw.npz"), tr_search.rw_beta)
        df_tab.to_csv(os.path.join(destination, f"{file}.tabulation.tsv"), 
            sep = '\t', index = True)
        df_breadth.to_csv(os.path.join(destination, f"{file}.breadth.tsv"),
            sep = '\t', index = True)
