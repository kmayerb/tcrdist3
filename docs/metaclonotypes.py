"""
2020-12-18
Seattle, WA
HLA CLASS I restricted meta-clonotype discovery.

This functions encapsulates a complete work flow for finding meta-clonotypes 
in antigen-enriched data.

prepare docker container:
apt-get install unzip
python3 -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
"""
import sys
import os
import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
from tcrdist.paths import path_to_base
from tcrdist.repertoire import TCRrep
from tcrdist.automate import auto_pgen
from tcrsampler.sampler import TCRsampler
from tcrdist.background import get_stratified_gene_usage_frequency
from tcrdist.background import sample_britanova
from tcrdist.background import make_gene_usage_counter, get_gene_frequencies
from tcrdist.background import calculate_adjustment, make_gene_usage_counter
from tcrdist.background import make_vj_matched_background, make_flat_vj_background
from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
from tcrdist.centers import calc_radii
from tcrdist.public import _neighbors_variable_radius
from tcrdist.public import _neighbors_sparse_variable_radius
from tcrdist.regex import _index_to_regex_str

def find_metaclonotypes(
    project_path = "tutorial48",
    source_path = os.path.join(path_to_base,'tcrdist','data','covid19'),
    antigen_enriched_file = 'mira_epitope_48_610_YLQPRTFL_YLQPRTFLL_YYVGYLQPRTF.tcrdist3.csv',
    ncpus = 4, 
    seed = 3434):
    """
    This functions encapsulates a complete 
    workflow for finding meta-clonotypes in antigen-enriched data.
    """
    np.random.seed(seed)
    if not os.path.isdir(project_path):
        os.mkdir(project_path)
    ############################################################################
    # Step 1: Select and load a antigen-enriched (sub)repertoire.           ####
    ############################################################################
    print(f"INITIATING A TCRrep() with {antigen_enriched_file}")
    assert os.path.isfile(os.path.join(source_path, antigen_enriched_file))
        # Read file into a Pandas DataFrame <df>
    df = pd.read_csv(os.path.join(source_path, antigen_enriched_file))
        # Drop cells without any gene usage information
    df = df[( df['v_b_gene'].notna() ) & (df['j_b_gene'].notna()) ]
        # Initialize a TCRrep class, using ONLY columns that are complete and unique define a a clone.
        # Class provides a 'count' column if non is present
        # Counts of identical subject:VCDR3 'clones' will be aggregated into a TCRrep.clone_df.
    from tcrdist.repertoire import TCRrep
    tr = TCRrep(cell_df = df[['subject','cell_type','v_b_gene', 'j_b_gene', 'cdr3_b_aa']], 
                organism = "human", 
                chains = ['beta'], 
                compute_distances = True)
    tr.cpus = ncpus
    ############################################################################
    # Step 1.1: Estimate Probability of Generation                          ####
    ############################################################################
    ### It will be useful later to know the pgen of each
    from tcrdist.automate import auto_pgen
    print(f"COMPUTING PGEN WITH OLGA (Sethna et al 2018)")
    print("FOR ANTIGEN-ENRICHED CLONES TO BE USED FOR SUBSEQUENT ANALYSES")
    auto_pgen(tr)

    # Tip: Users of tcrdist3 should be aware that by default a <TCRrep.clone_df> 
    # DataFrame is created out of non-redundant cells in the cell_df, and 
    # pairwise distance matrices automatically computed.
    # Notice that attributes <tr.clone_df>  and  <tr.pw_beta> , <tr.pw_cdr3_b_aa>, 
    # are immediately accessible.
    # Attributes <tr.pw_pmhc_b_aa>, <tr.pw_cdr2_b_aa>, and <tr.pw_cdr1_b_aa>  
    # are also available if <TCRrep.store_all_cdr> is set to True.
    # For large datasets, i.e., >15,000 clones, this approach may consume too much 
    # memory so <TCRrep.compute_distances> is automatically set to False. 
                                    
    ############################################################################
    # Step 2: Synthesize an Inverse Probability Weighted VJ Matched Background #
    ############################################################################
    # Generating an appropriate set of unenriched reference TCRs is important; for
    # each set of antigen-associated TCRs, discovered by MIRA, we created a two part
    # background. One part consists of 100,000 synthetic TCRs whose V-gene and J-gene
    # frequencies match those in the antigen-enriched repertoire, using the software
    # OLGA (Sethna et al. 2019; Marcou et al. 2018). The other part consists of
    # 100,000 umbilical cord blood TCRs sampled uniformly from 8 subjects (Britanova
    # et al., 2017). This mix balances dense sampling of sequences near the
    # biochemical neighborhoods of interest with broad sampling of TCRs from an
    # antigen-naive repertoire. Importantly, we adjust for the biased sampling by
    # using the V- and J-gene frequencies observed in the cord-blood data (see
    # Methods for details about inverse probability weighting adjustment). Using this
    # approach we are able to estimate the abundance of TCRs similar to a centroid
    # TCR in an unenriched background repertoire of ~1,000,000 TCRs, using a
    # comparatively modest background dataset of 200,000 TCRs. While this estimate
    # may underestimate the true specificity, since some of the neighborhood TCRs in
    # the unenriched background repertoire may in fact recognize the antigen of
    # interest, it is useful for prioritizing neighborhoods and selecting a radius
    # for each neighborhood that balances sensitivity and specificity.
    # Initialize a TCRsampler -- human, beta, umbilical cord blood from 8 people.
    print(f"USING tcrsampler TO CONSTRUCT A CUSTOM V-J MATCHED BACKGROUND")
    from tcrsampler.sampler import TCRsampler
    ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
    # Stratify sample so that each subject contributes similarly to estimate of 
    # gene usage frequency
    from tcrdist.background import get_stratified_gene_usage_frequency
    ts = get_stratified_gene_usage_frequency(ts = ts, replace = True) 
    # Synthesize an inverse probability weighted V,J gene background that matches 
    # usage in your enriched repertoire 
    df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = 'beta')
    # Get a randomly drawn stratified sampler of beta, cord blood from 
    # Britanova et al. 2016 
    # Dynamics of Individual T Cell Repertoires: From Cord Blood to Centenarians
    from tcrdist.background import  sample_britanova
    df_britanova_100K = sample_britanova(size = 100000)
    # Append frequency columns using, using sampler above
    df_britanova_100K = get_gene_frequencies(ts = ts, df = df_britanova_100K)
    df_britanova_100K['weights'] = 1
    df_britanova_100K['source'] = "stratified_random"
    # Combine the two parts of the background into a single DataFrame
    df_bkgd = pd.concat([df_vj_background.copy(), df_britanova_100K.copy()], axis = 0).\
        reset_index(drop = True)                                              
    # Assert that the backgrounds have the expected number of rows.
    assert df_bkgd.shape[0] == 200000
    # Save the background for future use
    background_outfile = os.path.join(project_path, f"{antigen_enriched_file}.olga100K_brit100K_bkgd.csv")
    print(f'WRITING {background_outfile}')
    df_bkgd.to_csv(background_outfile, index = False)
    # Load the background to a TCRrep without computing pairwise distances 
    # (i.e., compute_distances = False)
    tr_bkgd = TCRrep(
        cell_df = df_bkgd,
        organism = "human", 
        chains = ['beta'], 
        compute_distances = False)
    # Compute rectangular distances. Those are, distances between each clone in 
    # the antigen-enriched repertoire and each TCR in the background.
    # With a single 1 CPU and < 10GB RAM, 5E2x2E5 = 100 million pairwise distances, 
    # across CDR1, CDR2, CDR2.5, and CDR3 
    # 1min 34s ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each) 
    # %timeit -r 1 tr.compute_rect_distances(df = tr.clone_df, df2 = tr_bkdg.clone_df, store = False)
    ############################################################################
    # Step 4: Calculate Distances                                          #####
    ############################################################################
    print(f"COMPUTING RECTANGULARE DISTANCE")
    tr.compute_sparse_rect_distances(
        df = tr.clone_df, 
        df2 = tr_bkgd.clone_df,
        radius=50,
        chunk_size = 100)
    scipy.sparse.save_npz(os.path.join(project_path, f"{antigen_enriched_file}.rw_beta.npz"), tr.rw_beta)
        # Tip: For larger dataset you can use a sparse implementation: 
        # 30.8 s ± 0 ns per loop ; tr.cpus = 6
        # %timeit -r tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = tr_bkdg.clone_df,radius=50, chunk_size=85)
    ############################################################################
    # Step 5: Examine Density ECDFS                                        #####
    ############################################################################
        # Investigate the density of neighbors to each TCR, based on expanding 
        # distance radius.
    from tcrdist.ecdf import distance_ecdf, _plot_manuscript_ecdfs
    import matplotlib.pyplot as plt
        # Compute empirical cumulative density function (ecdf)
        # Compare Antigen Enriched TCRs (against itself).
    thresholds, antigen_enriched_ecdf = distance_ecdf(
        tr.pw_beta,
        thresholds=range(0,50,2))
        # Compute empirical cumulative density function (ecdf)
        # Compare Antigen Enriched TCRs (against) 200K probability 
        # inverse weighted background
    thresholds, background_ecdf = distance_ecdf(
        tr.rw_beta,
        thresholds=range(0,50,2),
        weights= tr_bkgd.clone_df['weights'], 
        absolute_weight = True)
        # plot_ecdf similar to tcrdist3 manuscript #
    antigen_enriched_ecdf[antigen_enriched_ecdf == antigen_enriched_ecdf.min()] = 1E-10
    f1 = _plot_manuscript_ecdfs(
        thresholds, 
        antigen_enriched_ecdf, 
        ylab= 'Proportion of Antigen Enriched TCRs', 
        cdr3_len=tr.clone_df.cdr3_b_aa.str.len(), 
        min_freq=1E-10)
    f1.savefig(os.path.join(project_path, f'{antigen_enriched_file}.ecdf_AER_plot.png'))
    f2 = _plot_manuscript_ecdfs(
        thresholds,
        background_ecdf,
        ylab= 'Proportion of Reference TCRs',
        cdr3_len=tr.clone_df.cdr3_b_aa.str.len(),
        min_freq=1E-10)
    f2.savefig(os.path.join(project_path, f'{antigen_enriched_file}.ecdf_BUR_plot.png'))
    ############################################################################
    # Step 6: Find optimal radii  (theta = 1E5                             #####
    ############################################################################
    # To ascertain which meta-clonotypes are likely to be most specific, 
    # take advantage of an existing function <bkgd_cntrl_nn2>.                                                                                                                                  
    #  d888   .d8888b.  8888888888     888888888  
    # d8888  d88P  Y88b 888            888        
    #   888  888    888 888            888        
    #   888  888    888 8888888        8888888b.  
    #   888  888    888 888                 "Y88b 
    #   888  888    888 888      888888       888 
    #   888  Y88b  d88P 888            Y88b  d88P 
    # 8888888 "Y8888P"  8888888888      "Y8888P"                                         
   
    level_tag = '1E5'
    from tcrdist.neighbors import bkgd_cntl_nn2
    centers_df  = bkgd_cntl_nn2(
        tr               = tr,
        tr_background    = tr_bkgd,
        weights          = tr_bkgd.clone_df.weights,
        ctrl_bkgd        = 10**-5, 
        col              = 'cdr3_b_aa',
        add_cols         = ['v_b_gene', 'j_b_gene'],
        ncpus            = 4,
        include_seq_info = True,
        thresholds       = [x for x in range(0,50,2)],
        generate_regex   = True,
        test_regex       = True,
        forced_max_radius = 36)

    ############################################################################
    # Step 6.2: (theta = 1E5) ALL meta-clonotypes .tsv file                   ##
    ############################################################################
    # save center to project_path for future use
    centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )
    
    # Many of meta-clonotypes contain redundant information. 
    # We can winnow down to less-redundant list. We do this 
    # by ranking clonotypes from most to least specific. 
        # <min_nsubject> is minimum publicity of the meta-clonotype,  
        # <min_nr> is minimum non-redundancy
    # Add neighbors, K_neighbors, and nsubject columns
    from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
    centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, radius_list = centers_df['radius'])
    centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
    # We determine how many <nsubjects> are in the set of neighbors 
    centers_df['nsubject']  = centers_df['neighbors'].\
            apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
    centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )

    from tcrdist.centers import rank_centers
    ranked_centers_df = rank_centers(
        centers_df = centers_df, 
        rank_column = 'chi2joint', 
        min_nsubject = 2, 
        min_nr = 1)
    ############################################################################
    # Step 6.3:  (theta = 1E5) NR meta-clonotypes .tsv file                  ###
    ############################################################################
    # Output, ready to search bulk data.
    ranked_centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )
    ############################################################################
    # Step 6.4: (theta = 1E5) Output Meta-Clonotypes HTML Summary            ###
    ############################################################################
    # Here we can make a svg logo for each NR meta-clonotype
    if ranked_centers_df.shape[0] > 0:
        from progress.bar import IncrementalBar
        from tcrdist.public import make_motif_logo
        cdr3_name = 'cdr3_b_aa'
        v_gene_name = 'v_b_gene'
        svgs = list()
        svgs_raw = list()
        bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        for i,r in ranked_centers_df.iterrows():
            bar.next()
            centroid = r[cdr3_name]
            v_gene   = r[v_gene_name]
            svg, svg_raw = make_motif_logo( tcrsampler = ts, 
                                            pwmat = tr.pw_beta,
                                            clone_df = tr.clone_df,
                                            centroid = centroid ,
                                            v_gene = v_gene ,
                                            radius = r['radius'],
                                            pwmat_str = 'pw_beta',
                                            cdr3_name = 'cdr3_b_aa',
                                            v_name = 'v_b_gene',
                                            gene_names = ['v_b_gene','j_b_gene'])
            svgs.append(svg)
            svgs_raw.append(svg_raw)
        bar.next();bar.finish()
        ranked_centers_df['svg']      = svgs
        ranked_centers_df['svg_raw'] = svgs_raw

        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        labels =['cdr3_b_aa','v_b_gene', 'j_b_gene', 'pgen',
                'radius', 'regex','nsubject','K_neighbors', 
                'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.html')
        # 888    888 88888888888 888b     d888 888      
        # 888    888     888     8888b   d8888 888      
        # 888    888     888     88888b.d88888 888      
        # 8888888888     888     888Y88888P888 888      
        # 888    888     888     888 Y888P 888 888      
        # 888    888     888     888  Y8P  888 888      
        # 888    888     888     888   "   888 888      
        # 888    888     888     888       888 88888888
        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")
    # To ascertain which meta-clonotypes are likely to be most specific, 
    # take advantage of an existing function <bkgd_cntrl_nn2>.       
    #  d888   .d8888b.  8888888888       .d8888b.  
    # d8888  d88P  Y88b 888             d88P  Y88b 
    #   888  888    888 888             888        
    #   888  888    888 8888888         888d888b.  
    #   888  888    888 888             888P "Y88b 
    #   888  888    888 888      888888 888    888 
    #   888  Y88b  d88P 888             Y88b  d88P 
    # 8888888 "Y8888P"  8888888888       "Y8888P" 
    ############################################################################
    # Step 6.5: Find optimal radii  (theta = 1E6)                            ###
    ############################################################################
    level_tag = '1E6'
    from tcrdist.neighbors import bkgd_cntl_nn2
    centers_df  = bkgd_cntl_nn2(
        tr               = tr,
        tr_background    = tr_bkgd,
        weights          = tr_bkgd.clone_df.weights,
        ctrl_bkgd        = 10**-6, 
        col              = 'cdr3_b_aa',
        add_cols         = ['v_b_gene', 'j_b_gene'],
        ncpus            = 4,
        include_seq_info = True,
        thresholds       = [x for x in range(0,50,2)],
        generate_regex   = True,
        test_regex       = True,
        forced_max_radius = 36)
    ############################################################################
    # Step 6.6: (theta = 1E6) ALL meta-clonotypes .tsv file                   ##
    ############################################################################
    # save center to project_path for future use
    centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )
    
    # Many of meta-clonotypes contain redundant information. 
    # We can winnow down to less-redundant list. We do this 
    # by ranking clonotypes from most to least specific. 
        # <min_nsubject> is minimum publicity of the meta-clonotype,  
        # <min_nr> is minimum non-redundancy
    # Add neighbors, K_neighbors, and nsubject columns
    from tcrdist.public import _neighbors_variable_radius, _neighbors_sparse_variable_radius
    centers_df['neighbors'] = _neighbors_variable_radius(pwmat=tr.pw_beta, radius_list = centers_df['radius'])
    centers_df['K_neighbors'] = centers_df['neighbors'].apply(lambda x : len(x))
    # We determine how many <nsubjects> are in the set of neighbors 
    centers_df['nsubject']  = centers_df['neighbors'].\
            apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
    centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )

    from tcrdist.centers import rank_centers
    ranked_centers_df = rank_centers(
        centers_df = centers_df, 
        rank_column = 'chi2joint', 
        min_nsubject = 2, 
        min_nr = 1)
    ############################################################################
    # Step 6.7:  (theta = 1E6) NR meta-clonotypes .tsv file                  ###
    ############################################################################
    # Output, ready to search bulk data.
    ranked_centers_df.to_csv( os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.tsv'), sep = "\t" )

    ############################################################################
    # Step 6.8: (theta = 1E6) Output Meta-Clonotypes HTML Summary            ###
    ############################################################################
    # Here we can make a svg logo for each meta-clonotype
    from progress.bar import IncrementalBar
    from tcrdist.public import make_motif_logo
    if ranked_centers_df.shape[0] > 0:
        cdr3_name = 'cdr3_b_aa'
        v_gene_name = 'v_b_gene'
        svgs = list()
        svgs_raw = list()
        bar = IncrementalBar('Processing', max = ranked_centers_df.shape[0])
        for i,r in ranked_centers_df.iterrows():
            bar.next()
            centroid = r[cdr3_name]
            v_gene   = r[v_gene_name]
            svg, svg_raw = make_motif_logo( tcrsampler = ts, 
                                            pwmat = tr.pw_beta,
                                            clone_df = tr.clone_df,
                                            centroid = centroid ,
                                            v_gene = v_gene ,
                                            radius = r['radius'],
                                            pwmat_str = 'pw_beta',
                                            cdr3_name = 'cdr3_b_aa',
                                            v_name = 'v_b_gene',
                                            gene_names = ['v_b_gene','j_b_gene'])
            svgs.append(svg)
            svgs_raw.append(svg_raw)
        bar.next();bar.finish()
        ranked_centers_df['svg']      = svgs
        ranked_centers_df['svg_raw'] = svgs_raw

        def shrink(s):
            return s.replace('height="100%"', 'height="20%"').replace('width="100%"', 'width="20%"')
        labels =['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'pgen', 'radius', 'regex','nsubject','K_neighbors', 'bkgd_hits_weighted','chi2dist','chi2re','chi2joint']
        
        output_html_name = os.path.join(project_path, f'{antigen_enriched_file}.ranked_centers_bkgd_ctlr_{level_tag}.html')
        # 888    888 88888888888 888b     d888 888      
        # 888    888     888     8888b   d8888 888      
        # 888    888     888     88888b.d88888 888      
        # 8888888888     888     888Y88888P888 888      
        # 888    888     888     888 Y888P 888 888      
        # 888    888     888     888  Y8P  888 888      
        # 888    888     888     888   "   888 888      
        # 888    888     888     888       888 88888888     
        with open(output_html_name, 'w') as output_handle:
            for i,r in ranked_centers_df.iterrows():
                #import pdb; pdb.set_trace()
                svg, svg_raw = r['svg'],r['svg_raw']
                output_handle.write("<br></br>")
                output_handle.write(shrink(svg))
                output_handle.write(shrink(svg_raw))
                output_handle.write("<br></br>")
                output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
                output_handle.write("<br></br>")




if __name__ == "__main__":
    # 888      888                               888             
    # 888      888                               888             
    # 888      888                               888             
    # 88888b.  888  8888b.     .d8888b   .d88b.  888888 .d8888b  
    # 888 "88b 888     "88b    88K      d8P  Y8b 888    88K      
    # 888  888 888 .d888888    "Y8888b. 88888888 888    "Y8888b. 
    # 888  888 888 888  888         X88 Y8b.     Y88b.       X88 
    # 888  888 888 "Y888888     88888P'  "Y8888   "Y888  88888P' 
    # December 2020 
    # Seattle, WA
    
    # If run as a script, this will find meta-clonotypes for the 18
    # SARS-CoV-2 MIRA sets with the most prior evidence of restriction
    # by as specific HLA-alleles.

    from collections import defaultdict 
    import os 
    import pandas as pd
    import re
    # <path> where files reside
    path = os.path.join(path_to_base,'tcrdist','data','covid19')
    # <all_files> list of all files
    all_files = [f for f in os.listdir(path) if f.endswith('.tcrdist3.csv')]
    # <restriction> list of tuples from Supporting Table S5
    restriction = \
    [('m_55_524_ALR','A*01'),
    ('m_1_8260_HTT','A*01'),
    ('m_45_689_SYF','A*01'),
    ('m_10_2274_LSP','B*07'),
    ('m_155_59_RAR','B*07'),
    ('m_133_102_NQK','B*15'),
    ('m_48_610_YLQ','A*02'),
    ('m_111_146_AEI','A*11'),
    ('m_53_532_NYL','A*24'),
    ('m_90_216_GYQ','C*07'),
    ('m_140_92_NSS','A*01'),
    ('m_55_524_ALR','B*08'),
    ('m_183_39_RIR','A*03'),
    ('m_10_2274_LSP','C*07'),
    ('m_99_191_QEC','B*40'),
    ('m_155_59_RAR','C*07'),
    ('m_185_39_ASQ','B*15'),
    ('m_147_73_DLF','B*08'),
    ('m_110_148_ELI','B*44'),
    ('m_51_546_AYK','A*03'),
    ('m_44_697_FPP','B*35'),
    ('m_118_136_ALN','A*11'),
    ('m_176_44_SST','A*11'),
    ('m_30_1064_KAY','B*57'),
    ('m_192_31_FQP','B*15'),
    ('m_70_345_DTD','A*01')]
    # <restrictions_dict> convert list to dictionary
    restrictions_dict = defaultdict(list)
    for k,v in restriction:
        restrictions_dict[k].append(v)
    # Loop through all files to construct a Dataframe of only those with strongest evidence of HLA-restriction
    cache = list()
    for f in all_files:
        rgs = re.search(pattern ='(mira_epitope)_([0-9]{1,3})_([0-9]{1,6})_([A-Z]{1,3})', 
            string = f).groups()    
        rgs4 = re.search(pattern ='(mira_epitope)_([0-9]{1,3})_([0-9]{1,6})_([A-Z]{1,4})', 
            string = f).groups()
        key3 = f'm_{rgs[1]}_{rgs[2]}_{rgs[3]}'
        key4 = f'm_{rgs4[1]}_{rgs4[2]}_{rgs4[3]}'
        setkey = f"M_{rgs[1]}"
        include = key3 in restrictions_dict.keys()
        alleles = restrictions_dict.get(key3)
        cache.append((setkey, key3, key4, f, int(rgs[1]), int(rgs[2]), include, alleles))
    # <hla_df>
    hla_df = pd.DataFrame(cache, columns = ['set', 'key3','key4','filename','set_rank','clones','hla_restricted','alleles']).\
        sort_values(['hla_restricted','clones'], ascending = True).\
        query('hla_restricted == True').reset_index(drop = True)
    
    from tcrdist.paths import path_to_base
    for ind,row in hla_df.iterrows():
        file = row['filename']
        print(file, ind)
        print(row)
        find_metaclonotypes(project_path = "hla_restricted_meta_clonotypes",
                            source_path = os.path.join(path_to_base,'tcrdist','data','covid19'),
                            antigen_enriched_file = file,
                            ncpus = 6, 
                            seed = 3434)


# 8888888888 888b    888 8888888b.  
# 888        8888b   888 888  "Y88b 
# 888        88888b  888 888    888 
# 8888888    888Y88b 888 888    888 
# 888        888 Y88b888 888    888 
# 888        888  Y88888 888    888 
# 888        888   Y8888 888  .d88P 
# 8888888888 888    Y888 8888888P"  

