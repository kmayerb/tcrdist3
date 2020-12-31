import os
import pandas as pd 
from tcrdist.adpt_funcs import import_adaptive_file
from tcrdist.repertoire import TCRrep
from tcrdist.tabulate import tabulate
import math
import time
"""
You can download the file with wget or follow the link
!wget https://www.dropbox.com/s/1zbcf1ep8kgmlzy/KHBR20-00107_TCRB.tsv
"""
bulk_file = 'KHBR20-00107_TCRB.tsv'
df_bulk = import_adaptive_file(bulk_file)
df_bulk = df_bulk[['cdr3_b_aa',
                    'v_b_gene',
                    'j_b_gene',
                    'templates',
                    'productive_frequency',
                    'valid_cdr3']].\
        rename(columns = {'templates':'count'})

df_bulk = df_bulk[(df_bulk['v_b_gene'].notna()) & (df_bulk['j_b_gene'].notna())].reset_index()
tr_bulk = TCRrep(cell_df = df_bulk,
                 organism = 'human',
                 chains = ['beta'],
                 db_file = 'alphabeta_gammadelta_db.tsv',
                 compute_distances = False)

search_file = os.path.join(
    'tcrdist',
    'data',
    'covid19',
    'mira_epitope_48_610_YLQPRTFL_YLQPRTFLL_YYVGYLQPRTF.tcrdist3.csv.ranked_centers_bkgd_ctlr_1E6.tsv')

df_search = pd.read_csv(search_file, sep = '\t')
df_search = df_search[['cdr3_b_aa','v_b_gene','j_b_gene','pgen','radius','regex']]
tr_search = TCRrep(cell_df = df_search,
                   organism = 'human',
                   chains = ['beta'],
                   db_file = 'alphabeta_gammadelta_db.tsv',
                   compute_distances = False)
tr_search.cpus = 1
tic = time.perf_counter()
tr_search.compute_sparse_rect_distances(df = tr_search.clone_df, df2 = tr_bulk.clone_df, chunk_size = 50, radius = 50) 
results = tabulate(clone_df1 = tr_search.clone_df, clone_df2 = tr_bulk.clone_df, pwmat = tr_search.rw_beta)
toc = time.perf_counter()
print(f"TABULATED IN {toc - tic:0.4f} seconds")
"""
Results tell us about abundance of each meta-clonotype feature in the bulk data.
'cdr3_b_aa' - meta-clonotype centroid CDR3
'v_b_gene'  - meta-clonotype centroid TRBV
'j_b_gene'  - meta-clonotype centroid TRBJ
'pgen'.     - meta-clonotype centroid TRBV,CDR3 pgen
'radius'    - meta-clonotype maximum neighbor TCRdist
'regex'     - pattern of conserved positions learned from the Antigen Enriched Dataset
'cdr1_b_aa' - meta-clonotype centroid CDR1
'cdr2_b_aa' - meta-clonotype centroid CDR2
'pmhc_b_aa' - meta-clonotype centroid CDR2.5
'bulk_sum_freq'    - In the bulk sample, frequency of TCRs within <radius> of centroid (RADIUS)
'bulk_sum_counts'  - In the bulk sample, total template TCRs within <radius> of centroid (RADIUS)
'bulk_seqs'        - In the bulk sample, CDR3 seqs ot TCRs within <radius> of centroid
'bulk_v_genes'     - In the bulk sample, TRBVs of TCRs within <radius> of centroid
'bulk_j_genes'     - In the bulk sample, TRBJs of TCRs within <radius> of centroid
'bulk_distances'   - In the bulk sample, distances of TCRs from centroid, for those  within <radius> of centroid
'bulk_counts'      - In the bulk sample, individual counts of each TCR within radius of centroid
'bulk_freqs'       - In the bulk sample, individual frequencies of each TCR within radius of centroid
'bulk_regex_match' - In the bulk sample, did CDR3 match regex motif pattern
'bulk_sum_freqs_regex_adj'  - In the bulk sample, sum of frequency of TCRs within <radius> of centroid and matching regex (RADIUS + MOTIF)
'bulk_sum_counts_regex_adj' - In the bulk sample, sum of counts of TCRs within <radius> of centroid and matching regex (RADIUS + MOTIF)
'bulk_sum_freqs_tcrdist0'   - In the bulk sample, sum of frequencies of TCRs within <radius> = 0 of centroid (EXACT)
'bulk_sum_counts_tcrdist0'  - In the bulk sample, sum of counts of TCRs within <radius> = 0 of centroid (EXACT)
"""



