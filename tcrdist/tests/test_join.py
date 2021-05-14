import pytest
import pandas as pd
v20df = pd.read_csv('https://www.dropbox.com/s/wmc5wc752t782kq/MIRA_v21_covid_diagnosed_sars_cov2_ci_epitope_specific_tcrs.tsv?dl=1', sep = "\t").head(1000)
v21df = pd.read_csv('https://www.dropbox.com/s/c3gfq1lu0xdefpy/MIRA_v20_covid_diagnosed_sars_cov2_ci_epitope_specific_tcrs.tsv?dl=1', sep = "\t").head(1000)

"""
Example 1, Edit Dist 1 join
"""
def test_tcr_join_edit1():
  import pandas as pd
  from tcrdist.join import join_by_dist
  from tcrdist.breadth import get_safe_chunk
  from tcrdist.rep_funcs import compute_pws_sparse
  import pwseqdist as pw

  my_metrics = { "cdr3_b_aa" : pw.metrics.nb_vector_editdistance}
  my_weights = { "cdr3_b_aa" : 1}
  my_kargs = {"cdr3_b_aa" :{'use_numba': True}}
  distances = compute_pws_sparse(
      df= v21df,
      df2 = v20df ,
      metrics = my_metrics,
      weights = my_weights,
      kargs = my_kargs,
      radius=1,
      cpu=2, 
      chunk_size=get_safe_chunk(v21df.shape[0], v20df.shape[0]), 
      store=False, 
      pm_pbar=True)
  csrmat = distances['tcrdist'] # it's called a tcrdist here, but was computed as edit distance 1
  left_join_df = join_by_dist(
      how = 'left',
      csrmat = csrmat,
      left_df = v21df,
      right_df = v20df,
      left_cols  = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      right_cols = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      left_suffix = '_x',
      right_suffix = '_y',
      max_n= 10,
      radius = 1)
  inner_join_df = join_by_dist(
      how = 'inner',
      csrmat = csrmat,
      left_df = v21df,
      right_df = v20df,
      left_cols  = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      right_cols = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      left_suffix = '_x',
      right_suffix = '_y',
      max_n= 10,
      radius = 1)
  outer_join_df = join_by_dist(
      how = 'outer',
      csrmat = csrmat,
      left_df = v21df,
      right_df = v20df,
      left_cols  = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      right_cols = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      left_suffix = '_x',
      right_suffix = '_y',
      max_n= 10,
      radius = 1)
  
  assert left_join_df.shape[0] > inner_join_df.shape[0]
  assert outer_join_df.shape[0] > inner_join_df.shape[0]
  assert outer_join_df.shape[0] > left_join_df.shape[0]
    
"""
2. Full Example using TCRdist on all CDRs on real data from two MIRA cohorts, 
Here we use TCRrep to infer CDR1,2,2.5 to compute a full tcrdist 
"""
def test_tcr_join_tcrdist():
  import pandas as pd
  from tcrdist.breadth import get_safe_chunk
  from tcrdist.repertoire import TCRrep
  from tcrdist.join import join_by_dist

  tr20 = TCRrep(cell_df = v20df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'bio_identity','protein_coordinate']].copy(),
      organism='human', 
      chains=['beta'], 
      compute_distances = False)
  tr21 = TCRrep(cell_df = v21df[['subject', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'bio_identity', 'protein_coordinate']].copy(),
      organism='human', 
      chains=['beta'], 
      compute_distances = False)
  tr21.cpus = 2
  tr21.compute_sparse_rect_distances(df = tr21.clone_df, df2 = tr20.clone_df, radius = 36, chunk_size = get_safe_chunk(tr21.clone_df.shape[0], tr20.clone_df.shape[0]))

  left_right_comparision = join_by_dist(
      how = 'inner',
      csrmat = tr21.rw_beta,
      left_df = v21df,
      right_df = v20df,
      left_cols  = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      right_cols = ['cdr3_b_aa','v_b_gene','j_b_gene','protein_coordinate','bio_identity','subject'],
      left_suffix = '_x',
      right_suffix = '_y',
      max_n= 10,
      radius = 24)

