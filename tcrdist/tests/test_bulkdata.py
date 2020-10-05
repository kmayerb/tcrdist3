import pytest

def test_quick_pipeline_with_fragmented_compute():

	"""
	How can I used tcrdist3 to test for TCRs that may HLA restricted. 

	
	"""

	import os
	import pandas as pd
	import numpy as np
	from scipy import sparse
	from tcrdist.repertoire import TCRrep
	from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory
	
	f = 'mira_epitope_67_382_APHGVVFL_APHGVVFLHV_GVVFLHVTY_VVFLHVTYV.tcrdist3.csv'
	f = os.path.join('tcrdist','data','covid19',f)
	assert os.path.isfile(f)

	df = pd.read_csv(f)
	df = df[['subject', 'cell_type', 'v_b_gene', 'j_b_gene', 'cdr3_b_aa', 'cdr3_b_nucseq',  'cohort', 'hla-a', 'hla-a_1','hla-b', 'hla-b_1']]
	tr = TCRrep(cell_df = df,               
				organism = 'human',
				chains = ['beta'],
				db_file = 'alphabeta_gammadelta_db.tsv',
				compute_distances = False,
				store_all_cdr = False)

	from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory
	
	S, fragments = compute_pw_sparse_out_of_memory(	tr = tr,
													row_size      = 100,
													pm_processes  = 2,
													pm_pbar       = True,
													max_distance  = 1000,
													matrix_name   = 'rw_beta',
													reassemble    = True,
													cleanup       = False)

	tr.clone_df['B07'] = (tr.clone_df['hla-b'].str.startswith("B*07") | tr.clone_df['hla-b_1'].str.startswith("B*07"))
	tr.clone_df['B07'] = ["B*07" if (x) else "NOTB*07 " for x in tr.clone_df['B07']]

	#sparse.save_npz("S.npz", S)
	from tcrdist.rep_funcs import  compute_n_tally_out_of_memory
	nn_tally_df_cohort = compute_n_tally_out_of_memory(fragments,
												matrix_name = "rw_beta",
												pm_processes  = 6,
												to_file = False,
												to_memory = True, 
												knn_radius = 25, 
												x_cols = ['B07'])

	from hierdiff.association_testing import cluster_association_test
	nn_associations = cluster_association_test(res = nn_tally_df_cohort, y_col='cmember', method='fishers')
	nn_associations = nn_associations.sort_values('pvalue', ascending = True)
	import ast 
	nn_associations['neighbors_i'] = nn_associations.neighbors.apply(lambda x: ast.literal_eval(x))

	from tcrdist.summarize import test_for_almost_subsets, filter_is, filter_gt
	nn_associations['mostly_unique'] = test_for_almost_subsets(nn_associations['neighbors_i'], thr = 5)
	nr_nn_associations = filter_is(nn_associations, 'mostly_unique', 1).copy()

	#nr_nn_associations = filter_gt(nr_nn_associations, 'K_neighbors', 25).copy()
	nr_nn_associations


	# MOTIF GENERATION
	from tcrsampler.sampler import TCRsampler
	t = TCRsampler()
	if  'olga_human_beta_t.sampler.tsv' not in t.currently_available_backgrounds():
		t.download_background_file('olga_sampler.zip')
	#t.download_background_file('olga_sampler.zip') # ONLY IF NOT ALREADY DONE
	tcrsampler_beta = TCRsampler(default_background = 'olga_human_beta_t.sampler.tsv')
	tcrsampler_beta.build_background(max_rows = 1000)

	"""SEE PALMOTIF DOCS (https://github.com/agartland/palmotif)"""
	from palmotif import compute_pal_motif, svg_logo
	from tcrdist.summarize import _select
	
	"""GENERATE SVG GRAPHIC FOR EACH NODE OF THE TREE"""
	#pwmat_str = 'pw_beta'
	cdr3_name = 'cdr3_b_aa'
	gene_names = ['v_b_gene','j_b_gene']
	svgs_beta = list()
	svgs_beta_raw = list()
	info_list = list()

	from tcrdist.rep_diff import member_summ
	summary = member_summ(  res_df = nr_nn_associations,
							clone_df = tr.clone_df,
							addl_cols=['cohort','hla-a', 'hla-a_1', 'hla-b', 'hla-b_1', 'subject'])

	nr_nn_associations = pd.concat([nr_nn_associations, summary], axis = 1).reset_index()

	for i,r in nr_nn_associations.head(25).iterrows():
		dfnode  = tr.clone_df.iloc[r['neighbors_i'],:].copy()
		# <pwnode> Pairwise Matrix for node sequences
		pwnode = S[r['neighbors_i'],:] [:,r['neighbors_i']].todense()
		if dfnode.shape[0] > 2:
			iloc_idx = pwnode.sum(axis = 0).argmin()
			centroid = dfnode[cdr3_name].to_list()[iloc_idx]
		else:
			centroid = dfnode[cdr3_name].to_list()[0]

		print(f"CENTROID: {centroid}")

		gene_usage_beta = dfnode.groupby(gene_names).size()
		sampled_rep = tcrsampler_beta.sample( gene_usage_beta.reset_index().to_dict('split')['data'],
			flatten = True, depth = max(100, 1000 // dfnode.shape[0]))

		sampled_rep  = [x for x in sampled_rep if x is not None]

		motif, stat = compute_pal_motif(
						seqs = _select(df = tr.clone_df,
									   iloc_rows = r['neighbors_i'],
									   col = cdr3_name),
						refs = sampled_rep,
						centroid = centroid)

		svgs_beta.append(svg_logo(motif, return_str= True))

		sampled_rep = sampled_rep.append(centroid)
		motif_raw, _ = compute_pal_motif(
					 seqs =_select(df = tr.clone_df,
									iloc_rows = r['neighbors_i'],
									col = cdr3_name),
					 centroid = centroid)
		svgs_beta_raw.append(svg_logo(motif_raw, return_str= True))
		info_list.append(r)


	def row_to_string(r, vals = ['ct_columns', 'val_0', 'ct_0', 'val_1', 'ct_1', 'val_2', 'ct_2','val_3', 'ct_3', 'levels', 'K_neighbors', 'R_radius', 'RR', 'OR', 'pvalue', 'FWERp','FDRq']):
		#d = {v:r[v] for v in vals}
		return "<br></br>".join([f"\t{v} : {r[v]}" for v in vals])

	def to_html_table(r, vals = ['ct_columns', 'hla-a', 'hla-a_1', 'hla-b', 'hla-b_1', 'val_0', 'ct_0', 'val_2', 'ct_2', 'K_neighbors', 'R_radius', 'pvalue', 'FDRq','cdr3_b_aa','v_b_gene', 'j_b_gene', 'cohort','subject']):
		return pd.DataFrame(r[vals]).transpose().to_html()

	def shrink(html_str):
		return html_str.replace('height="100%"',  'height="10%"').\
			replace('width="100%"', 'width="10%"')

	with open('svgs_in_line.html', 'w') as fh:
		fh.write(f"<html><body>\n")
		

		for svg, svg_raw, details in zip(svgs_beta, svgs_beta_raw, info_list):
			fh.write(f"{shrink(svg_raw)}{shrink(svg)}")
			try:
				fh.write(to_html_table(details))
			except:
				print("F")
			fh.write("<div></div>")
		fh.write(f"</html></body>\n")