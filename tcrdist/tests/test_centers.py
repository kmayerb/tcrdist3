"""
test_centers.py
"""
def test_calc_radii():
	import numpy as np
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index(drop = True)
	tr = TCRrep(cell_df = df.copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
	            compute_distances = True)

	df = pd.read_csv("dash.csv").query('epitope != "PA"').reset_index(drop = True)
	tr_bkgd = TCRrep(cell_df = df.copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
               compute_distances = False)

	from tcrdist.centers import calc_radii
	radii, thresholds, ecdfs = calc_radii(tr = tr, tr_bkgd = tr_bkgd, chain = 'beta', ctrl_bkgd = 10**-3, use_sparse = False, max_radius=50)
	from tcrdist.public import _neighbors_variable_radius
	# Compute neighbors <= variable radius in the background set and the foreground set
	neighbors     = _neighbors_variable_radius(pwmat = tr.pw_beta , radius_list = radii)
	background_neighbors = _neighbors_variable_radius(pwmat = tr.rw_beta , radius_list = radii)
	tr.clone_df['radius']               = radii
	tr.clone_df['neighbors']            = neighbors
	tr.clone_df['background_neighbors'] = background_neighbors
	tr.clone_df['nsubject']             = tr.clone_df['neighbors'].\
			apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
	tr.clone_df['qpublic']              = tr.clone_df['nsubject'].\
			apply(lambda x: x > 1)
	tr.clone_df


def test_calc_radii_if_big():
	import numpy as np
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	from tcrdist.centers import calc_radii
	from tcrdist.public import _neighbors_sparse_variable_radius
	df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index(drop = True)
	tr = TCRrep(cell_df = df.copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
	            compute_distances = False)
	# For a large matrix one can use compute_sparse_rect_distances() instead of
	# .compute_distances() for the pairwise square matrix
	tr.cpus = 2
	tr.compute_sparse_rect_distances(df = tr.clone_df, radius=50,chunk_size=100)
	tr.pw_beta = tr.rw_beta.copy()

	# 
	dfb = pd.read_csv("dash.csv").query('epitope != "PA"').reset_index(drop = True)
	tr_bkgd = TCRrep(cell_df = dfb.copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False)
	
	# Set rw_beta to none as it will be computed between target and background by calc_radii
	tr.rw_beta = None
	from tcrdist.centers import calc_radii
	radii, thresholds, ecdfs = calc_radii(tr = tr, tr_bkgd = tr_bkgd, chain = 'beta', ctrl_bkgd = 10**-3, use_sparse = True, max_radius=50)
	tr.clone_df['radius'] = radii

	from tcrdist.public import _neighbors_sparse_variable_radius
	# Compute neighbors <= variable radius in the background set and the foreground set
	neighbors     = _neighbors_sparse_variable_radius(csrmat = tr.pw_beta , radius_list = tr.clone_df['radius'])
	background_neighbors = _neighbors_sparse_variable_radius(csrmat = tr.rw_beta , radius_list = tr.clone_df['radius'])
	#tr.clone_df['radius']               = radii
	tr.clone_df['neighbors']            = neighbors
	tr.clone_df['background_neighbors'] = background_neighbors
	tr.clone_df['nsubject']             = tr.clone_df['neighbors'].\
			apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
	tr.clone_df['qpublic']              = tr.clone_df['nsubject'].\
			apply(lambda x: x > 1)
	tr.clone_df



def test_example_with_report():
	# Example that would work with a large bakcgournd
	import numpy as np
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	from tcrdist.background import sample_britanova
	from tcrdist.sample import _default_sampler
	"""A useful background for beta chain"""
	df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index(drop = True)

	tr = TCRrep(cell_df = df.copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
	            compute_distances = True)
	
	ts = _default_sampler(organism = "mouse", chain = "beta")()
	trb = TCRrep(cell_df = ts.ref_df.rename(columns = {'v_reps' : 'v_b_gene', 'j_reps': 'j_b_gene', 'cdr3': 'cdr3_b_aa'}).copy(), 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv',
	            compute_distances = False)

	tr.cpus = 2
	tr.compute_sparse_rect_distances(df = tr.clone_df, df2 = trb.clone_df, radius=50,chunk_size=100)
	
	from tcrdist.centers import calc_radii
	radii, thresholds, ecdfs = calc_radii(tr = tr, tr_bkgd = trb, chain = 'beta', ctrl_bkgd = 10**-5, use_sparse = True, max_radius=50)

	# Set a maximum radius of 26
	tr.clone_df['radius'] = radii
	tr.clone_df['radius'][tr.clone_df['radius'] > 26] = 26

	# Quick access to publicity
	from tcrdist.public import _neighbors_sparse_variable_radius, _neighbors_variable_radius
	tr.clone_df['neighbors'] = _neighbors_variable_radius(pwmat = tr.pw_beta, radius_list = tr.clone_df['radius'])
	tr.clone_df['background_neighbors'] = _neighbors_sparse_variable_radius(csrmat = tr.rw_beta, radius_list = tr.clone_df['radius'])
	tr.clone_df['nsubject']             = tr.clone_df['neighbors'].\
			apply(lambda x: tr.clone_df['subject'].iloc[x].nunique())
	tr.clone_df['qpublic']              = tr.clone_df['nsubject'].\
			apply(lambda x: x > 1)
	tr.clone_df

	# A Report
	from tcrdist.public import TCRpublic
	tp = TCRpublic(
	tcrrep = tr, 
	output_html_name = "quasi_public_clones.html")
	tp.fixed_radius = False
	rp = tp.report()





def alternative_actually_a_big_case():
	import os
	import numpy as np
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	from tcrdist.neighbors import compute_population_estimate_ecdf
	from tcrdist.ecdf import distance_ecdf
	project_path = os.path.join('tutorial')
	source_path = os.path.join('tcrdist','data','covid19')
	antigen_enriched_background_file = 'mira_epitope_55_524_ALRKVPTDNYITTY_KVPTDNYITTY.tcrdist3.csv'
	assert os.path.isfile(os.path.join(source_path, antigen_enriched_background_file))
		# Read file into a Pandas DataFrame <df>
	df = pd.read_csv(os.path.join(source_path, antigen_enriched_background_file))
		# Drop cells without any gene usage information
	df = df.query("v_b_gene.notna() & j_b_gene.notna()")
		# Provide a counts column if non is present
	df['count'] = 1
		# Initialize a TCRrep class, using ONLY columns that are complete and unique define a a clone.
		# Counts of identical 'clones' will be aggregated into a TCRrep.clone_df.
	tr = TCRrep(cell_df = df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'subject','count']], 
				organism = "human", 
				chains = ['beta'], 
				compute_distances = True,
				cpus = 6)

	df_bkgd = pd.read_csv(os.path.join(project_path, f"{antigen_enriched_background_file}.olga100K_brit100K_bkgd_2.csv"), sep = ",")
		# Load the background to a TCRrep without computing pairwise distances (i.e., compute_distances = False)
	tr_bkgd = TCRrep(cell_df = df_bkgd, organism = "human", chains = ['beta'], compute_distances = False)

	from tcrdist.centers import calc_radii

	radii = calc_radii(tr = tr, tr_bkgd = tr_bkgd,  ctrl_bkgd = 10**-5, use_sparse = True, max_radius=50, chunk_size=100)
	radii2 = calc_radii(tr = tr, tr_bkgd = tr_bkgd,  ctrl_bkgd = 10**-5, use_sparse = False, max_radius = 50)
	assert radii == radii2

	tr.clone_df['max_radii'] = radii


	






	# # Set to None so it is recomputed as 
	# tr.rw_beta = None
	# radii2, thresholds2, ecdfs2 = calc_radii(tr = tr, tr_bkgd = tr_bkgd, chain = 'beta', ctrl_bkgd = 10**-3, use_sparse = True, max_radius=50, chunk_size = 100)
	# assert len(radii) == tr.clone_df.shape[0]
	# assert np.all(radii == radii2)
	# assert np.all(thresholds == thresholds2)
	# assert np.all(ecdfs  == ecdfs2)
	# from tcrdist.public import _neighbors_sparse_variable_radius
	# background_neighbors2 = _neighbors_sparse_variable_radius(csrmat = tr.rw_beta, radius_list = radii2, maxd= 50)
	# np.all(background_neighbors == background_neighbors2)


















	# # Example with a report
	# import numpy as np
	# import pandas as pd
	# from tcrdist.repertoire import TCRrep
	# df = pd.read_csv("dash.csv").query('epitope == "PA"').reset_index(drop = True)
	# tr = TCRrep(cell_df = df.copy(), 
	#             organism = 'mouse', 
	#             chains = ['beta'], 
	#             db_file = 'alphabeta_gammadelta_db.tsv',
	#             compute_distances = True)

	# df = pd.read_csv("dash.csv").query('epitope != "PA"').reset_index(drop = True)
	# tr_bkgd = TCRrep(cell_df = df.copy(), 
	#             organism = 'mouse', 
	#             chains = ['beta'], 
	#             db_file = 'alphabeta_gammadelta_db.tsv',
	#             compute_distances = False)
	# radii, thresholds, ecdfs = calc_radii(tr = tr, tr_bkgd = tr_bkgd, chain = 'beta', ctrl_bkgd = 10**-5, use_sparse = False, max_radius=50)
	# tr.clone_df['radius']            = radii
	# from tcrdist.public import TCRpublic
	# tp = TCRpublic(
	# 	tcrrep = tr, 
	# 	output_html_name = "quasi_public_clones.html")
	# tp.fixed_radius = False
	# rp = tp.report()






	

	# # Find 













