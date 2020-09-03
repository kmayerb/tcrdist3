"""
These are unit tests
"""
import numpy as np
import re

flatten = lambda l: [item for sublist in l for item in sublist]

"""
These are not tests per say but useful illustrations
"""
def allele_01(genename):
	"""
	>>> allele_01('TRBV19*02')
	'TRBV19*01'
	"""
	g,a = re.match(pattern = '(.*)([0-9])$', string= genename).groups()
	return(f"{g}1")


def test_olga_sample_beta():
	np.random.seed(1)
	from tcrdist.pgen import OlgaModel
	olga_model_beta = OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta")
	result = olga_model_beta.gen_cdr3s(V = 'TRBV20-1*01', J = 'TRBJ1-2*01', n = 5)
	assert isinstance(result, list)
	assert len(result) == 5
	assert result == ['CSARQGLANYGYTF','CSARPSRGQDGYTF','CSARDQRTGQDGYTF','CSARDVSSSGGYYGYTF','CSAPEPLTSGRACNGYTF']

def test_olga_sample_alpha():
	np.random.seed(1)
	from tcrdist.pgen import OlgaModel
	olga_model_alpha = OlgaModel(recomb_type="VJ", chain_folder = "human_T_alpha")
	result = olga_model_alpha.gen_cdr3s(V = 'TRAV19*01', J = 'TRAJ37*01', n = 5)
	assert isinstance(result, list)
	assert len(result) == 5
	assert result == ['CALSEAPGNTGKLIF','CAPPSGNTGKLIF','CALAGNTGKLIF','CAQDNTGKLIF','CALRNTGKLIF']

def test_olga_sample_alphas_for_a_human_repertoire():
	import re
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	import palmotif

	from tcrdist.pgen import OlgaModel
	olga_model_alpha = OlgaModel(recomb_type="VJ", chain_folder = "human_T_alpha")

	from tcrdist.pgen import OlgaModel
	olga_model_beta  = OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta")


	df = pd.read_csv("dash_human.csv")
	tr = TCRrep(cell_df = df, 
	            organism = 'human', 
	            chains = ['alpha','beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv')

	rb = [olga_model_beta.gen_cdr3s(V = allele_01(r['v_b_gene']), J = allele_01(r['j_b_gene']), n = 1) for _,r in tr.clone_df[['v_b_gene', 'j_b_gene']].iterrows()]
	ra = [olga_model_alpha.gen_cdr3s(V = allele_01(r['v_a_gene']), J = allele_01(r['j_a_gene']), n = 1) for _,r in tr.clone_df[['v_a_gene', 'j_a_gene']].iterrows()]
	
	# assert that we covered 95% of dash_human wiht our olga sampler 
	#assert len([x for x in flatten(rb) if x is not None]) / len(flatten(rb)) > .95
	#assert len([x for x in flatten(ra) if x is not None]) / len(flatten(ra)) > .95


def test_olga_sample_alphas_for_a_large_repertoire():
	import re
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	import palmotif

	from tcrdist.pgen import OlgaModel
	olga_model_beta_mouse = OlgaModel(recomb_type="VDJ", chain_folder = "mouse_T_beta")


	df = pd.read_csv("dash.csv")
	tr = TCRrep(cell_df = df, 
	            organism = 'mouse', 
	            chains = ['beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv')

	rb = [olga_model_beta_mouse.gen_cdr3s(V = allele_01(r['v_b_gene']), J = allele_01(r['j_b_gene']), n = 1) for _,r in tr.clone_df[['v_b_gene', 'j_b_gene']].iterrows()]
	# assert that we covered 95% of dash_human wiht our olga sampler 
	# assert len([x for x in flatten(rb) if x is not None]) / len(flatten(rb)) > .95
	


def test_olga_sample():
	np.random.seed(310)
	from tcrdist.pgen import OlgaModel
	olga_model_beta = OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta")
	result = olga_model_beta.gen_cdr3(V = 'TRBV20-1*01', J = 'TRBJ1-2*01')
	# NOTE: seed is set, so we expect standard result
	#NOTE: .gen_cdr3() returns the full output tuple
	expected = ('TGCAGTGCTAGAGTAAGGGAAGCGGGAAGGACCTACACCTTC',
				 'CSARVREAGRTYTF',
				 29,
				 1,
				 {'V': 29,
				  'D': 2,
				  'J': 1,
				  'delV': 5,
				  'delJ': 15,
				  'delDl': 10,
				  'delDr': 7,
				  'insVD': 7,
				  'insDJ': 6})
	assert result == expected
	
	# NOTE: .gen_cdr3s() returns a list of CDR3s (amino acid only)
	np.random.seed(310)
	result = olga_model_beta.gen_cdr3s(V = 'TRBV20-1*01', J = 'TRBJ1-2*01', n = 4)
	expected = ['CSARVREAGRTYTF', 'CSAVPPGLPNYGYTF', 'CSARGPSQGYVRGLYGYTF', 'CSAQGLAGYGYTF']
	assert result == expected



def motif_creation_human_betas():
	import re
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	import palmotif

	from tcrdist.pgen import OlgaModel
	oma = OlgaModel(recomb_type="VJ", chain_folder = "human_T_alpha")

	from tcrdist.pgen import OlgaModel
	omb = OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta")


	df = pd.read_csv("dash_human.csv")
	tr = TCRrep(cell_df = df, 
	            organism = 'human', 
	            chains = ['alpha','beta'], 
	            db_file = 'alphabeta_gammadelta_db.tsv')

	from tcrdist.adpt_funcs import get_basic_centroids
	get_basic_centroids(tr, max_dist = 75)
	with open("test_3.svg", 'w') as oh:
		oh.write('<body>')
		for i,r in tr.centroids_df.iterrows():
			if len(r['neighbors']) < 5:	
				break
			seqs = tr.clone_df.iloc[r['neighbors'],]['cdr3_b_aa'].to_list()
			gene_usages = tr.clone_df.iloc[r['neighbors'],][['v_b_gene', 'j_b_gene']].value_counts().reset_index().to_dict('split')['data']
			depth = 3

			refs = flatten([omb.gen_cdr3s(allele_01(v),allele_01(j),i*depth) for v,j,i in combos_alpha])
			refs = [x for x in refs if x is not None]
			
			matrix, stats = palmotif.compute_pal_motif(seqs = seqs , refs = refs, centroid = r['cdr3_b_aa'])
			matrix_raw, _ = palmotif.compute_pal_motif(seqs = seqs , centroid = r['cdr3_b_aa'])
			refs.append(r['cdr3_b_aa'])
			matrix_bkgd, _ = palmotif.compute_pal_motif(seqs = refs, centroid = r['cdr3_b_aa'])
			
			svgs = [palmotif.svg_logo(matrix, 'test.svg', return_str = True), 
					palmotif.svg_logo(matrix_raw, 'test.svg', return_str = True),
					palmotif.svg_logo(matrix_bkgd, 'test.svg', return_str = True)]

			[oh.write(f"{s}<div></div>\n") for s in svgs]
			oh.write('<div></div>')
			oh.write(str(r))
			oh.write('<div></div>')
		
		oh.write('</body>')



def sim_all_cdr3_gen():
	import itertools
	def expand_grid(dct):
	    rows = itertools.product(*dct.values())
	    return pd.DataFrame.from_records(rows, columns=dct.keys())


	from tcrdist.pgen import OlgaModel
	omb = OlgaModel(recomb_type="VDJ", chain_folder = "human_T_beta")
	all_possible_beta = expand_grid(
	    {'V': omb.pgen_model.V_allele_names,
	     'J': omb.pgen_model.J_allele_names}
	)

	find_nones = list()
	results = list()
	for i,r in all_possible_beta.iterrows():
		e = omb.gen_cdr3s(V = r['V'], J = r['J'], n = 3)
		results.append(e)
		if e is None:
			find_nones.append( [r['V'],r['J'], e])
		print((r['V'],r['J'],e))



	from tcrdist.pgen import OlgaModel
	oma = OlgaModel(recomb_type="VJ", chain_folder = "human_T_alpha")
	all_possible_alpha = expand_grid(
	    {'V': oma.pgen_model.V_allele_names,
	     'J': oma.pgen_model.J_allele_names}
	)

	find_nones = list()
	results = list()
	for i,r in all_possible_alpha.iterrows():
		e = oma.gen_cdr3(V = r['V'], J = r['J'])
		results.append([r['V'],r['J'], e])
		if e is None:
			find_nones.append( [r['V'],r['J'], e])
		print((r['V'],r['J'],e))

	# Things we can't find:
	df = pd.DataFrame(results, columns = ['v','j','r'])
	df[df['r'].isna()][['v']].value_counts()
	df[df['r'].isna()][['j']].value_counts()





