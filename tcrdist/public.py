"""
kmayerb 
2020-11-12

Assist the user in finding quasi public meta-clonotypes
"""
import pandas as pd
import numpy as np
import os
import scipy.sparse
import warnings
from palmotif import compute_pal_motif, svg_logo
from tcrdist.summarize import _occurs_N_str, member_summ, filter_is, test_for_subsets
from progress.bar import IncrementalBar
from tcrdist.tree import _default_sampler, _default_sampler_olga, allele_01

__all__ = ['_neighbors_fixed_radius',
		   '_K_neighbors_fixed_radius',
		   '_neighbors_variable_radius',
		   '_K_neighbors_variable_radius',
		   '_neighbors_sparse_fixed_radius',
		   '_neighbors_sparse_variable_radius',
		   'make_motif_logo',
		   'make_motif_logo_from_index',
		   '_quasi_public_meta_clonotypes']

def _neighbors_fixed_radius(pwmat, radius):
	""" Returns the list of neighbor column indices if within the fixed radius """
	return [list(np.nonzero(pwmat[ii,]<=radius)[0]) for ii in range(pwmat.shape[0]) ]

def _K_neighbors_fixed_radius(pwmat, radius):
	""" Returns the number of neighbors (self-inclusive) if within the fixed radius"""
	return [len(list(np.nonzero(pwmat[ii,]<=radius)[0])) for ii in range(pwmat.shape[0]) ]

def _neighbors_variable_radius(pwmat, radius_list):
	""" Returns the list of neighbor column indices if within the fixed radius """
	return [list(np.nonzero(pwmat[ii,]<=radius)[0]) for ii, radius in zip(range(pwmat.shape[0]), radius_list) ]

def _K_neighbors_variable_radius(pwmat, radius_list):
	""" Returns the number of neighbors (self-inclusive) if within the fixed radius"""
	return [len(list(np.nonzero(pwmat[ii,]<=radius)[0])) for ii, radius in zip(range(pwmat.shape[0]), radius_list) ]

def _neighbors_sparse_fixed_radius(csrmat, radius):
	""" 
	Returns the list of neighbor column indices within the fixed radius

	Parameter
	--------- 
	csrmat : scipy.sparse.csr_matrix

	radius : int

	Return
	------
	list of lists of ints

	Example 
	-------

	>>> M = np.array([[-1,0,1,4],[1,-1,0,2],[0,10,-1,3],[0,0,2,-1]])
	>>> S = scipy.sparse.csr_matrix(M)
	>>> NN = _neighbors_sparse_fixed_radius(csrmat = S, radius = 4)
	>>> NN
	[[0, 2], [0, 1, 3], [2, 3], [2, 3]]
	>>> assert NN == [[0, 2], [0, 1, 3], [2, 3], [2, 3]]
	"""
	S = csrmat.copy()
	S[S > radius] = 0
	S.eliminate_zeros()
	return [list(x) for x in np.split(S.indices, S.indptr[1:-1])]


def _neighbors_sparse_variable_radius(csrmat, radius_list, maxd= 50):
	""" 
	Returns the list of neighbor column indices within the fixed radius

	Parameter
	--------- 
	csrmat : scipy.sparse.csr_matrix
	
	radius : int

	Return
	------
	list of lists of ints

	Example 
	-------
	>>> from tcrdist.public import _neighbors_sparse_variable_radius
	>>> import scipy.sparse
	>>> M = np.array([[-1,1,1,4],[1,-1,1,2],[1,10,-1,3],[5,1,2,-1]])
	>>> S = scipy.sparse.csr_matrix(M)
	>>> N = _neighbors_sparse_variable_radius(csrmat = S, radius_list = [1,1,20,1])
	>>> N
	[[0, 1, 2], [0, 1, 2], [0, 1, 2, 3], [1,3]]
	>>> assert N ==  [[0, 1, 2], [0, 1, 2], [0, 1, 2, 3], [1,3]]
	"""
	S = csrmat.copy()
	NN = list()
	for i,radius in enumerate(radius_list):
		if radius == 0:
			#row = np.asarray(S[i, :].todense())[0]
			#row[row == 0] = radius + 1
			#x = np.nonzero(row == -1)[0]
			ind = np.nonzero(S.data[S.indptr[i]:S.indptr[i+1]] == -1)[0]
			col_indices = S.indices[S.indptr[i]:S.indptr[i+1]][ind].tolist()
			NN.append(col_indices)
		else:
			ind = np.nonzero(S.data[S.indptr[i]:S.indptr[i+1]] <= radius)[0]
			col_indices = S.indices[S.indptr[i]:S.indptr[i+1]][ind].tolist()
			NN.append(col_indices)
	return NN 



class TCRpublic():
	"""

	Attributes

	tcrrep : tcrdist.repertoire.TCRrep

	clone_df : pd.DataFrame
		Clones information with standard tcrdist3 column names.
	pwmat : np.array
		Pairwise distances
	tcrsamper : tcrsampler.TCRsampler
		TCRSampler instance initialized with appropriate background
		set.
	cdr3_name : str
		Column name for amino acid CDR3 e.g., 'cdr3_d_aa'.
	v_gene_name : str
		Column name for TR[ABGD]V gene e.g., 'v_d_gene'.
	nr_filter : bool
		If True, sequqences with the exact same neighbors as another set will be
		dropped
	output_html_name : str
		Filename for the output html output.
	labels : list
		List of columns to display on html output beneath each logo plot. 
	fixed_radius : False
		If False, clone_df must have a column radius. 
		If True, argument radius will be used to define 
		maximum distance from centroid to neighboring TCR.
	radius : int or None
		Theshold distance (<=) for neighborhood membership.
		If int, then all centroids will be assigned the same
		radius. Alterntively radius can be provided for each
		centroid sequence by including radius as a numeric column
		in clone_df.
	query_str : str
		The string to include sequences in output. For instance
		'qpublic == True and K_neighbors > 3', implies that only 
		grouping of 4 or more TCRs from at leåast two individuals
		will be retained. Alternatively, 'nsubject > 1' or 
		'qpublic == True' could be used as true minimum requirements 
		for quasi-publicity.
	kargs_member_summ : dict
		kwargs 
	kargs_motif : dict
		kwargs for the motif genertation 
	"""
	def __init__(self, tcrrep, organism = None, chain = None, fixed_radius = True, output_html_name = "quasi_public_clones.html"):
		
		self.tcrrep = tcrrep
		
		if organism is None:
			self.organism = tcrrep.organism 
		else:
			assert organism in ['human', 'mouse'], "TCRPublic <organism> argument must be 'human' or 'mouse'"
			self.organism = organism 
		
		if chain is None:
			self.chain= tcrrep.chains[0]
		else:
			assert chain in ['alpha','beta', 'gamma','delta'], "TCRPublic <chain> argument must be 'alpha', 'beta', 'gamma' or 'delta'"
			self.chain = chain

		self.output_html_name = output_html_name

		self.sort_columns = ['nsubject','K_neighbors']
		self.sort_ascending = False
		
		# Get chain specific atributes
		self.pw_mat_str = {'alpha': 'pw_alpha',
						  'beta' : 'pw_beta', 
						  'gamma': 'pw_gamma',
						  'delta': 'pw_delta'}[self.chain]
		
		self.cdr3_name = {'alpha': 'cdr3_a_aa',
						  'beta' : 'cdr3_b_aa', 
						  'gamma': 'cdr3_g_aa',
						  'delta': 'cdr3_d_aa'}[self.chain]

		self.v_gene_name = {'alpha': 'v_a_gene',
						   'beta'  : 'v_b_gene', 
						   'gamma' : 'v_g_gene',
						   'delta' : 'v_d_gene'}[self.chain]

		self.j_gene_name = {'alpha': 'j_a_gene',
						   'beta'  : 'j_b_gene', 
						   'gamma' : 'j_g_gene',
						   'delta' : 'j_d_gene'}[self.chain]

		self.nr_filter = True
		
		self.output_html_name = output_html_name

		self.labels = ['clone_id',
						self.cdr3_name, 
						self.v_gene_name,
						self.j_gene_name,
						'radius',
						'neighbors',
						'K_neighbors',
						'nsubject',
						'qpublic',
						f'{self.cdr3_name}.summary', 
						f'{self.v_gene_name}.summary',
						f'{self.j_gene_name}.summary',
						f'{self.cdr3_name}.summary',
						'subject.summary']
		
		self.fixed_radius = fixed_radius
		
		self.radius = 18
		
		self.query_str = 'nsubject > 3'
		
		self.kargs_member_summ = {
			'key_col'   : 'neighbors', 
			'count_col' : 'count',
			'addl_cols' : ['subject'],
			'addl_n'    : 4}
		
		self.kargs_motif = {
			'pwmat_str'  : self.pw_mat_str,
			'cdr3_name'  : self.cdr3_name,
			'v_name'     : self.v_gene_name,
			'gene_names' : [self.v_gene_name, self.j_gene_name]}

		# Get the Default sampler
		try: 
			self.tcrsampler = _default_sampler_olga(organism = self.organism, chain = self.chain)()
		except KeyError:
			self.tcrsampler = _default_sampler(organism = self.organism, chain = self.chain)()

	def report(self):

		result = _quasi_public_meta_clonotypes(clone_df    = self.tcrrep.clone_df.copy(), 
									  pwmat       = getattr(self.tcrrep, self.pw_mat_str),
									  tcrsampler  = self.tcrsampler,
									  cdr3_name   = self.cdr3_name,
									  v_gene_name = self.v_gene_name,
									  nr_filter   = self.nr_filter,
									  output_html_name  = self.output_html_name, 
									  sort_columns = self.sort_columns,
									  sort_ascending = self.sort_ascending,
									  labels            = self.labels,
									  fixed_radius      = self.fixed_radius,
									  radius            = self.radius,
									  query_str         = self.query_str, 
									  kargs_member_summ = self.kargs_member_summ,
									  kargs_motif       = self.kargs_motif)
		
		print(f"WRITING: {self.output_html_name}")
		return result




def make_motif_logo(tcrsampler,
					clone_df,
					pwmat,
					centroid = 'CASSPDIEKYF',
					v_gene = 'TRBV7-9*01',
					radius = 24,
					pwmat_str = 'pw_delta',
					cdr3_name = 'cdr3_d_aa',
					v_name = 'v_d_gene',
					gene_names = ['v_d_gene','j_d_gene']):

	"""
	Make a motif from a tcrrep clone_df, pwmat, and a tcrsampler. 

	Parameters
	----------
	tcrsampler : tcrsamper.TCRsampler,
	clone_df : pd.DataFrame,
	pwmat : np.array,
	centroid : str
		e.g.,'CASSPDIEKYF',
	v_gene : str
		e.g. 'TRBV7-9*01',
	radius = int
		e.g., 26,
	pwmat_str : str
		e.g.,'pw_delta',
	cdr3_name : str
		e.g., 'cdr3_d_aa',
	v_name : str
		e.g., 'v_d_gene',
	gene_names : list
		eg., ['v_d_gene','j_d_gene']

	Returns 
	-------
	svg : str
	svg_raw : str

	Notes
	-----
	
	There is a safety first, efficiency loss involved 
	since we are relocating neighbors that 
	may already be know, but by looking 
	up the row index <irow> fisrst matching V,CDR3 this 
	function can be evoked without knowing 
	anything about the positions of 
	the neighbors ahead of time. This is particularly useful 
	since clone_df order is not stable after groupby 
	and deduplication. 
	"""
	irow = clone_df[(clone_df[cdr3_name] == centroid) & (clone_df[v_name ] == v_gene)].index[0]

	dfnode = clone_df[pd.Series(pwmat[irow , :]) <= radius].copy()
	dfnode[gene_names[0]] = dfnode[gene_names[0]].apply(lambda x : allele_01(x))
	dfnode[gene_names[1]] = dfnode[gene_names[1]].apply(lambda x : allele_01(x))

	gene_usage = dfnode.groupby(gene_names).size()
	
	sampled_rep = tcrsampler.sample( gene_usage.reset_index().to_dict('split')['data'],
					flatten = True, depth = 100)

	sampled_rep  = [x for x in sampled_rep if x is not None]

	motif, stat = compute_pal_motif(
					seqs = dfnode[cdr3_name],
					refs = sampled_rep,
					centroid = centroid)

	svg = svg_logo(motif, return_str= True)

	motif_raw, _ = compute_pal_motif(
				seqs     = dfnode[cdr3_name],
				centroid = centroid)

	svg_raw = svg_logo(motif_raw, return_str= True)

	return svg, svg_raw


def make_motif_logo_from_index(tcrsampler,
							   ind,
							   clone_df,
							   centroid,
							   cdr3_name = 'cdr3_b_aa',
							   v_name = 'v_b_gene',
							   gene_names = ['v_b_gene','j_b_gene']):

	"""
	make motif logo from a specific index
	"""
	dfnode = clone_df.iloc[ind,:].copy()

	dfnode[gene_names[0]] = dfnode[gene_names[0]].apply(lambda x : allele_01(x))
	dfnode[gene_names[1]] = dfnode[gene_names[1]].apply(lambda x : allele_01(x))

	gene_usage = dfnode.groupby(gene_names).size()
	
	sampled_rep = tcrsampler.sample( gene_usage.reset_index().to_dict('split')['data'],
					flatten = True, depth = 100)

	sampled_rep  = [x for x in sampled_rep if x is not None]

	motif, stat = compute_pal_motif(
					seqs = dfnode[cdr3_name],
					refs = sampled_rep,
					centroid = centroid)

	svg = svg_logo(motif, return_str= True)

	motif_raw, _ = compute_pal_motif(
				seqs     = dfnode[cdr3_name],
				centroid = centroid)

	svg_raw = svg_logo(motif_raw, return_str= True)

	return svg, svg_raw




def _quasi_public_meta_clonotypes(clone_df, 
								  pwmat,
								  tcrsampler, 
								  cdr3_name = 'cdr3_d_aa',
								  v_gene_name = 'v_d_gene',
								  nr_filter = True,
								  output_html_name = "quasi_public_clones.html",
								  sort_columns = ['nsubject','K_neighbors'],
								  sort_ascending = False,
								  labels = ['clone_id',
											'cdr3_d_aa', 
											'v_d_gene',
											'j_d_gene',
											'radius',
											'neighbors',
											'K_neighbors',
											#'cdr3s',
											'nsubject',
											'qpublic',
											'cdr3_d_aa.summary',
											'v_d_gene.summary',
											'j_d_gene.summary',
											'subject.summary'],
								  fixed_radius = False,
								  radius = None,
								  query_str = 'qpublic == True & K_neighbors > 1',
								  kargs_member_summ = {
									'key_col'   : 'neighbors', 
									'count_col' : 'count',
									'addl_cols' : ['subject'],
									'addl_n'    : 4},
								  kargs_motif = {
									'pwmat_str'  : 'pw_delta',
									'cdr3_name'  : 'cdr3_d_aa',
									'v_name'     : 'v_d_gene',
									'gene_names' : ['v_d_gene','j_d_gene']}):
	
	"""
	_quasi_public_meta_clonotypes


	Parameters
	----------
	clone_df : pd.DataFrame
		Clones information with standard tcrdist3 column names.
	pwmat : np.array
		Pairwise distances
	tcrsamper : tcrsampler.TCRsampler
		TCRSampler instance initialized with appropriate background
		set.
	cdr3_name : str
		Column name for amino acid CDR3 e.g., 'cdr3_d_aa'.
	v_gene_name : str
		Column name for TR[ABGD]V gene e.g., 'v_d_gene'.
	nr_filter : bool
		If True, sequqences with the exact same neighbors as another set will be
		dropped
	output_html_name : str
		Filename for the output html output.
	labels : list
		List of columns to display on html output beneath each logo plot. 
	fixed_radius : False
		If False, clone_df must have a column radius. 
		If True, argument radius will be used to define 
		maximum distance from centroid to neighboring TCR.
	radius : int or None
		Theshold distance (<=) for neighborhood membership.
		If int, then all centroids will be assigned the same
		radius. Alterntively radius can be provided for each
		centroid sequence by including radius as a numeric column
		in clone_df.
	query_str : str
		The string to include sequences in output. For instance
		'qpublic == True and K_neighbors > 3', implies that only 
		grouping of 4 or more TCRs from at leåast two individuals
		will be retained. Alternatively, 'nsubject > 1' or 
		'qpublic == True' could be used as true minimum requirements 
		for quasi-publicity.
	kargs_member_summ : dict
		kwargs 
	kargs_motif : dict
		kwargs for the motif genertation 

	Returns
	-------
	Returns DataFrames in a Dictionary. 
		nn_summary : pd.DataFrame 
			DataFrame matchign clone_df with summary measures added
		quasi_public_df pd.DataFrame
			Dataframe with only those rows that match the <query_str> and nr_filter.
	{'nn_summary': nn_summary : pd.DataFrame, 
	 'quasi_public_df': quasi_public_df : nn_summary : pd.DataFrame}

	Notes
	-----
	Importantly a html file is written displaying the quasi-public meta-clonotypes
	
	The easiest way to integrate this with existing nieghbor_diff 
	add 'neighbors' and 'K_neighbors' to the clone df. Other columns 
	could be added as well, and then displayed if added to the 
	lis of labels.

	nn_clone_df = pd.concat([tr.clone_df, ndif[['neighbors', 'K_neighbors','val_0','ct_0']] ], axis = 1)

	Examples
	--------
	"""
	if 'neighbors' not in clone_df.columns:
		if fixed_radius:
			clone_df['radius'] = radius
			clone_df['neighbors']   = _neighbors_fixed_radius(pwmat = pwmat, radius = radius)
		else:
			assert 'radius' in clone_df.columns, "if not using fixed_radius, the clone_df must have a numeric 'radius' columns"
			clone_df['neighbors']   = _neighbors_variable_radius(pwmat = pwmat, radius_list = clone_df.radius)

	if 'K_neighbors' not in clone_df.columns:
		if fixed_radius:
			clone_df['K_neighbors'] = _K_neighbors_fixed_radius(pwmat = pwmat, radius = radius)
		else: 
			clone_df['K_neighbors'] = _K_neighbors_variable_radius(pwmat = pwmat, radius_list = clone_df.radius)


	if 'nsubject' not in clone_df.columns:
		clone_df['nsubject']   = clone_df['neighbors'].\
			apply(lambda x: clone_df['subject'].iloc[x].nunique())

	if 'qpublic' not in clone_df.columns:
		clone_df['qpublic']     = clone_df['nsubject'].\
			apply(lambda x: x > 1)

	nn_summary = member_summ(res_df   = clone_df,
							clone_df  = clone_df,
							**kargs_member_summ)
	nn_summary = nn_summary.rename(columns = {k:f'{k}.summary' for k in nn_summary.columns})

	clone_df['cdr3s']    = clone_df['neighbors'].apply(lambda x: clone_df[cdr3_name].iloc[x].to_list())
	

	clone_df = pd.concat([clone_df, nn_summary], axis = 1)

	quasi_public_df = clone_df.query(query_str).\
		sort_values(sort_columns, ascending = sort_ascending).\
		reset_index(drop = True).\
		copy()
	
	if quasi_public_df.shape[0] == 0:
		raise ValueError("UNFORTUNATELY NO QUASI PUBLIC CLOONES WERE FOUND, CONSIDER YOUR QUERY STRINGENCY")

	quasi_public_df['unique_set'] = test_for_subsets(quasi_public_df['neighbors'])

	if nr_filter:
		quasi_public_df = filter_is(quasi_public_df, 'unique_set', 1).reset_index(drop = True)

	print(f"GENERATING {quasi_public_df.shape[0]} QUASI-PUBLIC MOTIFS SATISFYING {query_str}")
	bar = IncrementalBar('Processing', max = quasi_public_df.shape[0])
	svgs = list()
	svgs_raw = list()
	for i,r in quasi_public_df.iterrows():
		bar.next()
		centroid = r[cdr3_name]
		v_gene   = r[v_gene_name]
		svg, svg_raw = make_motif_logo( tcrsampler = tcrsampler, 
										pwmat = pwmat,
										clone_df = clone_df,
										centroid = centroid ,
										v_gene = v_gene ,
										radius = r['radius'],
										**kargs_motif)
		svgs.append(svg)
		svgs_raw.append(svg_raw)
	bar.next();bar.finish()

	quasi_public_df['svg'] = svgs
	quasi_public_df['svg_raw'] = svgs_raw

	def shrink(s):
		s = s.replace('height="100%"', 'height="20%"')
		s = s.replace('width="100%"', 'width="20%"')
		return s

	print(labels)

	with open(output_html_name, 'w') as output_handle:
		for i,r in quasi_public_df.iterrows():
			#import pdb; pdb.set_trace()
			svg, svg_raw = r['svg'],r['svg_raw']
			output_handle.write("<br></br>")
			output_handle.write(shrink(svg))
			output_handle.write(shrink(svg_raw))
			output_handle.write("<br></br>")
			output_handle.write(pd.DataFrame(r[labels]).transpose().to_html())
			output_handle.write("<br></br>")

	return {'nn_summary': nn_summary, 'quasi_public_df': quasi_public_df, 'clone_df': clone_df}
