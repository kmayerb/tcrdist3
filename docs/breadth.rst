.. _beradth:

Breadth Estimation
==================

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
can be specified in with the argument `breadth_cols`. Crucially, when 
running tcrdist.breadth.long_form_tabulation() the argument 
`search_cols` must include a column indicating the association between 
a metac-clonotype and a particular 'protein' or 'epitope'  
e.g., ['tag', 'protein', 'epitope', cdr3_b_aa', 'v_b_gene', 'j_b_gene', 
'pgen','regex', 'radius']

Note that breadth and depth follow the definitions in Synder et al. 2020

`Synder et al. 2020 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418734/>`_, 
See section: The breadth and depth of a disease-specific T-cell response. 

.. attention::

	REALLY IMPORTANT: 
	This is a live unit test, where `testing_only` is set to True, using only 
	a sample of metaclonotypes. However, 
	be sure to set `testing_only` to False when you 
	run your own meta-clonotypes 


.. tip::
   
   Notice that this example first shows how to do this for a single file. 
   Then in proceeds to an example showing 
   how to estmate breadth for a list of files in a loop. The default 
   is to write out the breadth estimates as flat files, but you could 
   store these in a list and pd.concat() them together


.. tip::
   
   The metaclonotype downloadable in this example come from 
   Mayer-Blackwell K, Schattgen S, Cohen-Lavi L, Crawford JC, Souquette A, 
   Gaevert JA, Hertz T, Thomas PG, Bradley PH, Fiore-Gartland A. 2020 
   TCR meta-clonotypes for biomarker discovery with tcrdist3: identification of public, 
   LA-restricted SARS-CoV-2 associated TCR features. bioRxiv (2020)
   `doi:10.1101/2020.12.24.424260 <https://www.biorxiv.org/content/10.1101/2020.12.24.424260v2>`_
   


.. literalinclude:: ../tcrdist/tests/test_breadth.py
    :linenos:
    :lines: 55-300
    :dedent: 4
    :language: python