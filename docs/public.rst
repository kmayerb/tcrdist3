.. _public:

(Quasi)Public Clones
====================


**Public TCRs** are shared clonotypes 
found in multiple individuals, arising 
from VDJ recombination biases and common selection pressures. 
Repertoire analyses often focuses on public clones; however
finding public antigen-specific TCRs is not always possible 
because TCR repertoires are characterized by extreme
diversity. As a consequence, only a small fraction of the repertoire 
can be assayed in a single sample,
making it difficult to reproducibly sample TCR clonotypes 
from an individual, let alone reliably detect shared clonotypes in a population. 

Enter, stage left, the **quasi-public TCRs** -- two or more TCRs, with a 
high degree of biochemical similarity -- that are found in two or more individuals. 
Identifying quasi public TCRs becomes useful when evaluating an
antigen enriched repertoire putatively recognizing the same epitope. 

Finding similar receptors from multiple individuals provides
stronger evidence of shared epitope recognition and
reveals mechanistic basis for CDR-peptide-MHC binding.

Moreover, meta-clonotypes are by definition more abundant 
than exact clonotype and thus can be more reliably be
detected in a single bulk unenriched sample,
facilitating more robust function comparisons across populations.


I am happy to use the defaults
------------------------------

For instance, you may want find all the (quasi)public
collections of TCRs within a fixed radius <= 18 
TCRdist units of each TCR in the antigen enriched 
input data.

.. literalinclude:: ../tcrdist/tests/test_example_17.py
    :linenos:
    :lines: 6-70
    :dedent: 4
    :language: python

In addition to the summary DataFrames returned, 
a HTML `quasi-publicity report <_static/quasi_public_clones.html>`_ is generated, 
allowing for the inspection of logo-motifs formed from highly similar 
antigen-enriched TCR sequences found in multiple subjects.


I'd like to tweak a default parameter
-------------------------------------

If you want to add or subtract information from the report you can do 
so relatively easily. For instance, suppose you want to 
summarize cohort information and add that to the report.

Just like elsewhere in tcrdist3, the python object `TCRpublic`
stores all options as attributes:

.. code-block:: python

	In [2]: tp.__dict__
	Out[2]:
	{'tcrrep': <tcrdist.repertoire.TCRrep at 0x13d5b9310>,
	 'organism': 'human',
	 'chain': 'beta',
	 'output_html_name': 'quasi_public_clones.html',
	 'pw_mat_str': 'pw_beta',
	 'cdr3_name': 'cdr3_b_aa',
	 'v_gene_name': 'v_b_gene',
	 'j_gene_name': 'j_b_gene',
	 'nr_filter': True,
	 'labels': ['clone_id',
	  'cdr3_b_aa',
	  'v_b_gene',
	  'j_b_gene',
	  'radius',
	  'neighbors',
	  'K_neighbors',
	  'nsubject',
	  'qpublic',
	  'cdr3_b_aa.summary',
	  'v_b_gene.summary',
	  'j_b_gene.summary',
	  'cdr3_b_aa.summary',
	  'subject.summary'],
	 'fixed_radius': False,
	 'radius': None,
	 'query_str': 'qpublic == True & K_neighbors > 5',
	 'kargs_member_summ': {'key_col': 'neighbors',
	  'count_col': 'count',
	  'addl_cols': ['subject'],
	  'addl_n': 4},
	 'kargs_motif': {'pwmat_str': 'pw_beta',
	  'cdr3_name': 'cdr3_b_aa',
	  'v_name': 'v_b_gene',
	  'gene_names': ['v_b_gene', 'j_b_gene']},
	 'tcrsampler': <tcrsampler.sampler.TCRsampler at 0x13d5b9850>}


The default only summarizes subjects 'addl_cols': ['subject'], so adding an additional 
categorical variable to include in the summary is as easy as:

.. code-block:: python

	tp.kargs_member_summ['addl_cols'] = ['subject', 'cohort']
	tp.labels.append("cohort.summary")


You can also specify your standard for publicity. Instead of 
'qpublic == True & K_neighbors > 5' you can ask to
find super public meta-clonotypes, returning only those 
groups that satisfy : 'nsubject > 8'  

.. code-block:: python

	tp.query_str = 'nsubject > 8'  


Here's the a full example:

.. literalinclude:: ../tcrdist/tests/test_example_18.py
    :linenos:
    :lines: 6-70
    :dedent: 4
    :language: python


As you can see in this new a html `quasi-publicity report <_static/quasi_public_clones2.html>`_ , 
the report has a new column for summarizing the percentage of TCRs coming from 
each cohort in the study and the number of meta-clonotypes are fewer, since 
only those with TCRs drawn from more than 8 subject are reported. 


I want my search radius to be sequence specific
-----------------------------------------------

The radius applied to each centroid can be specified 
in a column of the clone_df.

.. literalinclude:: ../tcrdist/tests/test_example_19.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python

Notice that radius varies by row in this `quasi-publicity report <_static/quasi_public_clones3.html>`_ , 

Working from neighbor_diff output
+++++++++++++++++++++++++++++++++

.. literalinclude:: ../tcrdist/tests/test_example_21.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python


I hate OOP just show me the functions
-------------------------------------

`TCRpublic` is for convenience. You can customize a lot including 
the background tcrsampler; but the power users may want to work with the 
underlying functions directly. Here are some examples of how:


I just want to quickly find neighbors and (quasi)public clones
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. literalinclude:: ../tcrdist/tests/test_example_20.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python



I have neighbors and radii already, I want logos
++++++++++++++++++++++++++++++++++++++++++++++++

Suppose you want to specify exactly what to include in a 
motif logo report. This example is slightly different 
then those above because we are going to use two 
inputs files. The first input file includes all
of the TCRs in antigen enriched repertoire. The
second file is a subset of the first, 
specifying exactly the TCRs centroids to 
report. Remember that any element of a 
clone_df can be included/excluded from the 
HTML report. Those fields to include can be 
specified as labels.


.. literalinclude:: ../tcrdist/tests/test_example_22.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python


Will this work with sparse matrix options?
++++++++++++++++++++++++++++++++++++++++++

tcrdist3 has a memory efficient options for larger datasets 
that produce scipy.sparse rather than dense representations 
of distance relationships. 

Currently you can't call TCRpublic() on this sparse representation. 
However, here is an example of how you can achieve similar results
via a script, reporting (quasi)Public meta-clonotypes from a sparse 
format.

.. literalinclude:: ../tcrdist/tests/test_example_23.py
    :linenos:
    :lines: 6-150
    :dedent: 4
    :language: python

For more on sparse matrices in tcrdist3 see the tab on 'Working With Bulk Data'.











