.. tip::
   
   NOTE: tcrdist3 is in a pre-release phase and these docs are currently being drafted. 
   All code snippets are actual working unit and integration tests, but some 
   described features may not be available until the 1.0.0 release (Expected August 31st, 2020)

tcrdist3
========

`tcrdist3 <https://github.com/kmayerb/tcrdist3>`_ is an open-source python package that enables a broad array of T cell receptor sequence analyses.

In addition to the docstrings found throughout the code, we provide an example in a minimal TCR analysis workflow. A gallery of runnable examples demonstrates the suite of features and more complext analyses.


Installation
------------
.. code-block:: none

   pip install git+https://github.com/kmayerb/tcrdist3.git@0.1.3


Brief background
----------------

A human alpha/beta T cell repertoire is huge, comprised of an estimated 10^11 to 10^13 total T cells, with the process of V(D)J recombination capable of generating as many as 10^15 to 10^20 unique receptor sequences. With the wide adoption of high-throughput Adaptive Immune Receptor Repertoire sequencing (AIRR-seq) and single-cell sequencing technologies, immunologists are witnessing an explosion in the availability of TCR repertoire data.

TCR functional diversity is concentrated in complementarity defining regions (CDRs) that contact antigen-presenting molecules and enable receptor docking and antigenic specificity. Although software tools exist for TCR repertoire analysis we have developed a modern Python package that improves upon a biochemically aware TCR distance metric that can effectively reduce the dimensionality of repertoire data, enable statistical analyses that efficiently handle the scale of and molecular diversity in AIRR-seq datasets. Tcrdist3 is designed as an interactive set of functional operations that can be interfaced with existing tools and used to construct flexible hypothesis-driven analysis pipelines within batch computing frameworks.

Loading a TCR Dataset
---------------------

Loading a TCR dataset from most standard formats, including AIRR, MIXCR and Adaptive, and preparing it for analysis with tcrdist3 is simple. Tcrdist3 operates with a `pandas.DataFrame`, which can be created from CSV/TSV and many other tabular formats. Once the DataFrame is loaded its a matter of renaming columns to match tcrdist or telling tcrdist which columns to use.

[example of loading from a `.csv file <https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/tcrdist/test_files_compact/dash.csv>`_ ].

.. literalinclude:: ../tcrdist/tests/test_introduction_1.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python


For analyses that require the sequences of the CDRs encoded by the V genes, 
these sequences can be inferred from gene/allele names based on IMGT nomenclature 
(the nomenclature used by Adaptive Biotechnology is also supported as shown in 
an example on :ref:`loading_adaptive_biotechnology`).

Gene Level Analyses
-------------------

Loading a TCR dataset will immediately enable gene-level analyses of repertoires, 
for example Sankey plots of V/D/J gene frequency or statistical tests for differential 
gene enrichment in two or more conditions:


.. literalinclude:: ../tcrdist/tests/test_introduction_3.py
    :linenos:
    :lines: 16-28
    :dedent: 4
    :language: python


PA and NP gene usage diagrams. Click to enlarge.

|pic1| |pic2|

.. |pic1| image:: PA_gene_usage_plot.svg
   :width: 49%

.. |pic2| image:: NP_gene_usage_plot.svg
   :width: 49%


.. literalinclude:: ../tcrdist/tests/test_introduction_3.py
    :linenos:
    :lines: 30-50
    :dedent: 4
    :language: python



+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
|     | xcol    | xval | ycol     | yval       | X+Y+  | X+Y-  | X-Y+  | X-Y-  | X_marg   | Y_marg   | X|Y+     | X|Y-     | Y|X+     | Y|X-     | OR        | pvalue   |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
| 0   | epitope | NP   | j_b_gene | TRBJ1-6*01 | 48.0  | 255.0 | 2.0   | 963.0 | 0.238959 | 0.039432 | 0.960000 | 0.209360 | 0.158416 | 0.002073 | 90.635294 | 0.000001 |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
| 1   | epitope | NP   | j_b_gene | TRBJ2-7*01 | 52.0  | 251.0 | 146.0 | 819.0 | 0.238959 | 0.156151 | 0.262626 | 0.234579 | 0.171617 | 0.151295 | 1.162146  | 0.414406 |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
| ..  | ...     | ...  | ...      | ...        | ...   | ...   | ...   | ...   | ...      | ...      | ...      | ...      | ...      | ...      | ...       | ...      |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
| 103 | epitope | PB1  | v_b_gene | TRBV3*01   | 10.0  | 631.0 | 3.0   | 624.0 | 0.505521 | 0.010252 | 0.769231 | 0.502789 | 0.015601 | 0.004785 | 3.296355  | 0.090888 |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+
| 104 | epitope | PB1  | v_b_gene | TRBV30*01  | 5.0   | 636.0 | 0.0   | 627.0 | 0.505521 | 0.003943 | 1.000000 | 0.503563 | 0.007800 | 0.000000 | inf       | 0.062083 |
+-----+---------+------+----------+------------+-------+-------+-------+-------+----------+----------+----------+----------+----------+----------+-----------+----------+



.. toctree::
   :caption: Workflow
   :maxdepth: 1

   workflow
   

.. toctree::
   :caption: Gallery
   :maxdepth: 1

   tcrdistances
   visualizing
  
   

.. toctree::
   :caption: More Details
   :maxdepth: 1

   inputs
   adaptive
   

   


