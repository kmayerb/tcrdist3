.. _join:


Tabulating Meta-Clonotypes
==========================

Step-by-Step Example
++++++++++++++++++++

This page offers a step-by-step explanation of how to tabulate 
meta-clonotype conformant clones in a bulk TCRb Repertoire using :py:func:`tcrdist.join.join_by_dist`.
For copying into a new example, all the code discussed below is contained in a single block (:ref:`full_example`).

.. code-block:: python
    :emphasize-lines: 8

    import multiprocessing
    import numpy as np
    import os 
    import pandas as pd
    from tcrdist.setup_tests import download_and_extract_zip_file
    from tcrdist.repertoire import TCRrep
    from tcrdist.breadth import get_safe_chunk
    from tcrdist.join import join_by_dist
    import re


Meta-clonotypes can be learned from antigen-associated TCR data
collected via tetramer sorting or activation marker enrichment. 
However, a single TCR may be conformant  multiple meta-clonotypes, 
which should be a consideration when applying a set of meta-clonotypes together.
For instance, tallying conformant TCRs in a repertoire should avoid counting 
a TCR more than once. This example illustrates an efficient approach to 
tabulating meta-clonotypes conformant sequences in bulk repertoires. 

Determine the appropriate number of CPUS based on your system

.. code-block:: python

    ncpus = min(multiprocessing.cpu_count(), 6)


This block of code downloads the files that will be used in this example.

* 2 bulk repertoires that have been preprocessed to include valid IMGT names and 
* Meta-clonotypes described in (Mayer-Blackwell et al. 2021) from the data v2.0 release of SARS-CoV-2 peptide associated CD8 TCRs (Nolan et. al 2020)
 

.. code-block:: python

    files = ['1588BW_20200417_PBMC_unsorted_cc1000000_ImmunRACE_050820_008_gDNA_TCRB.tsv.tcrdist3.tsv']
    if not np.all([os.path.isfile(f) for f in files]):
        download_and_extract_zip_file("ImmunoSeq_MIRA_matched_tcrdist3_ready_2_files.zip", source = "dropbox", dest = ".")
    if not os.path.isfile("bioRxiv_v2_metaclonotypes.tsv.zip"):
        download_and_extract_zip_file('bioRxiv_v2_metaclonotypes.tsv.zip', source = "dropbox", dest = ".")


We will consider meta-clonotypes as the search sequences, for instance:

+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+
| feature                                                               | v_b_gene   | cdr3_b_aa        | radius | regex                            | protein_coordinate | protein |
+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+
| M_55_1E6+TRBV28*01+CASSLRTDHYEQYF+22+(S[RLMF][RK][ST][ND].YEQ)        | TRBV28*01  | CASSLRTDHYEQYF   | 22     | (S[RLMF][RK][ST][ND].YEQ)        | ORF1ab 1316:1330   | ORF1ab  |
+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+
| M_55_1E6+TRBV28*01+CASSLRSDSYEQYF+10+(SL[RK][ST][ND]SYEQ)             | TRBV28*01  | CASSLRSDSYEQYF   | 10     | (SL[RK][ST][ND]SYEQ)             | ORF1ab1316:1330    | ORF1ab  |
+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+
| M_55_1E6+TRBV5-5*01+CASSPGQGAFTDTQYF+22+(S.G[QE]G[AS]F[ST]DTQ)        | TRBV5-5*01 | CASSPGQGAFTDTQYF | 22     | (S.G[QE]G[AS]F[ST]DTQ)           | ORF1ab 1316:1330   | ORF1ab  |
+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+
| M_55_1E6+TRBV28*01+CASSLKTDSYEQYF+16+(S[RLMF][RGKS][ST][ND][ANQS]YEQ) | TRBV28*01  | CASSLKTDSYEQYF   | 16     | (S[RLMF][RGKS][ST][ND][ANQS]YEQ) | ORF1ab 1316:1330   | ORF1ab  |
+-----------------------------------------------------------------------+------------+------------------+--------+----------------------------------+--------------------+---------+

.. code-block:: python

    df_search = pd.read_csv("bioRxiv_v2_metaclonotypes.tsv", sep = "\t")
    df_bulk = pd.read_csv(files[0], sep = "\t")
    
    df_bulk = df_bulk.sort_values('count').reset_index(drop = True)
    df_bulk['rank'] = df_bulk.index.to_list()


.. note::

    'rank' refers to the rank abundance of each clone and ensures that we get an accurate
    estimate of clonal breadth without collapsing sequences with identical at the amino acid 
    level 


.. code-block:: python
    
    from tcrdist.repertoire import TCRrep
    
    tr = TCRrep(
        cell_df = df_search, 
        organism = "human",
        chains = ['beta'], 
        compute_distances= False)
    tr.cpus = ncpus
    
    tr_bulk = TCRrep(
        cell_df = df_bulk, 
        organism = "human",
        chains = ['beta'], 
        compute_distances= False)  
    

This is the only computationally demanding step, where the TCRdist is computed 
between 4500 TCRs in the search DataFrame and thousands of TCRs in the 
bulk repertoire. This computation is spread across multiple cores and 
only distances within the TCRdist specified by the radius argument are 
retained to avoid a massive memory demand. 

.. note:: 
   
    The chunk_size is the 
    number of rows to compute on any one node at a time. The helper function
    get_safe_chunk will in most instances choose a chunk size that will avoid
    memory demand in excess of 1-2GB per core. By using multiple cores 
    this relatively large computational task can be completed quickly 
    and without the massive memory demand that would be required to 
    store a massive number of pairwise distances.

.. code-block:: python

    from tcrdist.breadth import get_safe_chunk
    chunk_size = get_safe_chunk(tr.clone_df.shape[0], tr_bulk.clone_df.shape[0])
    tr.compute_sparse_rect_distances(
        df = tr.clone_df, 
        df2 = tr_bulk.clone_df, 
        radius = 36,
        chunk_size = chunk_size)
        

Sequences in the meta-clonotype search set are joined by distance
with similar sequences in the bulk repertoire. Similar to database joins on a common key, 
we can perform either a "left","inner", "outer" join operation (see :py:func:`tcrdist.join.join_by_dist` for details)

.. code-block:: python
  
    from tcrdist.join import join_by_dist
    df_join = join_by_dist(
        how = 'inner',
        csrmat = tr.rw_beta,
        left_df = tr.clone_df,
        right_df = tr_bulk.clone_df,
        left_cols  = tr.clone_df.columns.to_list(),
        right_cols = tr_bulk.clone_df.columns.to_list(),
        left_suffix = '_search',
        right_suffix = '_bulk',
        max_n= 1000,
        radius = 36)
        

A a dataframe <df_join> is returned, where
for every search sequences in the LEFT hand side
the <max_n> closest neighbors sequences in the RIGHT 
hand side dataframe within the <radius>
argument are returned.

.. warning::
    Note that this includes potentially redundant hits, since a TCR clone in the 
    bulk (RIGHT) dataframe could be a neighbor of multiple 
    search sequences in the meta-clonotype (LEFT) dataframe. We 
    will be able to address this redundancy in the next few steps. 
    
First, we create new columns recording which sequence pairs meet certain meta-clonotype definitions:
    
* We apply a value of true for each sequence that is within its meta-clonotype search radius.   
* We use regex matching to determine if the each joined sequence pair matches the motif-constraint, assigning True to a column MOTIF to indicate a motif-conformant CDR3.
* We apply a value of True to the column RADIUSANDMOTIF to indicate that both the RADIUS and MOTIF condition is met for a given joined sequence pair.

.. code-block:: python

    df_join['RADIUS'] = df_join.apply(lambda x: x['radius_search'] <= x['dist'], axis = 1)
    import re
    df_join['MOTIF'] = df_join.apply(lambda x: re.search(string = x['cdr3_b_aa_bulk'], pattern = x['regex_search']) is not None, axis = 1)
    df_join['RADIUSANDMOTIF'] =  df_join['RADIUS'] & df_join['MOTIF']
    df_join['unique_clones'] = 1


.. tip::
    
   Finding number of unique clones (breadth) can be achieved by adding a column of ones for unique clones 
   that when summed will represent breadth of unique clones per ORF


Suppose you are interested in all RADIUS + MOTIF conforming sequences. 

* First, we query those joined pairs and sort them in descending order by distance.
* Then, because it is possible that a sequence could show up as a match to multiple search meta-clonotypes, we use the Pandas groupby and head  methods with the 'rank' variable. This ensures that a given bulk clone can only appear once in the final result even if it matches two meta-clnootypes.

.. code-block:: python

    df_join.query('RADIUSANDMOTIF').\
        sort_values('dist', ascending = True).\
        groupby(['rank_bulk']).\
        head(1)

Now suppose we are interested tabulating the number 
of unique RADIUS+MOTIF conformant 
clones per ORF. Note that 'count_bulk'
has stored the productive frequency 
whereas templates_bulk stores raw read count 
that might be useful for counts regression 
modeling.

.. code-block:: python

    df_join.query('RADIUSANDMOTIF').\
        sort_values('dist', ascending = True).\
        groupby(['rank_bulk']).\
        head(1).\
        groupby('protein_search').\
        sum()[['count_bulk', 'templates_bulk']]
    

.. code-block:: 

                    count_bulk  templates_bulk  unique_clones
    protein_search
    E                 0.000013               1              1
    M                 0.000401              31             23
    N                 0.001371             106             48
    ORF1ab            0.002368             183            154
    ORF3a             0.000725              56             42
    ORF6              0.000013               1              1
    ORF7a             0.000039               3              2
    ORF7b             0.000013               1              1
    ORF8              0.000026               2              2
    S                 0.001203              93             42


Alternatively, now suppose we are interested tabulating the number 
of unique RADIUS+MOTIF conformant 
clones per epitope in the spike protein. 
Note that in the final groupby method
we use 'protein_coordinate_search' instead 
of 'protein_search'


.. code-block:: python

    df_join.query('RADIUSANDMOTIF').\
        sort_values('dist', ascending = True).\
        groupby(['rank_bulk']).\
        head(1).\
        query('protein_search == "S"').\
        groupby('protein_coordinate_search').\
        sum()[['count_bulk', 'templates_bulk']]

.. code-block::

                               count_bulk  templates_bulk
    protein_coordinate_search
    S 1056:1069                  0.000091               7
    S 265:278                    0.000789              61
    S 55:69                      0.000259              20
    S 83:98                      0.000013               1
    S 860:875                    0.000052               4





:py:func:`tcrdist.join.join_by_dist` can be used in multiple contexts
where it may be useful for finding non-exact matches between sets of TCRs. 


.. _full_example:

Full Example
++++++++++++

.. code-block:: python
    :emphasize-lines: 47

    import multiprocessing
    import numpy as np
    import os 
    import pandas as pd
    from tcrdist.setup_tests import download_and_extract_zip_file
    from tcrdist.repertoire import TCRrep
    from tcrdist.breadth import get_safe_chunk
    from tcrdist.join import join_by_dist
    import re
    
    ncpus = min(multiprocessing.cpu_count(), 6)
    files = ['1588BW_20200417_PBMC_unsorted_cc1000000_ImmunRACE_050820_008_gDNA_TCRB.tsv.tcrdist3.tsv']
    if not np.all([os.path.isfile(f) for f in files]):
        download_and_extract_zip_file("ImmunoSeq_MIRA_matched_tcrdist3_ready_2_files.zip", source = "dropbox", dest = ".")
    if not os.path.isfile("bioRxiv_v2_metaclonotypes.tsv.zip"):
        download_and_extract_zip_file('bioRxiv_v2_metaclonotypes.tsv.zip', source = "dropbox", dest = ".")
    
    df_search = pd.read_csv("bioRxiv_v2_metaclonotypes.tsv", sep = "\t")
    f = '1588BW_20200417_PBMC_unsorted_cc1000000_ImmunRACE_050820_008_gDNA_TCRB.tsv.tcrdist3.tsv'
    df_bulk = pd.read_csv(f, sep = "\t")
    # When one want to track each clone indivually regardless of identical TRBV-CDR3-TRBJ
    df_bulk = df_bulk.sort_values('count').reset_index(drop = True)
    
    df_bulk['rank'] = df_bulk.index.to_list()
    
    from tcrdist.repertoire import TCRrep
    tr = TCRrep(
        cell_df = df_search, 
        organism = "human",
        chains = ['beta'], 
        compute_distances= False)
    tr.cpus = ncpus
    
    tr_bulk = TCRrep(
        cell_df = df_bulk, 
        organism = "human",
        chains = ['beta'], 
        compute_distances= False)  
    
    chunk_size = get_safe_chunk(tr.clone_df.shape[0], tr_bulk.clone_df.shape[0])
    tr.compute_sparse_rect_distances(
        df = tr.clone_df, 
        df2 = tr_bulk.clone_df, 
        radius = 36,
        chunk_size = chunk_size)
    
    df_join = join_by_dist(
        how = 'inner',
        csrmat = tr.rw_beta,
        left_df = tr.clone_df,
        right_df = tr_bulk.clone_df,
        left_cols  = tr.clone_df.columns.to_list(),
        right_cols = tr_bulk.clone_df.columns.to_list(),
        left_suffix = '_search',
        right_suffix = '_bulk',
        max_n= 1000,
        radius = 36)
    
    df_join['RADIUS'] = df_join.apply(lambda x: x['radius_search'] <= x['dist'], axis = 1)
    import re
    df_join['MOTIF'] = df_join.apply(lambda x: re.search(string = x['cdr3_b_aa_bulk'],
        pattern = x['regex_search']) is not None, axis = 1)
    
    df_join['RADIUSANDMOTIF'] =  df_join['RADIUS'] & df_join['MOTIF']
    df_join['unique_clones'] = 1
    
    df_join.query('RADIUSANDMOTIF').\
        sort_values('dist', ascending = True).\
        groupby(['rank_bulk']).\
        head(1).\
        groupby('protein_search').\
        sum()[['count_bulk', 'templates_bulk','unique_clones']]

    """
                    count_bulk  templates_bulk  unique_clones
    protein_search
    E                 0.000013               1              1
    M                 0.000401              31             23
    N                 0.001371             106             48
    ORF1ab            0.002368             183            154
    ORF3a             0.000725              56             42
    ORF6              0.000013               1              1
    ORF7a             0.000039               3              2
    ORF7b             0.000013               1              1
    ORF8              0.000026               2              2
    S                 0.001203              93             42
    """



.. automodule:: tcrdist.join

.. autofunction:: tcrdist.join.join_by_dist