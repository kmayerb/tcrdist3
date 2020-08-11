.. _tcrdistances:

TCR Distances
=============

The following short examples use data `dash.csv <https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/tcrdist/test_files_compact/dash.csv>`_ (375KB). 

I am happy to use the defaults
------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_1.py
    :linenos:
    :lines: 6-28
    :dedent: 4
    :language: python

I'd like to tweak a default parameter
-------------------------------------

.. literalinclude:: ../tcrdist/tests/test_example_2.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python

I want complete control
-----------------------

.. literalinclude:: ../tcrdist/tests/test_example_3.py
    :linenos:
    :lines: 6-120
    :dedent: 4
    :language: python

.. tip::

    * metrics - tcrdist3 allows specification of a unique metric to each CDR. You can even provide your own
        * By default -  tcrdist3 uses pwseqdist >= 0.2.0, by the same authors, for various distance metrics
            * `pw.metrics.nb_vector_tcrdist` - a Numba accelerated metric that matches the approach in Dash et al. (2018)
            * `pw.metrics.nw_hamming_metric` - a pairwise Needleman-Wunsch alignment followed by hamming distance (with multiprocessing enabled for faster run-time)
            * `pw.metrics.nw_metric` -  a reciprocal Needleman-Wunsch alignment score based dissimilarity using BLOSUM62
        * The field - commonly uses editdistance to compare CDRs. For instance, the Cython enabled packages like python-Levenshtein metrics can be used seamlessly with tcrdist3
        * Dictionary keys must match appropriate columns of the clone_df DataFrame.

.. tip::

    * weights - typically the CDR3 is given a higher weight because it is in direct contact with antigen peptide.
        * Dictionary keys must match appropriate columns of the clone_df DataFrame.

.. tip::

    * kargs -  tuneable parameters
        * Tunable parameters are passed to metric as keyword-arguments. Some important parameters for  `pw.metrics.nb_vector_tcrdist` are
            * `use_numba` - must be set to true for metrics starting with 'nb'.
            * `ctrim` - number of amino acids to trim of c-terminal end of a CDR. (i.e. ctrim = 3 CASSQDFEQ-YF consider only positions after CAS-SQDFEQ-YF)
            * `ntrim` - number of amino acids to trim of n-terminal end of a CDR. (i.e. ntrim = 2 CASSQDFEQ-YF consider only positions before "YF", SQDFEQ-YF)
            * `gap_penalty` - the penalty accrued for gaps in the pairwise alignment.
            * 'fixed_gappos' - When sequences are of different length, there will be a gap position. If 'fixed_gappos' is False, then the metric inserts a single gap at an optimal position based on a BLOSUM62 scoring matrix. This is recommended for the CDR3, but it is not necessary when the CDR1, 2, and 2.5 are already imgt_aligned and of a fixed length.
         * The default parameters match those used in Dash et al. (2018) *Quantifiable predictive features define epitope-specific T cell receptor repertoires
         * Dictionary keys must match appropriate columns of the clone_df DataFrame.


.. tip::

  Protip: CDR2.5 alpha and CDR2.5 beta, the pMHC-facing loop between CDR2
  and CDR3, are referred to in tcrdist2 as pmhc_a and phmc_b, respectively.




I just want to count mismatches
-------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_4.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python

I want to use my own distance metric
------------------------------------


.. literalinclude:: ../tcrdist/tests/test_example_5.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python

I want tcrdistances, but I hate OOP
-----------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_7.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python

I hate OOP, and I only want distances for the CDR3 
--------------------------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_6.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python

I want tcrdistances but I want to keep my variable names
--------------------------------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_8.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python


I want to use TCRrep but I want to keep my variable names
---------------------------------------------------------
.. literalinclude:: ../tcrdist/tests/test_example_9.py
    :linenos:
    :lines: 6-100
    :dedent: 4
    :language: python


I want something else
---------------------
Let us know. We are happy to help. Get in touch with *kmayerbl AT fredhutch DOT org*