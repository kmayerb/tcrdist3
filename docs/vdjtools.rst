.. _vdjtools:

VDJtools Data
=============

An example of importing from VDJtools formated input. 

Import VDJtools formated input ``.tsv`` or ``.tsv.gz file``. 

.. literalinclude:: ../tcrdist/tests/test_vdj_funcs.py
    :linenos:
    :lines: 6-27
    :dedent: 4
    :language: python
 

.. code:: none

	count      freq        cdr3_b_aa     v_b_gene    j_b_gene                                  cdr3_b_nucseq  valid_v  valid_j  valid_cdr3
	38957  0.122460   CASSPRAGKGEQFF   TRBV7-9*01  TRBJ2-1*01     TGTGCCAGCAGCCCAAGGGCAGGGAAGGGTGAGCAGTTCTTC     True     True        True
	24130  0.075852    CASSFWTPYEQYF  TRBV12-4*01  TRBJ2-7*01        TGTGCCAGCAGTTTTTGGACACCCTACGAGCAGTACTTC     True     True        True
	19711  0.061961   CASSAPAGVGEQYF   TRBV7-9*01  TRBJ2-7*01     TGTGCCAGCAGCGCCCCAGCGGGGGTCGGCGAGCAGTACTTC     True     True        True
	19285  0.060622   CASSPPGQHNEQFF   TRBV5-6*01  TRBJ2-1*01     TGTGCCAGCAGCCCGCCGGGACAGCACAATGAGCAGTTCTTC     True     True        True
	 8710  0.027380  CASSLDEPTDNEQFF   TRBV5-1*01  TRBJ2-1*01  TGCGCCAGCAGCTTGGACGAGCCTACCGACAATGAGCAGTTCTTC     True     True        True



Input file must have columns ``['count', 'freq', 'cdr3aa', 'v', 'j','cdr3nt']``. 

This is a powerful place to start because there are many tools for converting 
outputs of popular tools to this format. 
See section on 'Formats supported for conversion' in the 
`VDJtools DOCS <https://vdjtools-doc.readthedocs.io/en/master/input.html#vdjtools-format>`_ 
for more details on conversion to this format from 

* MiTCR,
* MiGEC,
* IgBlast (MIGMAP), 
* ImmunoSEQ, 
* VDJdb, 
* Vidjil, and
* MiXCR
