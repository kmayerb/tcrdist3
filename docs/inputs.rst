.. _Inputs:

Inputs
======

Data
++++

The tcrdist3 standard input is a Pandas DataFrame.

The header and first line of a typical input for a beta-chain analysis would look like this:

+------------+------------+------------+------------+------------+----------------------+-----------------------------------------------------------+
| subject    | epitope    | count      | v_b_gene   | j_b_gene   | cdr3_b_aa            | cdr3_b_nucseq                                             |
+============+============+============+============+============+======================+===========================================================+
| s1         |   NP       |   1        | TRBV1*01   | TRBJ1-1*01 | CACDSLGDKSSWDTRQMFF  | TGTGCCTGTGACTCGCTGGGGGATAAGAGCTCCTGGGACACCCGACAGATGTTTTTC |
+------------+------------+------------+------------+------------+----------------------+-----------------------------------------------------------+			


Column names reflect the chain under investigation. 
- a : alpha
- b : beta
- g : gamma
- d : delta

One or more of the following columns are required and are case-sensitive  
    - 'v_a_gene', 'v_b_gene', 'v_g_gene', or 'v_d_gene' 
    - 'j_a_gene', 'j_b_gene', 'j_g_gene', or 'j_d_gene' 
    - 'cdr3_a_aa', 'cdr3_b_aa', 'cdr3_g_aa', or  'cdr3a_d_aa'
    - 'cdr3_a_nucseq', 'cdr3_b_nucseq, 'cdr3_g_nucseq', or 'cdr3a_d_nucseq' 


For v_x_gene, include the full IMGT gene name and allele (e.g., TRBV1*01). If you don’t know the allele, use *01. But an allele must be present to infer v_b_genes based on matching one of the id rows in this `table <https://github.com/kmayerb/tcrdist3>`_.


.. tip::

    Two of each can be supplied for paired analysis. tcrdistances can be calculated 
    without nucleotide sequences, but some other features require them.


The following is required.
    - 'count'

The following are optional:
    - 'epitope`
    - 'subject'


The following are usually inferred from germline reference v-gene but can be supplied by the user in some advanced use-cases only!
    -  'cdr1_a_aa', 'cdr1_b_aa',  'cdr1_g_aa',  or 'cdr1_d_aa'
    -  'cdr2_a_aa', 'cdr2_b_aa',  'cdr2_g_aa',  or 'cdr2_d_aa'
    -  'pmhc_a_aa', 'pmhc_a_aa',  'pmhc_a_aa',  or 'pmhc_a_aa' (pmhc = cdr 2.5)

.. tip::

  CDR2.5, the pMHC-facing loop between CDR2 and CDR3, are referred to in tcrdist3 as pmhc_a and phmc_b, respectively.




Arguments
+++++++++

chain(s)
-------

Most classes and functions in tcrdist3 require specification of the appropriate t cell receptor 
chains:
    - ['alpha'], ['beta'], ['gamma'], or ['delta'] for single-chain analysis, 
    - ['alpha', 'beta'] or ['gamma', 'delta'] for paired-chain analyis 

organism
--------

Most classes and functions in tcrdist3 require specification of an appropriate host organism. 
Currently only 'human' or 'mouse' are supported. This is required because reference TCR genes
are organism specific. 

db_file
-------

The `db_file` is used by tcrdist3 to supply updated information about reference TCR germline sequences. 


.. tip:: 

    Getting new database files:
    Reference json  https://github.com/repseqio/library-imgt/releases
    `Data coming from IMGT server may be used for academic research only, provided that it is referred to IMGT®, and cited as "IMGT®, the international ImMunoGeneTics information system® http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France)."`
