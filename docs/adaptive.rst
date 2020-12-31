.. _adaptive:

Adaptive ImmunoSEQ Data
=======================

Adaptive uses a distinct naming convention to `IMGT Nomenclature <http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.html>`_. This poses a formatting challenge when using ImmunoSEQ files as inputs to tcrdist3. According to Adaptive's technical team: "Adaptive's nomenclature is more expanded to both facilitate alphanumeric sorting, and also specify the precision of the identification." Adaptive's technical team further explained to us the difference between naming systems, which we paraphrase here: 

Both naming systems follow a [locus and family]-[gene]*[allele] convention, where IMGT naming prioritizes brevity, opting for "a single letter or number where possible" (except for alleles). IMGT also leaves out gene-level information when there is only one gene in the family. For instance, IMGT drops the gene-level info in naming TRBV15*02. By contrast, Adaptive uses the following three possible names:

* A gene with allele-level identification: TCRBV15-01*02
* Gene-level identification: TCRBV15-01
* Family-level only: TCRBV15

Adaptive's output files can contain gene-level names within the 'bioidentity' field like TCRBV15-X, when there is ambiguity about the gene-level assignment. 

tcrdist3 uses IMGT gene names throughout, so the first step to working with ImmunoSEQ files is name conversion. To avoid losing lots of CDR3 data, when the V gene may not be fully resolved we often use Adaptive 'bioidentity' gene-level calls and replace allele with *01. Depending on your project's goals, you may want to do this cleaning by hand, so let's first take a look at how to convert Adaptive's `v_gene` into its IMGT*01 equivalent:

.. code-block:: python
    :emphasize-lines: 5

    !wget https://raw.githubusercontent.com/kmayerb/tcrdist3/master/Adaptive2020.tsv
    import pandas as pd
    from tcrdist.swap_gene_name import adaptive_to_imgt
    adpt_input = pd.read_csv('Adaptive2020.tsv', sep = '\t')
    adpt_input['v_b_gene'] = adpt_input['v_gene'].apply(lambda x : adaptive_to_imgt['human'].get(x))



.. _loading_adaptive_biotechnology:

Cleaning Adaptive ImmunoSEQ Files
---------------------------------

We also have a one line conversion function that works with recent ImmunoSEQ files
containing the 'bioidentity' field, as shown here:

import_adaptive_file
++++++++++++++++++++

.. literalinclude:: ../tcrdist/tests/test_introduction_2.py
    :linenos:
    :lines: 6-10
    :dedent: 4
    :language: python


.. automodule:: tcrdist.adpt_funcs

.. autofunction:: import_adaptive_file


After conversion, the data as a Pandas DataFrame can be directly imported to tcrdist3.

Loading Adaptive ImmunoSEQ Files
------------------------------------

.. literalinclude:: ../tcrdist/tests/test_introduction_2.py
    :linenos:
    :lines: 6-16
    :dedent: 4
    :language: python

Look Up Adaptive Conversion
---------------------------

.. literalinclude:: ../tcrdist/tests/test_introduction_2.py
    :linenos:
    :lines: 17-50
    :dedent: 4
    :language: python
