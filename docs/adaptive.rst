.. _adaptive:

Adaptive ImmunoSEQ Data
=======================

The primary challenge in using ImmunoSEQ files is the different naming convention for TRV and TRJ genes. Adatptive's technical team explained the difference clearly, which I will paraphrase here: Adaptive uses a distinct naming convention to `IMGT Nomenclature <http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.html>`_. Both follow a [locus and family]-[gene]*[allele] convention, where IMGT naming prioritizes brevity, opting for "a single letter or number where possible" (except for alleles). IMGT also leaves out gene-level information when there is only one gene in the family. For instance, the gene-level info is dropped in naming TRBV15*02. According to Adaptive's technical team: "Adaptive's nomenclature is more expanded to both facilitate alphanumeric sorting, and also specify the precision of the identification."

* A gene with allele-level identification: TRBV15-01*01
* Gene-level identification: TRBV15-01
* Family-level only: TRBV15

When Gene-level resolution is missing, we have found a that some of Adaptive's output files can contain gene-level names within the bioidentiy field like TRBV15-X, when there is ambiguity about the gene-level assignment.

tcrdist3 uses IMGT gene names throughout, so the first step to working with ImmunoSEQ files is name conversion.  To avoid losing lots of CDR3 data, when V gene may not be full resolved we often use Adaptive gene-level calls and replace allele with *01. You may want to do this cleaning by hand and include genes resolved at the allele level, so let's take a look at how to convert Adaptive's `v_gene` into it's IMGT*01 equivalent:

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

import_adaptive_file
++++++++++++++++++++

.. literalinclude:: ../tcrdist/tests/test_introduction_2.py
    :linenos:
    :lines: 6-10
    :dedent: 4
    :language: python


.. automodule:: tcrdist.adpt_funcs

.. autofunction:: import_adaptive_file



Loading Adaptive Biotechnology Files
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
