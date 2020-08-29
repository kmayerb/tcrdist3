.. _simulate_cdr3s:

Simulate CDRs
=============

Simulation of 1 Million CDRs for both Human and Mouse (beta only).

.. literalinclude:: ../tcrdist/tests/longtest_simulate_cdr3_w_olga.py
    :linenos:
    :lines: 1-100
    :language: python

These files can be downloaded directly using:

.. code-block:: python

    from tcrdist.setup_tests import download_and_extract_zip_file
    download_and_extract_zip_file('olga_T_alpha_beta_1000K_simulated_cdr3.zip')

Simulation of CDRs rely on the OLGA by Sethna and colleages (2019).

References 
----------

Zachary Sethna, Yuval Elhanati, Curtis G Callan, Aleksandra M Walczak, Thierry Mora
`Bioinformatics (2019) <https://doi.org/10.1093/bioinformatics/btz035>`_ 
OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs
