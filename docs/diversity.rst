.. _diversity_example:

Repertoire Diversity
====================

The diversity of a TCR repertoire can be an important indicator: loss of diversity could mean antigen driven clonal expansion or could be a sign of a successful enrichment experiment. There are many ways to compute diversity/clonality, many borrowing from the field of ecology. The Simpson's Diversity index is one of the more commonly used measures of diversity and was recently reviewed and expanded upon for community analysis by 'Grabchak et al. (2017) <https://doi.org/10.1371/journal.pone.0173305>`_. We've added code for the generalized Simpson's diversity index in tcrdist3 to allow for easy application to TCR repertoires, including numba-compiled code to quickly compute confidence intervals and compare the diversity of two repertoires with a statistical test.

We also added code to compute a distance-aware "fuzzy" diversity measure, originally termed TCRdiv in Dash et al (2017). While SDI order=2 is defined as (one minus) the probability of sampling two identical TCRs from a repertoire, the fuzzy index is (one minus) the probability of sampling two TCRs that are within a specified distance of each other. Using a fuzzy metric can help reveal subtle changes in diversity in repertoires that may have very few identical clones, but that may have many clones with a similar TCR sequence. The example below demonstrates how to use these diversity measures on TCR repertoire data.

.. literalinclude:: ../tcrdist/tests/test_diversity.py
    :linenos:
    :lines: 23-63
    :dedent: 4
    :language: python
