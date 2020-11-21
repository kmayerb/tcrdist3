.. _radius:


Radius
======

tcrdist3 uses a flexible distance based metric to find biochemically similar, non identical, TCR chains.
**What distance is appropriate?** A conservative approach is to use a small fixed radius such as 16 
TCRdistance units (TDUs), this corresponds with 1-3 amino acid changes between CDR3s, depending on 
biochemical similarity of the substituted residues. 

Rather than use a fixed radius for finding similar sequences, on can be optimized
based on how common a given query TCR is within an un-enriched bulk
background repertoire. Due to biases in VDJ recombination some TCRs 
will have more neighbors, whereas other longer TCRs with more 
randomly inserted nucleotides will naturally have fewer neighbors, absent 
strong convergent selective pressure, such as from immune maturation, vaccination, infection, or cancer.

tcrdist3 allows a user to pick a radius based on an estimated rate of 
TCR-neighbor discovery in a control repertoire of your choosing. 

The following example shows how:

.. literalinclude:: ../tcrdist/tests/test_example_24_centers.py
    :linenos:
    :lines: 6-115
    :dedent: 4
    :language: python


|pic1| |pic2|

.. |pic1| image:: _static/PA1.png  
   :width: 49%

.. |pic2| image:: _static/PA2.png
   :width: 49%


See `PA QuasiPublic Meta Clonotypes <_static/quasi_public_clones_PA.html>`_.









