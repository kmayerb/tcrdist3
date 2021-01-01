.. _motif_gallery:

CDR3 Motifs
===========


We've built flexibility into tcrdist3's motif discovery process. This code makes the following 
`Interactive PA-PB1 Epitope Hierdiff Tree <_static/hierdiff_example_PA_v_PB1.html>`_.

.. literalinclude:: ../tcrdist/tests/test_gallery_hdiff.py
    :linenos:
    :lines: 6-200
    :dedent: 4
    :language: python


.. tip ::
    This example introduces features that are implemented in two 
    stand-alone pip installable python pakcages, by tcrdist3's authors, 
    `tcrsampler <https://pypi.org/project/tcrsampler/>`_
    and `palmotif <https://pypi.org/project/palmotif/>`_. 
    They are more extensively `documented <https://github.com/kmayerb/tcrsampler/blob/master/docs/tcrsampler.md>`_ on their own project pages, 
    but we introduce their use here with some basic illustrative examples in the 
    :ref:`modules` section of this page.


.. _modules:

Modules
-------

tcrsampler
++++++++++

Suppose you want to compare a cluster of 3 biochemically similar TCRs 
to a sample of TCRs with the same V and J gene usage from a background 
set, to detect selective pressure beyond what might arise by 
natural V(D)J recombination biases.

.. code-block:: python

    from tcrsampler.sampler import TCRsampler

    t = TCRsampler()
    #t.download_background_file("ruggiero_mouse_sampler.zip")
    t = TCRsampler(default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv')
    
    df = pd.DataFrame( { 
        "cdr3_b_aa": ['CASSPVRAGDTQYF', 'CASSPIRVGDTQYF', 'CASSPVRLGDTQYF'],
        "v_b_gene":['TRBV29*01', 'TRBV29*01','TRBV29*01'],
        "j_b_gene":['TRBJ2-7*01','TRBJ2-7*01','TRBJ2-5*01']
        })
    gene_usage = df.groupby(['v_b_gene','j_b_gene']).size()
    t.sample( gene_usage.reset_index().to_dict('split')['data'], 
        flatten = True, 
        depth = 2,
        seed = 1)
    """
    ['CASSPGHNQDTQYF',
    'CASSPGGGQDTQYF',
    'CASSPLGQGSYEQYF',
    'CASSFGQGADEQYF',
    'CASSHGEQYF',
    'CASSPGQSSYEQYF']
    """


.. tip :: 
    You can use a default background (the example above used 
    data from 'High-resolution analysis of the human T-cell receptor repertoire'
    by Ruggiero and colleagues (2015), or see tcrsampler docs 
    for info on creating your own background sampler from the
    TCR data of your choice).



palmotif
++++++++

For collections of variable length CDR3s, 
palmotif aligns, computes and plots sequence logos, 
with output as SVG or matplotlib. 
One can use the background set generated above with a 
TCRsampler to compare our cluster to a background, given 
the same frequency of V and J gene usage.

.. code-block:: python

    from palmotif import compute_pal_motif, svg_logo
    motif, stats = \
        compute_pal_motif(
            centroid = 'CASSPVRAGDTQYF',
            seqs = ['CASSPVRAGDTQYF', 'CASSPRVGDTQYF', 'CASSPVRLGDTQYF'], 
            refs =  t.sample( 
                gene_usage.reset_index().to_dict('split')['data'], 
                flatten = True, 
                depth = 10,
                seed = 1)
        )

    """return a str if return_str = True"""
    svg = svg_logo(motif, return_str = True)

    """Or directly output svg to a file"""
    svg_logo(motif, 'test.svg')

.. image:: _static/test.svg
   :target: _static/test.svg

One can also plot the raw logo without considering background, 
simply by omiting the `refs` argument:

.. code-block:: python

    from palmotif import compute_pal_motif, svg_logo
    motif, stats = \
        compute_pal_motif(
            centroid = 'CASSPVRAGDTQYF',
            seqs = ['CASSPVRAGDTQYF', 'CASSPRVGDTQYF', 'CASSPVRLGDTQYF']
        )
    """return a str if return_str = True"""
    svg = svg_logo(motif, return_str = True)
    
    """Or directly output svg to a file"""
    svg_logo(motif, 'testraw.svg')

.. image:: _static/testraw.svg
   :target: _static/testraw.svg