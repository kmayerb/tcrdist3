.. _tree:

Trees
=====

Interactive tree diagrams can be easily produced in tcrdist3. 
To automate the processes decsribed in more detail on the :ref:`motif_gallery` page, 
initiate a `TCRtree` class as shown below. The result is a 
html page with an `Interactive Hierdiff Tree <_static/dash.mouse.b.tree.html>`_.
Hovering on the nodes reveals information about each node including a sequence logo.


I am happy to use the defaults
------------------------------

.. literalinclude:: ../tcrdist/tests/test_example_13.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python


I'd like to tweak a default parameter
-------------------------------------

There are three core processes executed by TCRtree.build_tree():
* hcluster_diff - hierarchically clusters TCRs, and tallies how many TCRs in each node of the hiearchical tree
have particular catagorical label values (i.e. number that are CD8+ vs CD4+, or number coming from a pre-vaccine vs. 
post-vaccine sample).

* member_summ - summarizes cluster meta data such as the % of TCRs with a given V gene
* plot_hclust (part of hierdiff package) - make the D3 interactive tree

When invoked via `.build_tree()`, each of these functions is controlled by keyword arguments ('kwargs')
dictionaries stored as TCRree default attributes. 
Attribute values for each can be found in the docstrings for `?TCRtree`. 
Here we mention some of the most important. First 'x_cols' passed to 
`hcluster_diff`


.. literalinclude:: ../tcrdist/tests/test_example_14.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python


I want more control
-------------------

If you want full control, you can pass an entire dictionary to TCRtree.default attributes
prior to running `build_tree().` For instance remove any of 'tootips_cols' to simplify 
what is displayed when one hovers over a tree node.

.. literalinclude:: ../tcrdist/tests/test_example_15.py
    :linenos:
    :lines: 6-80
    :dedent: 4
    :language: python


I want a paired alpha and beta tree 
-----------------------------------

This is accomodated. But you may have to decide which chain's distance matrix you wish to use for the purpose of clustering. 

.. literalinclude:: ../tcrdist/tests/test_example_16.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python




Changing the background
-----------------------

Perhaps you have a larger dataset and you don't need all the SVG logos. 





I prefer to write my own tree from scratch
------------------------------------------











