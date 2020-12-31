.. _welcome:

Hello tcrdist3 
==============

tcrdist3 can be run interatively. The best way to learn more is 
to try a few examples. Examples use the data file (`dash.csv <https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/tcrdist/test_files_compact/dash.csv>`_).


Hello tcrdist3 
--------------

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube.com/embed/podmo6F5_-4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

This is a hello world video for tcrdist3. It's an easy test to make sure your installation is working and 3,686,400 pairwise paired chain distances in less than a few seconds isn't a bad place to start.

Try it:

.. literalinclude:: ../tcrdist/tests/test_example_1.py
    :linenos:
    :lines: 6-28
    :dedent: 4
    :language: python

For details on options on computing distances, see :ref:`tcrdistances`.


Sparse Representation
---------------------

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube-nocookie.com/embed/3G39TcD361w" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


For large datasets, you may want to set compute_distances to False and then use a sparse implementation. First, set `tr.cpus` to the appropriate number of cpus available. When computing distances with the sparse implementation, the argument radius is the maximum distance to be stored. All distances greater than `radius` will be converted to zero, reducing the memory required. The argument `chunk_size` tells tcrdist3 how many rows to compute at a time. For instance, if you have 100,000 x 100,000 clones, then a chunk size of 100 will compute distances 100x100,000 on each node and save each of the 1000 intermediate results in sparse format before recombining them it a single sparse matrix. Larger chunk sizes will result in less overhead, but chunk size should be tuned based on available memory. The results are object attributes `rw_beta` and `rw_alpha`, which store as scipy.sparse.csr matrices all distances less than radius. True zeros are represented as -1.


Try it: 

.. literalinclude:: ../docs/sparse.py
    :linenos:
    :lines: 1-28
    :language: python


