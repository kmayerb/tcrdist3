.. _sparsity:

Sparse Representation
=====================

For large datasets, you may want to set compute_distances to False and then use a sparse implementation. First, set `tr.cpus` to the appropriate number of cpus available. When computing distances with the sparse implementation, the argument radius is the maximum distance to be stored. All distances greater than `radius` will be converted to zero, reducing the memory required. The argument `chunk_size` tells tcrdist3 how many rows to compute at a time. For instance, if you have 100,000 x 100,000 clones, then a chunk size of 100 will compute distances 100x100,000 on each node and save each of the 1000 intermediate results in sparse format before recombining them it a single sparse matrix. Larger chunk sizes will result in less overhead, but chunk size should be tuned based on available memory. The results are object attributes `rw_beta` and `rw_alpha`, which store as scipy.sparse.csr matrices all distances less than radius. True zeros are represented as -1. The techniques for customizing the distance metric such as changing trims, gap when using the sparse implementation.

.. literalinclude:: ../docs/sparse.py
    :language: python

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube-nocookie.com/embed/3G39TcD361w" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

