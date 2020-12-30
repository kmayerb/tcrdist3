.. _weclome:

Video Welcome
=============

tcrdist3 can be run interatively. These video tutorials offer a quick way to test your installation or running docker container and get started.


Hello tcrdist3 
--------------

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube.com/embed/podmo6F5_-4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

This is a hello world video for tcrdist3. It's an easy test to make sure your installation is working and 3,686,400 pairwise paired chain distances in less than a few seconds isn't a bad place to start.

Sparse Representation
---------------------

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube-nocookie.com/embed/3G39TcD361w" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


For large datasets, you may want to set compute_distances to False and then use a sparse implementation. First, set `tr.cpus` to the appropriate number of cpus available. When computing distances with the sparse implementation, the argument radius is the maximum distance to be stored. All distances greater than `radius` will be converted to zero, reducing the memory required. The argument `chunk_size` tells tcrdist3 how many rows to compute at a time. For instance, if you have 100,000 x 100,000 clones, then a chunk size of 100 will compute distances 100x100,000 on each node and save each of the 1000 intermediate results in sparse format before recombining them it a single sparse matrix. Larger chunk sizes will result in less overhead, but chunk size should be tuned based on available memory. The results are object attributes `rw_beta` and `rw_alpha`, which store as scipy.sparse.csr matrices all distances less than radius. True zeros are represented as -1.