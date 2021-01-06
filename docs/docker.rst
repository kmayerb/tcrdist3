.. _docker:

Docker
======


The fastest way to try out tcrdist3 is to use our Docker image. 

If you run Linux, we assume you know how to use Docker. If you 
use a Mac or Windows machine, you will need to download 
`Docker Desktop <https://www.docker.com/products/docker-desktop>`_ 
in order to run a Docker container. 

Getting the Docker Image
------------------------

Pull the image:

.. code-block:: none
    
    docker pull quay.io/kmayerb/tcrdist3:0.1.9


Running a Container
-------------------

Run a container iteractively (flag `-it`):


.. code-block:: none
    
    docker run -it quay.io/kmayerb/tcrdist3:0.1.9


The tcrdist3 image runs Python 3.8 and IPython, so the fastest way to try things out
is to run IPython and copy and paste the code snippets on this doc page. It's easy to copy 
a file from inside a running Docker container to your host machine, but remember once you 
exit a container your files don't persist. 

.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube.com/embed/JBZtF5sc8QE" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    


Mounting Files
--------------

.. tip::
    Mounting your docker container to a working director is easy. 
    (For more on this see `Docker for the busy researcher <http://erick.matsen.org/2018/04/19/docker.html>`_)
    
.. code-block:: none

    docker run -v ${HOME}/mydata/:/yourdata/ -it quay.io/kmayerb/tcrdist3:0.1.9
    ls data 
 
    
.. raw:: html

	<iframe width="672" height="378" src="https://www.youtube.com/embed/78m5WXr895w" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
    
    
Copying Files
-------------

You can also copy a file out of running container if you forget to mount a shared folder.

.. code-block:: none
    
    """
    find the container ID
    """
    docker ps 
    """
    CONTAINER ID        IMAGE                            COMMAND
    3a51162eb585        quay.io/kmayerb/tcrdist3:0.1.9   "/bin/bash"
    """
    docker cp 3a51162eb585:hierdiff_example_PA_v_PB1.html hierdiff_example_PA_v_PB1.html
