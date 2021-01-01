.. _downloads:

Downloads
=========


.. automodule:: tcrdist.setup_tests

All of the code snippets shown in these docs are part of tcrdist3's 
suite of integration tests. If you want to re-run any of these 
snippets with the same input data files:

* Get the data with :py:func:`tcrdist.setup_tests.download_and_extract_zip_file`.
* See available zip files using:  :py:func:`tcrdist.setup_tests.list_available_zip_files`.

.. literalinclude:: ../tcrdist/tests/test_downloads.py
    :linenos:
    :lines: 6-50
    :dedent: 4
    :language: python



.. autofunction:: list_available_zip_files

.. autofunction:: download_and_extract_zip_file