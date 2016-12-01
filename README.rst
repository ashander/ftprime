ftprime
======

Important files:

-  [devel/algorithm-notes.md]: The writeup of algorithmic considerations.
-  [devel/test_msprime.py]: An exploration of the capabilities of msprime to digest various notions of a tree sequence.


Development
-----------


Install in locally editable (``-e``) mode and run the tests:

.. code-block:: console

    pip install -e .[test]
    pytest

    
Style
----

We'll try to use PEP8, which you can enforce using these vim plugins:

.. code-block:: console

    Plugin'klen/python-mode'
    Plugin 'scrooloose/syntastic'
