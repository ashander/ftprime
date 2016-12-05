ftprime
======

Important files:

-  [devel/forwards_algorithm.md](devel/forwards_algorithm.md): Description of the algorithm for outputting a valid tree sequence from a forwards-time simulation.
-  [devel/pedrecorder.py](devel/pedrecorder.py): Implementation of the forwards algorithm.
-  [devel/meiosis.py](devel/meiosis.py): Wraps the storage mechanism in pedrecorder.py into a tagger that can be used in simuPOP.
-  [devel/with-simupop.py](devel/with-simupop.py): Example of using the simuPOP interface.
-  [devel/wf.py](devel/wf.py): Very simple forwards-time Wright-Fisher simulation that uses the pedigree recorder.
-  [devel/algorithm-notes.md](devel/algorithm-notes.md): A writeup of algorithmic considerations.
-  [devel/test_msprime.py](devel/test_msprime.py): An exploration of the capabilities of msprime to digest various notions of a tree sequence.



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


Notes
-----

Since we aren't tracking genotypes explicitly,
we can take our individuals to be diploid without even telling simuPOP about it.
