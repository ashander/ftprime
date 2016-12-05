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

**Diploidy:** Chromosomes are labeled via a 1-to-2 map from individual IDs,
so that IDs for msprime are the floor of the ID from simuPOP divided by two.

**Division of labor:** We've split the work out into two parts: 

1. one part that interfaces with simuPOP
to assign IDs, choose recombination breakpoints, and say which chromsome is which other's parent;

2. and the second part (pedrecorder.py) the push all this into msprime's coalescence records.

**The clock:** I haven't figured out how to get simuPOP to say what time it is in during-reproduction operators,
so we've just got an external clock that we have to update (not a big deal).

To-do
=====

1. Change recombination to be Poisson, and **check** we are not missing a factor of 2 in the definition of Morgans.
