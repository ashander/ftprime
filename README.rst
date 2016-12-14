ftprime
======

Contents
--------

This package does two things: provides code to simulate the outcome of mating of diploid individuals, in terms of stretches of chromosome inherited;
and to keep track of these in a forwards-time simulation in a way appropriate for input into [msprime](https://github.com/jeromekelleher/msprime).
This also provides a class to facilitate doing this with [simuPOP](https://github.com/BoPeng/simuPOP).

Since from the point of view of the ARG it is natural for the unit of populations to be *chromsomes* rather than (diploid) individuals,
we deal with two types of IDs: *individual* and *chromosome*; the population simulation software deals with individuals;
while the interface to msprime deals with chromosome IDs, **but** as it doesn't know about higher-level individuals,
refers to these as "individuals".

[Core module:](ftprime/)

-  [ftprime/argrecorder.py](ftprime/argrecorder.py): Provides `ARGrecorder`, really just an ordered dict whose keys are chromosome IDs
    and whose values are ordered lists of nonoverlapping coalescence records.  This must be initialized using `add_individual`; and
    kept up to date using `add_record`; at the end samples are designated using `add_samples` and the msprime tree sequence is output 
    with `tree_sequence`.  Note that this needs to be given upper bounds on the total number of generations and the total number of samples
    at the start.  The work of adding a newly inherited segment to the list of coalescence records is done by `merge_records`.

-  [ftprime/meiosis.py](ftprime/meiosis.py): Provides `MeiosisTagger`, which can be used as an IdTagger in simuPOP
    with the side effect of simulating recombination events and storing everything in an `ARGrecorder`.

Tests and examples:

-  [tests/test_ftprime_with_simuPOP.py](tests/test_ftprime_with_simuPOP.py): Example of using the simuPOP interface.
-  [tests/wf/](test/wf/__init__.py): Very simple forwards-time Wright-Fisher simulation that uses the underlying machinery to the ARGrecorder.
-  [tests/test_merge_records_with_wf.py](tests/test_merge_records_with_simuPOP.py): Example of using the wf interface.

[Documentation of the problem and the methods:](writeups/)

Since there are many different ways to store an ARG as a set of coalescence records,
a good deal of this is devoted to describing and verifying msprime's requirements
for such a set, and thinking about different ways to do it.

-  [writeups/forwards_algorithm.md](writeups/forwards_algorithm.md): Description of the algorithm for outputting a valid tree sequence from a forwards-time simulation.
-  [writeups/algorithm-notes.md](writeups/algorithm-notes.md): A writeup of algorithmic considerations.
-  [devel/test_msprime.py](devel/test_msprime.py): An exploration of the capabilities of msprime to digest various notions of a tree sequence.



Development
-----------


Install in locally editable (``-e``) mode and run the tests:

.. code-block:: console

    pip install -e .[test]
    pytest


To-do
=====

1. Change recombination to be Poisson, and **check** we are not missing a factor of 2 in the definition of Morgans.
