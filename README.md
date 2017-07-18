ftprime
======

This code was mentioned during the Evolution 2017 talk, Ashander _et al_ (2017) "Using Genomic Data to Inform Population Viability in a Long-Lived Endangered Vertebrate". Note the code is pre-production but if you're building on the ideas here please cite this repository [![DOI](https://zenodo.org/badge/72480698.svg)](https://zenodo.org/badge/latestdoi/72480698)


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

-  [ftprime/recomb_collector.py](ftprime/recomb_collector.py): Provides `RecombCollector`, whose function `collect_recombs` can be used
    as output for simuPOP's `Recombinator` operator.


Tests and examples:

-  [simupop_example.py](simupop_example.py): A simple walk-through of the steps of using with simuPOP.
-  [tests/wf/](test/wf/__init__.py): Very simple forwards-time Wright-Fisher simulation that uses the underlying machinery to the ARGrecorder.
-  [tests/test_arg_recorder_with_wf.py](tests/test_arg_recorder_with_wf.py): Example of using the wf interface.
-  [tests/test_recomb_collector.py](tests/test_recomb_collector.py): Uses a RecombCollector with simuPOP.


Legacy interface:

-  [attic/meiosistagger.py](attic/meiosistagger.py): Provides `MeiosisTagger`, which can be used as an IdTagger in simuPOP
    with the side effect of simulating recombination events and storing everything in an `ARGrecorder`.
-  [tests/dont_test_ftprime_with_simuPOP.py](tests/dont_test_ftprime_with_simuPOP.py): Example of using the legacy simuPOP interface.

[Documentation of the problem and the methods:](writeups/)

Since there are many different ways to store an ARG as a set of coalescence records,
a good deal of this is devoted to describing and verifying msprime's requirements
for such a set, and thinking about different ways to do it.

-  [writeups/forwards_algorithm.md](writeups/forwards_algorithm.md): Description of the algorithm for outputting a valid tree sequence from a forwards-time simulation.
-  [writeups/algorithm-notes.md](writeups/algorithm-notes.md): A writeup of algorithmic considerations.
-  [devel/test_msprime.py](devel/test_msprime.py): An exploration of the capabilities of msprime to digest various notions of a tree sequence.



Development
-----------

Test status and code coverage:

[![CircleCI](https://circleci.com/gh/ashander/ftprime/tree/master.svg?style=svg)](https://circleci.com/gh/ashander/ftprime/tree/master) [![Coverage Status](https://coveralls.io/repos/github/ashander/ftprime/badge.svg?branch=master)](https://coveralls.io/github/ashander/ftprime?branch=master) [![codecov](https://codecov.io/gh/ashander/ftprime/branch/master/graph/badge.svg)](https://codecov.io/gh/ashander/ftprime)


Clone:

    git clone https://github.com/ashander/ftprime.git
    cd ftprime

For best results, use [miniconda](https://conda.io/miniconda.html),
which provides the command line dependency manager `conda`.
Once you have it installed, make a new environment to do development:

    conda config --add channels conda-forge
    conda env create -f environment.yml -n ftprime python=3.5
    source activate ftprime  # Enter the development environment

Install ``tortoisim`` in locally editable (``-e``) mode and run the tests.
After the ``pip`` command you should see a bunch of messages about requirements
already satisfied (because you've installed them with ``conda``, above):

    pip install -e .[dev]  # Don't need the [dev] if you used conda above
    pytest

To-do
=====

1. Change recombination to be Poisson, and **check** we are not missing a factor of 2 in the definition of Morgans.
