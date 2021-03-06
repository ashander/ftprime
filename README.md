ftprime
======

An earlier version of this code was mentioned during the Evolution 2017 talk:
> Ashander, McCartney-Melstad, Ralph, Shaffer (2017) "Using Genomic Data to Inform Population Viability in a Long-Lived Endangered Vertebrate". 

Note the code is pre-production but if you're building on the ideas here please cite this repository using the DOI: [![DOI](https://zenodo.org/badge/72480698.svg)](https://zenodo.org/badge/latestdoi/72480698)

The current version of this package (on `HEAD`) works with msprime version 0.6.1 (current as of Oct 17 2018).
For the earlier version, used in [the paper](https://www.biorxiv.org/content/early/2018/06/07/248500), get [tag 0.0.6](https://github.com/ashander/ftprime/tree/0.0.6).


Contents
--------

The purpose of this package is to provide python code to easily store ancestry information from a forwards-time, 
individual-based simulation, using [msprime](https://github.com/jeromekelleher/msprime)'s "tree sequence" data format, so that after the simulation,
we can

1. have the entire tree sequence (equivalent to the ARG) of the final generation:w
2. put down neutral mutations on the tree sequence afterwards without carrying them along in the simulation, and
3. use [msprime](https://github.com/jeromekelleher/msprime) to efficiently store results and quickly compute statistics.

We also provide a class to facilitate doing this with [simuPOP](https://github.com/BoPeng/simuPOP).

[Core module:](ftprime/)

-  [ftprime/argrecorder.py](ftprime/argrecorder.py): Provides `ARGrecorder`, which records haploid parentage and crossing over events.

-  [ftprime/recomb_collector.py](ftprime/recomb_collector.py): Provides `RecombCollector`, which does the bookkeeping to use ARGrecorder
    with a diploid simulation with discrete loci, whose function `collect_recombs` can be used as output for simuPOP's `Recombinator` operator.


Tests:

-  [tests/wf/](test/wf/__init__.py): Very simple forwards-time Wright-Fisher simulation that uses the underlying machinery to the `ARGrecorder`.
-  [tests/test_with_wf.py](tests/test_with_wf.py): Example of using the Wright-Fisher interface.
-  [examples/](examples/): more complex examples using `RecombCollector`.


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

Earlier work
------------

Legacy interface:

-  [attic/meiosistagger.py](attic/meiosistagger.py): Provides `MeiosisTagger`, which can be used as an IdTagger in simuPOP
    with the side effect of simulating recombination events and storing everything in an `ARGrecorder`.

[Documentation of the problem and the methods:](writeups/)

Since there are many different ways to store an ARG as a set of coalescence records,
a good deal of this is devoted to describing and verifying msprime's requirements
for such a set, and thinking about different ways to do it.

-  [writeups/forwards_algorithm.md](writeups/forwards_algorithm.md): Description of the algorithm for outputting a valid tree sequence from a forwards-time simulation.
-  [writeups/algorithm-notes.md](writeups/algorithm-notes.md): A writeup of algorithmic considerations.
-  [devel/test_msprime.py](devel/test_msprime.py): An exploration of the capabilities of msprime to digest various notions of a tree sequence.
