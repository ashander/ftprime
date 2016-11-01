import pytest
from ftprime import Population
from ftprime.chromosome import default_chromosome
from numpy.random import randint
import msprime


def test_creation():
    with pytest.raises(Warning):
        Population(size=1.2, verbose=True)


def test_creation_int():
    pops = randint(low=100, high=10000, size=10)
    for n in pops:
        p = Population(size=2)
        realsize = len([i for i in p])
        numchroms = len([c for c in p.chromosomes()])
        assert realsize == p.size()
        assert numchroms == p.size() * 2


def test_initialize():
    __import__("numpy").random.seed(1221)
    p = Population(size=2)
    for k, seglist in p.unmerged_records():
        for seg in seglist:
            assert seg.node is default_chromosome['rparent']
            assert seg.right is default_chromosome['right_end']
            assert seg.left is default_chromosome['left_end']


def test_working_generation():
    __import__("numpy").random.seed(1221)
    pops = randint(low=100, high=1000, size=10)
    for sz in pops:
        p = Population(size=sz)
        for i in p:
            print(i)
            pass
        for c in p.chromosomes():
            print(c)
            pass
        p.generation(2, True)
        print(p._nc)
        p.finalize()
        with open('working.tsv', 'w') as f:
            p.write_records(f)
        msprime.load_txt('working.tsv')
