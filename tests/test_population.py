import pytest
from ftprime import Population
from numpy.random import randint


def test_creation():
    with pytest.raises(Warning):
        Population(size=1.2, verbose=True)


def test_creation_int():
    pops = randint(low=100, high=10000, size=10)
    for n in pops:
        p = Population(size=2)
        realsize = len([i for i in p.individuals()])
        numchroms = len([c for c in p.chromosomes()])
        assert realsize == p.size()
        assert numchroms == p.size() * 2


def test_generation():
    p = Population(size=5)
    for i in p.individuals():
        print(i)
    for c in p.chromosomes():
        print(c)
    assert False

    for o in p.generation(15):
        print(o)

