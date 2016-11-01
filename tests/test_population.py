import pytest
from ftprime import Population
from ftprime.chromosome import default_chromosome
from numpy.random import randint, choice, seed
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
    seed(1221)
    p = Population(size=2)
    for k, seglist in p.unmerged_records():
        for seg in seglist:
            assert seg.node is default_chromosome['rparent']
            assert seg.right is default_chromosome['right_end']
            assert seg.left is default_chromosome['left_end']


def test_working_generation():
    seed(1221)
    pop_size = randint(low=100, high=500, size=20)
    for i, sz in enumerate(pop_size):
        p = Population(size=sz)
        p.generation(2, True)
        p.finalize()
        with open('working.tsv', 'w') as f:
            p.write_records(f)

        tr = msprime.load_txt('working.tsv')
        all_samples = tr.get_samples()
        print('total samples:', len(all_samples))
        # subset requires child of commit 21be37b in msprime
        with pytest.raises(Exception) as e_info:
            ts = tr.subset(choice(all_samples, size=10, replace=False))
        print(i)
        if i in (0, 2, 3, 5, 12, 15, 19):
            assert 'LibraryError' in str(e_info)
        elif i in (1, 4, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18):
            assert 'IndexError' in str(e_info)
