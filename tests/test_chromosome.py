import pytest
from ftprime.chromosome import (
    Chromosome,
    default_chromosome,
    Segment,
    default_segment
)


def test_fails_wo_required():
    with pytest.raises(Exception) as e:
        Chromosome(identity=1)
        assert type(e) is TypeError

    with pytest.raises(Exception) as e:
        Chromosome(breakpoint=.5)
        assert type(e) is TypeError


def test_defaults():
    c = Chromosome(1)
    print(c)
    assert c.identity == 1
    assert c.breakpoint == default_chromosome['breakpoint']
    assert c.lparent == default_chromosome['lparent']
    assert c.rparent == default_chromosome['rparent']
    assert c.left_end == default_chromosome['left_end']
    assert c.right_end == default_chromosome['right_end']

    s = Segment()
    print(s)
    ds = default_segment
    assert s.left == ds['left']
    assert s.right == ds['right']
    assert s.time == ds['time']
    assert s.population == ds['population']
    assert s.node == ds['node']


def test_position_args():
    parg = (101, 1.9, 1, 5, 5, 2)
    c = Chromosome(*parg)
    print(c)
    for i, j in zip(c, parg):
        assert i == j
    parg2 = (0.0, 0.9, 1, (5,), 1.0, 1)
    s = Segment(*parg2)
    print(s)


def test_varying_args():
    varg = (101, 1.9, 1, 5, 5, 2)
    c = Chromosome(
        breakpoint=varg[1],
        lparent=varg[2],
        rparent=varg[3],
        identity=varg[0],
    )
    print(c)

    assert c.breakpoint == varg[1]
    assert c.lparent == varg[2]
    assert c.rparent == varg[3]
    assert c.identity == varg[0]
    assert c.left_end == default_chromosome['left_end']
    assert c.right_end == default_chromosome['right_end']
