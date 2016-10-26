import pytest
from ftprime import Chromosome, default_chromosome_values


def test_fails_wo_required():
    with pytest.raises(Exception) as e:
        Chromosome(identity=1)
        assert type(e) is TypeError

    with pytest.raises(Exception) as e:
        Chromosome(breakpoint=.5)
        assert type(e) is TypeError


def test_defaults():
    c = Chromosome(1, 0.5)
    print(c)
    assert c.identity == 1
    assert c.breakpoint == 0.5
    assert c.lparent == default_chromosome_values[0]
    assert c.rparent == default_chromosome_values[1]
    assert c.left_end == default_chromosome_values[2]
    assert c.right_end == default_chromosome_values[3]



def test_position_args():
    parg = (101, 1.9, 1, 5, 5, 2)
    c = Chromosome(*parg)
    print(c)
    for i, j in zip(c, parg):
        assert i == j


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
    assert c.left_end == default_chromosome_values[2]
    assert c.right_end == default_chromosome_values[3]
