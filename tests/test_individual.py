import pytest
from ftprime import Individual, Chromosome


def test_making_and_using():
    with pytest.raises(Exception) as e_info:
        Individual((1))
    with pytest.raises(Exception) as e_info:
        Individual((None, None))

    c = Chromosome(1, .5, 0, 2)
    with pytest.raises(Exception) as e_info:
        Individual((c, ))
    with pytest.raises(Exception) as e_info:
        Individual((c, 1))

    i = Individual((c, c))
    i.age()
    print([ch for ch in i])
    print([g for g in i.raw_gametes(5)])
    with pytest.raises(Exception) as e_info:
        print([g for g in i.raw_gametes(4.5)])


def test_raw_gametes_and_gametes():
    p1_id = 1
    p2_id = 3
    c1 = Chromosome(p1_id, .5, 0, 2)
    c2 = Chromosome(p2_id, .5, 0, 2)
    i = Individual((c1, c2))
    i_rgams = [g for g in i.raw_gametes(5)]
    assert len(i_rgams) == 5
    for g in iter(i_rgams):
        assert len(g) == 3
        # the chrome has both parents (no full-parents in this sim)
        assert p2_id in g
        assert p1_id in g
        ch_args = (5, ) + g
        Chromosome(*ch_args)

    ids = (101, 151, 2015, 201, 1510, 1029)
    i_gams = [c for c in i.gametes(5, tagger=iter(ids))]
    assert len(i_gams) == 5
    for c, idnum in zip(iter(i_gams), iter(ids)):
        pars = (c.lparent, c.rparent)
        assert type(c) is Chromosome
        assert p2_id in pars
        assert p1_id in pars
        assert c.breakpoint <= 1.0 and c.breakpoint >= 0.0
        assert c.identity == idnum
