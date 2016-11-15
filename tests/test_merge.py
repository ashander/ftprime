from ftprime.merge import merge_records, CoalescenceRecord
from ftprime import merge

c1 = CoalescenceRecord(0.0, 0.5, 3, (1, 2), 0.0, 0)
c2 = CoalescenceRecord(0.3, 1.0, 3, (1, 2), 1.0, 0)
c5 = CoalescenceRecord(0.4, 1.0, 3, (1, 2, 7), 0.0, 0)
c6 = CoalescenceRecord(0.0, 0.2, 3, (1, 2, 7), 1.0, 0)
c3 = CoalescenceRecord(0.0, 0.3, 3, (1, ), 0.0, 0)
c4 = CoalescenceRecord(0.4, 1.0, 3, (2, ), 0.0, 0)
crecs = (c1, c2, c5, c6)
irecs = (c3, c4)


def test_ll_flatten():
    childlist = [(1, 2), (1, )]
    assert merge._flatten(childlist) == (1, 2)
    childlist = [(1, 2), (1, ), (7, 8, 9, 20)]
    assert merge._flatten(childlist) == (1, 2, 7, 8, 9, 20)
    childlist = [(1, ), ]
    assert merge._flatten(childlist) == (1, )


def test_ll_prep_records_overlapping():
    ep, ch, t, add, ind =  merge._prepare_records_to_merge([c1, c2])
    assert ep == (0.0, 0.3, 0.5, 1.0)
    assert len(ch) == 4
    for c in ch:
        assert c == (1,2)
    assert t == (0.0, 1.0, 0.0, 1.0)
    assert add == (True, True, False, False)
    assert ind == (0, 1, 0, 1)


def test_ll_prep_records_nonoverlapping_sibs():
    ep, ch, t, add, ind =  merge._prepare_records_to_merge([c5, c6])
    assert ep == (0.0, 0.2, 0.4, 1.0)
    assert len(ch) == 4
    for c in ch:
        assert c == (1,2, 7)
    assert t == (1.0, 1.0, 0.0, 0.0)
    assert add == (True, False, True, False)
    assert ind == (1, 1, 0, 0)


def test_ll_prep_records_nonoverlapping_orphans():
    ep, ch, t, add, ind =  merge._prepare_records_to_merge([c3, c4])
    assert ep == (0.0, 0.3, 0.4, 1.0)
    assert len(ch) == 4
    assert (1, ), (1, ) == ch[:2]
    assert (2, ), (2, ) == ch[2:]
    assert t == (0.0, 0.0, 0.0, 0.0)
    assert add == (True, False, True, False)
    assert ind == (0, 0, 1, 1)


def test_ll_prep_records_overlapping_children_vary():
    ep, ch, t, add, ind =  merge._prepare_records_to_merge([c1, c5])
    assert ep == (0.0, 0.4, 0.5, 1.0)
    assert len(ch) == 4
    assert (1, 2), (1, 2, 7) == ch[:2]
    assert (1, 2), (1, 2, 7) == ch[2:]
    assert t == (0.0, 0.0, 0.0, 0.0)
    assert add == (True, True, False, False)
    assert ind == (0, 1, 0, 1)

    ep, ch, t, add, ind =  merge._prepare_records_to_merge([c3, c6])
    assert ep == (0.0, 0.0, 0.2, 0.3)
    assert len(ch) == 4
    assert (1, ), (1, 2, 7) == ch[:2]
    assert (1, ), (1, 2, 7) == ch[2:]
    assert t == (0.0, 1.0, 1.0, 0.0)
    assert add == (True, True, False, False)
    assert ind == (0, 1, 1, 0)


def test_ll_in_or_complete_records():
    comp = [c for c in merge._complete_records(crecs + irecs)]
    inc = [i for i in merge._incomplete_records(crecs + irecs)]
    assert len(inc) == len(irecs)
    assert len(comp) == len(crecs)
    for c in crecs:
        assert c in comp
    for i in irecs:
        assert i in inc
    # just a couple
    comp = [c for c in merge._complete_records(crecs)]
    inc = [i for i in merge._incomplete_records(crecs)]
    assert len(inc) == 0
    assert len(comp) == len(crecs)
    for c in crecs:
        assert c in comp
    # just a couple
    comp = [c for c in merge._complete_records(irecs)]
    inc = [i for i in merge._incomplete_records(irecs)]
    assert len(inc) == len(irecs)
    assert len(comp) == 0
    for i in irecs:
        assert i in inc


def test_singleton_non_overlapping():
    comp, inc = merge_records([c3, c4], debug=True)
    assert len(comp) == 0
    assert len(inc) == 2
    assert c3 in inc
    assert c4 in inc


def test_sibs_non_overlapping():
    c1 = CoalescenceRecord(0.0, 0.3, 3, (1, 2), 0.0, 0)
    c2 = CoalescenceRecord(0.4, 1.0, 3, (1, 2), 0.0, 0)
    comp, inc = merge_records([c1, c2], debug=True)
    assert len(comp) == 2
    assert len(inc) == 0  # completeness is >= 2 children
    assert c1 in comp
    assert c2 in comp


def test_simple_overlapping():
    c1 = CoalescenceRecord(0.0, 0.5, 3, (1, 2), 0.0, 0)
    c2 = CoalescenceRecord(0.3, 1.0, 3, (1, 2), 0.0, 0)
    comp, inc = merge_records([c1, c2], debug=True)
    assert len(comp) == 1
    assert len(inc) == 2
    assert c1 in inc
    assert c2 in inc
