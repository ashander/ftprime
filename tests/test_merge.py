from ftprime.merge import merge_records, CoalescenceRecord

def test_singleton_non_overlapping():
    c1 = CoalescenceRecord(0.0, 0.3, 3, (1, ), 0.0, 0)
    c2 = CoalescenceRecord(0.4, 1.0, 3, (2, ), 0.0, 0)
    comp, inc = merge_records([c1, c2], debug=True)
    assert len(comp) == 0
    assert c1 in inc
    assert c2 in inc


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
