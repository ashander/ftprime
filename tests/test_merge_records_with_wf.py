from ftprime import merge_records
import msprime
import random
from wf import wf

def check_tables(args):
    nodes = args.node_table()
    assert(nodes.num_rows == args.num_nodes)
    edgesets = args.edgeset_table()
    # check edgesets are in order and all parents are recorded
    node_times = nodes.time
    last_time = 0.0
    for p in edgesets.parent:
        assert(node_times[p] >= last_time)
        last_time = node_times[p]
        assert(p < args.num_nodes)
    for ch in edgesets.children:
        assert(ch < args.num_nodes)

def test_1():

    new_rec=msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(22,))
    old_recs=[
            msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19,)),
            msprime.Edgeset(left=0.6, right=1.0, parent=18, children=(19, 20))]
    merge_records(new_rec,old_recs)

    old_recs

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19,)),
        msprime.Edgeset(left=0.6, right=0.8, parent=18, children=(19,20)),
        msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(19, 20, 22)),
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')

def test_2():

    new_rec = msprime.Edgeset(left=0.5, right=1.0, parent=4, children=(13,))
    old_recs = [ msprime.Edgeset(left=0.0, right=0.9, parent=4, children=(7, 12)),
                 msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7,)) ]
    merge_records(new_rec,old_recs)

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.5, parent=4, children=(7, 12)),
        msprime.Edgeset(left=0.5, right=0.9, parent=4, children=(7, 12, 13)),
        msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7, 13))
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')

def test_3():

    new_rec = msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,))
    old_recs = [msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,))]
    merge_records(new_rec,old_recs)

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,)),
        msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,))
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')


def test_simulation_runs():

    random.seed(123)
    records = wf(N=5,ngens=5,nsamples=5,survival=0.5)

    check_tables(records)

    for x in records:
        print(x, records[x])

    print(records.edgeset_table())
    print(records.node_table())

    ts = records.tree_sequence()

    for t in ts.trees():
        print(t)

    print("Mean pairwise diversity:",ts.get_pairwise_diversity())
    print("(should be zero)")
    assert ts.get_pairwise_diversity() == 0.0
