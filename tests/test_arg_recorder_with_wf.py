import pytest
import random

from wf import wf


def check_tables(args):
    nodes = args.nodes
    assert(nodes.num_rows == args.num_nodes)
    edgesets = args.edgesets
    # check edgesets are in order and all parents are recorded
    node_times = nodes.time
    last_time = 0.0
    for p in edgesets.parent:
        assert(node_times[p] >= last_time)
        last_time = node_times[p]
        assert(p < args.num_nodes)
    for ch in edgesets.children:
        assert(ch < args.num_nodes)

@pytest.mark.parametrize(('N', 'gen', 'nsamples'), [
    (2, 2, 2),
    (5, 5, 5),
    (10, 10, 5),
])
def test_simulation_runs(N, gen, nsamples):

    random.seed(123)
    records = wf(N=N, ngens=gen, nsamples=nsamples, survival=0.5, debug=True)

    check_tables(records)

    for x in records:
        print(x, records[x])

    print(records.edgesets())
    print(records.nodes())

    ts = records.tree_sequence()

    for t in ts.trees():
        print(t)

    print("Mean pairwise diversity:",ts.get_pairwise_diversity())
    print("(should be zero)")
    assert ts.get_pairwise_diversity() == 0.0
