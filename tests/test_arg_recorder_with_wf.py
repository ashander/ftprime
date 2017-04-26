import pytest
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

@pytest.mark.parametrize(('N', 'gen', 'samples'), [
    (2, 2, 2),
    (5, 5, 5),
    (10, 10, 5),
])
def test_simulation_runs(N, gen, samples):

    random.seed(123)
    records = wf(N=N, ngens=gen, nsamples=samples, survival=0.5)

    # check_tables(records)

    # for x in records:
    #     print(x, records[x])

    # print(records.edgeset_table())
    # print(records.node_table())

    # ts = records.tree_sequence()

    # for t in ts.trees():
    #     print(t)

    # print("Mean pairwise diversity:",ts.get_pairwise_diversity())
    # print("(should be zero)")
    # assert ts.get_pairwise_diversity() == 0.0

if __name__ == '__main__':
    test_simulation_runs(1000, 200, 100)
