import ftprime
import msprime
import random
import unittest

from tests import *

from .wf import wf


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

class WfTestCase(FtprimeTestCase):

    def run_wf(self, N, ngens, nsamples, survival=0.0):
        records = wf(N=N, ngens=ngens, nsamples=nsamples, survival=survival,
                     debug=True, seed=self.random_seed)
        return records

    def test_get_nodes(self):
        N = 2
        ngens = 2
        records = self.run_wf(N=N, ngens=ngens, nsamples=N)
        # this should be the input IDs for final gen
        final_gen = random.sample([N*ngens + x for x in range(N)], N)
        records.check_ids(final_gen)
        final_nodes = records.get_nodes(final_gen)
        flags = records.nodes.flags
        for k in range(N):
            self.assertTrue(final_nodes[k] < records.nodes.num_rows)
            self.assertEqual(records.node_ids[final_gen[k]], final_nodes[k])
            self.assertEqual(flags[final_nodes[k]], msprime.NODE_IS_SAMPLE)


    def test_sample_ids(self):
        records = self.run_wf(N=2, ngens=2, nsamples=2)
        sample_ids = records.sample_ids()
        sample_nodes = [records.node_ids[k] for k in sample_ids]
        flags = records.nodes.flags
        print(records)
        print("samples:", sample_ids)
        print("sample nodes:", sample_nodes)
        for k in range(len(flags)):
            if k in sample_nodes:
                self.assertEqual(flags[k], msprime.NODE_IS_SAMPLE)
            # # currently, simplify doesn't reset sample flag
            # else:
            #     self.assertEqual(flags[k], 0)

    @unittest.skip
    def test_nonoverlapping_generations(self):

        N = 3
        ngens = 2
        nsamples = 3
        records = wf(N=N, ngens=ngens, nsamples=nsamples, survival=0.0, 
                     debug=True, seed=self.random_seed)

        check_tables(records)

        print(records)

        samples = records.sample_ids()
        print("samples:", samples)
        ts = records.tree_sequence()

        for t in ts.trees():
            print(t)
