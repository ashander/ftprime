import ftprime
import msprime
import numpy as np
import unittest

from tests import *

from .wf import wf_call


class WfTestCase(FtprimeTestCase):

    def run_wf(self, N, ngens, nsamples, survival=0.0, simplify_interval=10,
               mutation_rate=0.0):
        records = wf_call(N=N, ngens=ngens, nsamples=nsamples, survival=survival,
                          debug=False, simplify_interval=simplify_interval,
                          seed=self.random_seed, mutation_rate=mutation_rate)
        return records

    def check_tables(self, records):
        nodes = records.tables.nodes
        edges = records.tables.edges
        # check edges are in order and all parents are recorded
        node_times = nodes.time
        last_time = 0.0
        for p in edges.parent:
            self.assertTrue(node_times[p] >= last_time)
            last_time = node_times[p]
            self.assertTrue(p < records.tables.nodes.num_rows)
        for ch in edges.child:
            self.assertTrue(ch < records.tables.nodes.num_rows)

    def test_runs(self):
        N = 10
        ngens = 20
        records = self.run_wf(N=N, ngens=ngens, nsamples=N)
        self.check_tables(records)

    def test_simplify_interval(self):
        # since all randomness is in wf, should get *exactly the same trees*
        # running with different simplify_intervals.
        N = 5
        ngens = 20
        for mut_rate in [0.0, 0.1]:
            print("-------------------------\n")
            print("Mut rate:" + str(mut_rate) + "\n")
            records_a = self.run_wf(N=N, ngens=ngens, nsamples=N, simplify_interval=ngens,
                                    mutation_rate=mut_rate)
            records_b = self.run_wf(N=N, ngens=ngens, nsamples=N, simplify_interval=2,
                                    mutation_rate=mut_rate)
            records_c = self.run_wf(N=N, ngens=ngens, nsamples=N, simplify_interval=100,
                                    mutation_rate=mut_rate)
            self.assertEqual(records_a.num_simplifies, 1+1)
            self.assertEqual(records_b.num_simplifies, 10+1)
            self.assertEqual(records_c.num_simplifies, 0+1)
            sample_ids = [N*ngens + x for x in range(N)]
            self.check_trees(records_a.tree_sequence(sample_ids),
                             records_b.tree_sequence(sample_ids))
            self.check_trees(records_a.tree_sequence(sample_ids),
                             records_c.tree_sequence(sample_ids))

    def test_get_nodes(self):
        N = 10
        ngens = 20
        records = self.run_wf(N=N, ngens=ngens, nsamples=N)
        # this should be the input IDs for final gen
        final_gen = np.random.choice([N*ngens + x for x in range(N)], N,
                                     replace=False)
        records.check_ids(final_gen)
        final_nodes = records.get_nodes(final_gen)
        flags = records.tables.nodes.flags
        for k in range(N):
            self.assertTrue(final_nodes[k] < records.tables.nodes.num_rows)
            self.assertEqual(records.node_ids[final_gen[k]], final_nodes[k])
            self.assertEqual(flags[final_nodes[k]], msprime.NODE_IS_SAMPLE)


    def test_sample_ids(self):
        records = self.run_wf(N=11, ngens=20, nsamples=5)
        sample_ids = records.sample_ids()
        sample_nodes = [records.node_ids[k] for k in sample_ids]
        flags = records.tables.nodes.flags
        print(records)
        print("samples:", sample_ids)
        print("sample nodes:", sample_nodes)
        for k in range(len(flags)):
            if k in sample_nodes:
                self.assertEqual(flags[k], msprime.NODE_IS_SAMPLE)
            # # currently, simplify doesn't reset sample flag
            # else:
            #     self.assertEqual(flags[k], 0)

    def test_overlapping_generations(self):
        records = self.run_wf(N=11, ngens=20, nsamples=5, survival=0.5)

        self.check_tables(records)

        print(records)

        samples = records.sample_ids()
        print("samples:", samples)
        ts = records.tree_sequence()

        for t in ts.trees():
            print(t)
