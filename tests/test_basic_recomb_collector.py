import ftprime
import msprime
import six
import random
import math

from tests import FtprimeTestCase


class RecombCollectorTest(FtprimeTestCase):

    def simple_ex(self):
        # this will begin with a single diploid indiv
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       0           0               1.00000000000000
        1       1           1               0.00000000000000
        2       1           2               0.00000000000000
        """)
        edges = six.StringIO("""\
        id      left            right           parent  child
        0       0.00000000      3.00000000      0       1
        1       0.00000000      3.00000000      0       2
        """)
        init_ts = msprime.load_text(nodes=nodes, edges=edges)
        # diploid 0 maps initially to haploids 1 and 2 in init_ts:
        node_ids = {(0,0):1, (0,1):2}
        locus_position = [0.0, 1.0, 2.0, 3.0]
        rc = ftprime.RecombCollector(ts=init_ts, node_ids=node_ids, 
                                     locus_position=locus_position,
                                     benchmark=True)
        assert rc.mode == 'text'
        rc2 = ftprime.RecombCollector(ts=init_ts, node_ids=node_ids, 
                                     locus_position=locus_position,
                                     benchmark=True, mode='binary')
        assert rc2.mode == 'binary'
        return rc, node_ids

    def medium_ex(self):
        """
        below we add two new offspring with no recombinations
        this adds 2 nodes and 3 edges
        """
        rc, node_ids = self.simple_ex()
        self.assertListEqual(rc.locus_position,
                             [0.0, 1.0, 2.0, 3.0])
        # Input is pairs of
        #     offspringID parentID startingPloidy rec1 rec2 ....
        # in pairs of offpsring chromosomes
        lines_list = ["""
        1   0   1
        1   0   0   1
        """]
        for lines in lines_list:
            rc.increment_time()
            rc.collect_recombs(lines)
        rc.args.update_times()
        expected_nodes = 3 + 2
        expected_edges = 2 + 3
        return rc, expected_edges, expected_nodes

    def bigger_ex(self):
        """
        Below we have this situation ('indiv' is the diploid ID)
            indiv   chromosome  parent:[left,right) ...
            # at t=1.0
            1       0           (0,1):[0,3)
            1       1           (0,0):[0,3)
            2       0           (0,1):[0,0.5) (0,0):[0.5,3)
            2       1           (0,0):[0,1.5) (0,1):[1.5,3)
            3       0           (0,0):[0,0.5) (0,1):[0.5,3)
            3       1           (0,1):[0,1.5) (0,0):[1.5,3)
            # at t=2.0
            4       0           (2,0):[0,0.5) (2,1):[0.5,1.5)
            4       1           (1,1):[0,0.5) (1,0):[0.5,3)
            5       0           (1,1):[0,0.5) (1,0):[0.5,3)
            5       1           (2,0):[0,0.5) (2,1):[0.5,1.5) (2,0):[1.5,2.5) (2,1):[2.5,3)
        and since since
            locus_position = [0.0, 1.0, 2.0, 3.0]
        we'll say recombination "at 1" above is somewhere before 1.0, i.e.
        after locus 0.

        This has trees (with ** the ancestor and omitting indiv 3)
                   **                      **                        **                       **
                  /  \                    /  \                      /  \                     /  \
              (0,0) (0,1)              (0,0) (0,1)              (0,0) (0,1)              (0,0) (0,1)
             /   |    |  \            /   | \     \            /   |   |   \            /   |   |   \
        (1,1) (2,1) (2,0) (1,0)  (1,1) (2,1) (2,0) (1,0)  (1,1) (2,0) (2,1) (1,0)  (1,1) (2,0) (2,1) (1,0)
           | \        |  \            /  |        / |          /   |       / |          /     /     / |
        (4,1) (5,0) (4,0) (5,1)  (4,0) (5,1) (4,1) (5,0)  (4,0) (5,1) (4,1) (5,0)  (4,0) (5,1) (4,1) (5,0)
                               0.5                     1.5

        ... and above 2.5 is the same as 0.5 to 1.5

        After simplifying to [4,5] should be:
                   9.                      9.                        9.                       9.
                  /  \                    /  \                      /  \                     /  \
              ---/    \                   |   \___                 |    \__                 /  --8--
             /        |                   |       \               /        \               /    /   \
        --5--       --6--              --7--       --4--       --6--       --4--         /     /   --4--
           | \        |  \            /  |        / |          /  |        / |          /    /      / |
        --1-- --2-- --0-- --3--  --0-- --3-- --1-- --2--  --0-- --3-- --1-- --2--  --0-- --3-- --1-- --2--
                               0.5                     1.5
        """
        rc, node_ids = self.simple_ex()
        self.assertListEqual(rc.locus_position, 
                             [0.0, 1.0, 2.0, 3.0])
        # Input is pairs of
        #     offspringID parentID startingPloidy rec1 rec2 ....
        # in pairs of offpsring chromosomes
        lines_list = ["""
        1   0   1
        1   0   0
        2   0   1   0
        2   0   0   1
        3   0   0   0
        3   0   1   1
        """, """
        4   2   0   0 1
        4   1   1   0
        5   1   1   0
        5   2   0   0 1 2
        """]
        for lines in lines_list:
            rc.increment_time()
            rc.collect_recombs(lines)
        rc.args.update_times()
        return rc

    def check_node_ids(self, rc, haploid_node_ids):
        rc, node_ids = self.simple_ex()
        for h in haploid_node_ids:
            u = rc.i2c(h[0], h[1])
            self.assertTrue(u in rc.args.node_ids)
            self.assertEqual(rc.args.node_ids[u], haploid_node_ids[h])

    def test_init(self):
        rc, node_ids = self.simple_ex()
        self.assertEqual(rc.sequence_length, 3.0)
        self.assertEqual(rc.args.nodes.num_rows, 3)
        self.assertEqual(rc.time, 0.0)
        self.check_node_ids(rc, node_ids)

    def test_simple_i2c(self):
        rc, node_ids = self.simple_ex()
        self.assertEqual(rc.i2c(0, 0), 0)
        self.assertEqual(rc.i2c(0, 1), 1)
        self.assertEqual(rc.i2c(1, 0), 2)
        self.assertEqual(rc.i2c(1, 1), 3)
        self.assertEqual(rc.i2c(2, 0), 4)
        self.assertEqual(rc.i2c(2, 1), 5)
        with self.assertRaises(ValueError):
            rc.i2c(2, -1)
        with self.assertRaises(ValueError):
            rc.i2c(2, 2)

    def test_simple_i2n(self):
        rc, node_ids = self.simple_ex()
        self.assertEqual(rc.i2n(0, 0), 1)
        self.assertEqual(rc.i2n(0, 1), 2)
        with self.assertRaises(ValueError):
            rc.i2n(2, -1)
        with self.assertRaises(ValueError):
            rc.i2n(2, 2)

    def test_increment_time(self):
        rc, node_ids = self.simple_ex()
        self.assertEqual(rc.time, 0.0)
        rc.increment_time()
        self.assertEqual(rc.time, 1.0)

    def test_collect_recombs_simple(self):
        rc, ee, en = self.medium_ex()
        self.assertEqual(rc.time, 1.0)
        print(rc.args.node_ids)
        print(rc.args)
        nodes = rc.args.nodes
        print(nodes)
        self.assertEqual(nodes.num_rows, en)
        edges = rc.args.edges
        self.assertArrayEqual(nodes.time,
                              [2.0, 1.0, 1.0] + [0.0, 0.0])
        print(edges)
        self.assertEqual(edges.num_rows, ee)
        # Node IDs:
        # {0: 1, 1: 2, 2: 3, 3: 4}
        # Nodes:
        # id	flags	population	time
        # 0	    0	    0		    2.00000000000000
        # 1	    1	    1		    1.00000000000000
        # 2	    1	    2		    1.00000000000000
        # 3	    1	    -1		    0.00000000000000
        # 4	    1	    -1		    0.00000000000000
        # Edges:
        # id	left		right		parent	child
        # 0	0.00000000	3.00000000	0	1
        # 1	0.00000000	3.00000000	0	2
        # 2	0.00000000	3.00000000	2	3
        # 3	0.00000000	3.00000000	1	4

    def test_collect_recombs(self):
        rc = self.bigger_ex()
        self.assertEqual(rc.time, 2.0)
        # check nodes
        nodes = rc.args.nodes
        print(nodes)
        self.assertEqual(nodes.num_rows, 3+10)
        self.assertArrayEqual(nodes.time,
                [3.0, 2.0, 2.0]
                 + [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        # check edges
        edges = rc.args.edges
        print(edges)
        self.assertEqual(edges.num_rows, 2 + 21)
        dip_parents = [[(0,1)],
                       [(0,0)],
                       [(0,1), (0,0)],
                       [(0,0), (0,1)],
                       [(0,0), (0,1)],
                       [(0,1), (0,0)],
                       [(2,0), (2,1), (2,0)],
                       [(1,1), (1,0)],
                       [(1,1), (1,0)],
                       [(2,0), (2,1), (2,0), (2,1)]]
        true_parents = [0, 0] + [rc.i2n(*y) for x in dip_parents for y in x]
        dip_children = [[(1,0)],
                        [(1,1)],
                        [(2,0), (2,0)],
                        [(2,1), (2,1)],
                        [(3,0), (3,0)],
                        [(3,1), (3,1)],
                        [(4,0), (4,0), (4,0)],
                        [(4,1), (4,1)],
                        [(5,0), (5,0)],
                        [(5,1), (5,1), (5,1), (5,1)]]
        true_children = [1, 2] + [rc.i2n(*y) for x in dip_children for y in x]
        self.assertArrayEqual(edges.parent, true_parents)
        self.assertArrayEqual(edges.child, true_children)
        obs_lefts = [math.floor(z) for z in edges.left]
        obs_rights = [math.ceil(z) for z in edges.right]
        true_lefts = [0, 0] + ([0] +
                            [0] +
                            [0, 0] +
                            [0, 1] +
                            [0, 0] +
                            [0, 1] +
                            [0, 0, 1] +
                            [0, 0] +
                            [0, 0] +
                            [0, 0, 1, 2])
        true_rights = [3, 3] + ([3] +
                             [3] +
                             [1, 3] +
                             [2, 3] +
                             [1, 3] +
                             [2, 3] +
                             [1, 2, 3] +
                             [1, 3] +
                             [1, 3] +
                             [1, 2, 3, 3])
        print(true_rights)
        print(obs_rights)
        self.assertArrayEqual(true_lefts, obs_lefts)
        self.assertArrayEqual(true_rights, obs_rights)

    def test_i2c(self):
        rc = self.bigger_ex()
        self.assertEqual(rc.i2c(0, 0), 0)
        self.assertEqual(rc.i2c(0, 1), 1)
        self.assertEqual(rc.i2c(1, 0), 2)
        self.assertEqual(rc.i2c(1, 1), 3)
        self.assertEqual(rc.i2c(2, 0), 4)
        self.assertEqual(rc.i2c(2, 1), 5)

    def test_i2n(self):
        rc = self.bigger_ex()
        self.assertEqual(rc.i2n(0, 0), 1)
        self.assertEqual(rc.i2n(0, 1), 2)
        self.assertEqual(rc.i2n(1, 0), 3)
        self.assertEqual(rc.i2n(1, 1), 4)
        self.assertEqual(rc.i2n(2, 0), 5)
        self.assertEqual(rc.i2n(2, 1), 6)

    def test_add_locations(self):
        rc = self.bigger_ex()
        dip_indivs = range(1,6)
        dip_locations = [2+x for x in random.sample(dip_indivs, len(dip_indivs))]
        rc.add_locations(input_ids=dip_indivs, locations=dip_locations)
        true_locations = [dip_locations[k-1] for k in dip_indivs for p in [0,1]]
        obs_locations = [rc.args.nodes.population[rc.i2n(k,p)] for k in dip_indivs for p in [0,1]]
        print(true_locations)
        print(obs_locations)
        self.assertArrayEqual(true_locations, obs_locations)

    def test_simple_simplify(self):
        rc, node_ids = self.simple_ex()
        rc.simplify([0])
        self.assertArrayEqual(rc.args.nodes.time, [0,0,1])
        self.assertArrayEqual(rc.args.edges.parent, [2, 2])
        self.assertArrayEqual(rc.args.edges.child, [0,1])

    def test_simplify(self):
        rc = self.bigger_ex()
        dip_indivs = range(1,6)
        dip_locations = [2+x for x in random.sample(dip_indivs, len(dip_indivs))]
        rc.add_locations(input_ids=dip_indivs, locations=dip_locations)
        def f(x):
            if x == 0.0 or x == 3.0:
                return x
            else:
                return (math.floor(x) + math.ceil(x))/2.0
        # round breakpoints so we know what's supposed to happen
        round_left = [f(x) for x in rc.args.edges.left]
        round_right = [f(x) for x in rc.args.edges.right]
        rc.args.edges.set_columns(left=round_left, right=round_right,
                                     parent=rc.args.edges.parent,
                                     child=rc.args.edges.child)

        print(rc.args)
        rc.simplify([4,5])
        print(rc.args)
        # this should remove (3,0) and (3,1) and (0,0) (see above)
        nodes = rc.args.nodes
        self.assertEqual(nodes.num_rows, 10)
        true_times = [3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        true_times.reverse()
        self.assertArrayEqual(nodes.time, true_times)
        # check edges
        edges = rc.args.edges
        self.assertEqual(edges.num_rows, 19)
        true_parents = [4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 9, 9, 9]
        self.assertArrayEqual(true_parents, edges.parent)
        true_children = [1, 2,
                         1, 2,
                         0, 0, 3, 3,
                         0, 3,
                         3, 4,
                         0, 4,
                         5, 6, 6, 7, 8]
        print(true_children)
        print(list(edges.child))
        self.assertArrayEqual(true_children, edges.child)
