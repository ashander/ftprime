import ftprime
import msprime
import six
import unittest

from tests import FtprimeTestCase

class BasicTestCase(FtprimeTestCase):
    """
    Test basic operations.
    """
    nodes = six.StringIO("""\
    id      is_sample   population      time
    0       0           -1              1.00000000000000
    1       1           -1              0.20000000000000
    2       1           -1              0.00000000000000
    """)
    edges = six.StringIO("""\
    id      left            right           parent  child
    0       0.00000000      1.00000000      0       1
    1       0.00000000      1.00000000      0       2
    """)
    init_ts = msprime.load_text(nodes=nodes, edges=edges)
    init_map = {0:1, 1:2}

    def test_init(self):
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        for input_id in self.init_map:
            node_id = self.init_map[input_id]
            self.assertEqual(records.nodes.time[node_id], 
                             self.init_ts.node(node_id).time)
            self.assertEqual(records.node_ids[input_id], node_id)
            self.assertEqual(records.edges.num_rows, self.init_ts.num_edges)

    def test_add_individual(self):
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        records.add_individual(5, 2.0, population=2)
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+1)
        self.assertEqual(records.nodes.num_rows, 4)
        self.assertEqual(records.nodes.time[records.node_ids[5]], 2.0)
        self.assertEqual(records.nodes.population[records.node_ids[5]], 2)
        self.assertRaises(ValueError,
                          records.add_individuals_dbg, (1, ), (1.5, ))

        # multi
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        records.add_individuals([5, 6], [2.0, 2.5], populations=[2, 1])
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        self.assertEqual(records.nodes.num_rows, 5)
        self.assertEqual(records.nodes.time[records.node_ids[5]], 2.0)
        self.assertEqual(records.nodes.time[records.node_ids[6]], 2.5)
        self.assertEqual(records.nodes.population[records.node_ids[5]], 2)
        self.assertRaises(ValueError,
                          records.add_individuals_dbg, (1, 2), (1.5, 2.5))

    def test_add_record(self):
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        records.add_individual(4, 2.0, population=2)
        records.add_individual(5, 2.0, population=2)
        # adding edges should not change number of nodes
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        records.add_record(0.0, 0.5, 0, (4,5))
        records.add_record(0.5, 1.0, 0, (4,))
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        print(records)
        self.assertEqual(records.edges.num_rows, 5)  # initial 2 + 3 added above
        self.assertEqual(records.edges.parent[2], records.node_ids[0])
        self.assertEqual(records.edges.child[2], records.node_ids[4])
        self.assertEqual(records.edges.child[3], records.node_ids[5])
        self.assertEqual(records.edges.child[4], records.node_ids[4])
        # try adding record with parent who doesn't exist
        self.assertRaises(ValueError, records.add_record, 0.0, 0.5, 8, (0,1))

        # multi
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        records.add_individuals([4, 5], [2.0, 2.0], populations=[2, 2])
        # adding edges should not change number of nodes
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        records.add_records([0.0, 0.5], [0.5, 1.0], [0, 0],
                            [(4,5), (4, )])
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        print(records)
        self.assertEqual(records.edges.num_rows, 5)  # initial 2 + 3 added above
        self.assertEqual(records.edges.parent[2], records.node_ids[0])
        self.assertEqual(records.edges.child[2], records.node_ids[4])
        self.assertEqual(records.edges.child[3], records.node_ids[5])
        self.assertEqual(records.edges.child[4], records.node_ids[4])
        # try adding record with parent who doesn't exist
        self.assertRaises(ValueError, records.add_records_dbg, (0.0, ),
                          (0.5, ), (8, ), ((0,1), ))

    def test_update_times(self):
        records_a = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        # check doing update_times along the way doesn't change things
        records_a.update_times()
        records_b = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        for r in (records_a, records_b):
            r.add_individuals([4, 5], [2.0, 2.0], populations=[2,2])
            r.add_records((0.0, 0.5), (0.5, 1.0), (0, 0), ((4,5), (4, )))
        records_a.update_times()
        records_b.update_times()
        self.assertArrayEqual(records_a.nodes.time, records_b.nodes.time)
        # check update_times is idempotent
        records_b.update_times()
        self.assertArrayEqual(records_a.nodes.time, records_b.nodes.time)
        # and check is right answer
        self.assertArrayEqual(records_a.nodes.time, [3, 2.2, 2, 0, 0])

    def test_get_nodes(self):
        # init
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)

        print(records)
        records.add_individuals([4, 5], [2.0, 2.0], populations=[2,2])
        records.add_records((0.0, 0.5), (0.5, 1.0), (0, 0), ((4,5), (4, )))
        self.assertEqual(records.nodes.num_rows, self.init_ts.num_nodes+2)
        self.assertEqual(records.edges.num_rows, 5)
        print(records)
        final_nodes = records.get_nodes(list(self.init_map.keys()) + [4, 5])
        print(final_nodes)
        for input_id in records.node_ids:
            node_id = records.node_ids[input_id]
            assert node_id in final_nodes
            node_id < self.init_ts.num_nodes + 2
        records.check_ids(list(self.init_map.keys()) + [4, 5])

    def test_simplify(self):
        # test that we get the same tree sequence by doing tree_sequence
        # and simplify -> tree_sequence
        records = ftprime.ARGrecorder(ts=self.init_ts, node_ids=self.init_map)
        records.add_individuals((4, 5), (2.0, 2.0), populations=(2,2))
        records.add_records((0.0, 0.5), (0.5, 1.0), (0, 0), ((4,5), (4, )))
        print(records)
        tsa = records.tree_sequence([4, 5])
        print("---------------- sequence a -----------")
        print(tsa.dump_tables())
        records.simplify([4, 5])
        tsb = records.tree_sequence([4, 5])
        print("---------------- sequence b -----------")
        print(tsb.dump_tables())
        self.check_trees(tsa, tsb)

    def test_simplify2(self):
        # test that nonsensical sequence_length gets caught
        self.assertRaises(ValueError, ftprime.ARGrecorder, ts=self.init_ts, 
                          node_ids=self.init_map, sequence_length=0.5)



class ExplicitTestCase(FtprimeTestCase):
    """
    An explicit test case.

    With `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:

    1. Begin with an individual `a` (and another anonymous one) at `t=0`.
    2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
    3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
    4. `(d,e,0.7)->f` at `t=3`
    5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
    6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(c,h,0.4)->k` at `t=5`.
    7. We sample `i`, `j` and `k`.


    Here are the trees:
    ```
    t                  |              |              |             |             |             |             |             |             |            
                                                                                                                                                      
    0       --a--      |     --a--    |     --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--   
           /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
    1     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b   |   c 
          |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
    2     | d   e |    |   | d   e |  |   | d   e |  |  | d   e |  |  | d   e    |  | d   e    |    d   e    |    d   e    |    d   e    |    d   e   
          | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
    3     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |    | f      |    | f     
          | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
    4     | g   h |    |   | g   h |  |   | g   h |  |  | g   h |  |  | g   h    |  | g   h    |    g   h    |    g   h    |    g   h    |    g   h   
          |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
    5     i   j   k    |   i   j   k  |   i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k 
                                                                                                                                                      
                       |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 
    ```

    and a labeling of the lineages
    ```
    t                  |              |              |             |             |             |             |             |             |            
                                                                                                                                                      
    0       --a--      |     --a--    |     --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--   
           /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
    1     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b   |   c 
          |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
    2     | d   e |    |   | d   e |  |   | d   e |  |  | d   e |  |  | d   e    |  | d   e    |    d   e    |    d   e    |    d   e    |    d   e   
          | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
    3     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |    | f      |    | f     
          | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
    4     | g   h |    |   | g   h |  |   | g   h |  |  | g   h |  |  | g   h    |  | g   h    |    g   h    |    g   h    |    g   h    |    g   h   
          |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
    5     i   j   k    |   i   j   k  |   i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k 
                                                                                                                                                      
                       |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 
    ```
    """

    def f(self, arg, lparent, rparent, breakpoint, child, btime):
        arg.add_individual(self.ids[child], btime)
        if breakpoint>0.0:
            arg.add_record(0.0, breakpoint, self.ids[lparent], (self.ids[child],))
        if breakpoint<1.0:
            arg.add_record(breakpoint, 1.0, self.ids[rparent], (self.ids[child],))

    # the correct tree sequence, unsimplified
    nodes = six.StringIO("""\
    id      is_sample   population      time
    0       0           -1              5.00000000000000  # a
    1       0           -1              4.00000000000000  # b
    2       0           -1              4.00000000000000  # c
    3       0           -1              3.00000000000000  # d
    4       0           -1              3.00000000000000  # e
    5       0           -1              2.00000000000000  # f
    6       0           -1              1.00000000000000  # g
    7       0           -1              1.00000000000000  # h
    8       1           -1              0.00000000000000  # i
    9       1           -1              0.00000000000000  # j
    10      1           -1              0.00000000000000  # k
    """)
    edges = six.StringIO("""\
    id      left            right           parent  child
    0       0.40000000      0.50000000      7       10
    0       0.50000000      1.00000000      7       9
    0       0.50000000      1.00000000      7       10
    0       0.00000000      0.50000000      6       9
    0       0.60000000      1.00000000      6       8
    0       0.00000000      0.20000000      5       6
    0       0.20000000      0.80000000      5       6
    0       0.20000000      0.80000000      5       7
    0       0.80000000      1.00000000      5       7
    0       0.00000000      0.20000000      4       7
    0       0.70000000      1.00000000      4       5
    0       0.00000000      0.70000000      3       5
    0       0.80000000      1.00000000      3       6
    0       0.00000000      0.10000000      2       10
    0       0.10000000      0.40000000      2       4
    0       0.10000000      0.40000000      2       10
    0       0.40000000      1.00000000      2       4
    0       0.00000000      0.60000000      1       3
    0       0.00000000      0.60000000      1       8
    0       0.60000000      0.90000000      1       3
    0       0.00000000      0.10000000      0       1
    0       0.00000000      0.10000000      0       2
    0       0.00000000      0.10000000      0       4
    0       0.10000000      0.90000000      0       1
    0       0.10000000      0.90000000      0       2
    0       0.90000000      1.00000000      0       1
    0       0.90000000      1.00000000      0       2
    0       0.90000000      1.00000000      0       3
    """)
    true_ts = msprime.load_text(nodes=nodes, edges=edges)
    true_tss = true_ts.simplify()

    ids = dict( [ (y,x) for x,y in enumerate(['a','b','c','d','e','f','g','h','i','j','k']) ] )
    true_times = [0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5]
    sample_ids = ('i', 'j', 'k')
    sample_input_ids = [8, 9, 10]

    def test_build_ts(self):
        # build initial tree sequence with just a, b, c
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       0           -1              1.00000000000000
        1       1           -1              0.00000000000000
        2       1           -1              0.00000000000000
        """)
        edges = six.StringIO("""\
        id      left            right           parent  child
        0       0.00000000      1.00000000      0       1
        1       0.00000000      1.00000000      0       2
        """)
        init_ts = msprime.load_text(nodes=nodes, edges=edges)

        first_gen = {self.ids[k] : v for k, v in [('a', 0), ('b', 1), ('c', 2)]}
        arg = ftprime.ARGrecorder(ts=init_ts, node_ids=first_gen, time=1.0)
        # 1. Begin with an individual `a` (and another anonymous one) at `t=0`.
        # taken care of in init_ts
        # arg.add_individual(self.ids['a'], 0.0)
        # # 2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
        # self.f(arg, 'a', 'z', 1.0, 'b', 1.0)
        # self.f(arg, 'a', 'z', 1.0, 'c', 1.0)
        # 3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
        self.f(arg, 'b', 'a', 0.9, 'd', 2.0)
        self.f(arg, 'a', 'c', 0.1, 'e', 2.0)
        # 4. `(d,e,0.7)->f` at `t=3`
        self.f(arg, 'd', 'e', 0.7, 'f', 3.0)
        # 5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
        self.f(arg, 'f', 'd', 0.8, 'g', 4.0)
        self.f(arg, 'e', 'f', 0.2, 'h', 4.0)
        # 6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(c,h,0.4)->k` at `t=5`.
        self.f(arg, 'b', 'g', 0.6, 'i', 5.0)
        self.f(arg, 'g', 'h', 0.5, 'j', 5.0)
        self.f(arg, 'c', 'h', 0.4, 'k', 5.0)
        # 7. We sample `i`, `j` and `k`.
        arg.mark_samples(samples=self.sample_input_ids)
        arg.update_times()

        arg_ids = {k:arg.node_ids[self.ids[k]] for k in self.ids}
        self.assertEqual(arg.nodes.num_rows, len(self.ids))
        self.assertEqual(arg.max_time, 5.0)
        for x in self.ids:
            self.assertEqual(arg.nodes.time[arg_ids[x]], 5.0 - self.true_times[self.ids[x]])
            if x in self.sample_ids:
                self.assertEqual(arg.nodes.flags[arg_ids[x]], msprime.NODE_IS_SAMPLE)
            else:
                self.assertEqual(arg.nodes.flags[arg_ids[x]], 0)

        tss = arg.tree_sequence(self.sample_input_ids)

        self.check_trees(tss, self.true_tss)

    def test_node_times_stable(self):
        # build initial tree sequence with just a, b, c
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       0           -1              1.00000000000000
        1       1           -1              0.00000000000000
        2       1           -1              0.00000000000000
        """)
        edges = six.StringIO("""\
        id      left            right           parent  child
        0       0.00000000      1.00000000      0       1
        1       0.00000000      1.00000000      0       2
        """)
        init_ts = msprime.load_text(nodes=nodes, edges=edges)
        first_gen = {self.ids[k] : v for k, v in [('a', 0), ('b', 1), ('c', 2)]}
        arg = ftprime.ARGrecorder(ts=init_ts, node_ids=first_gen, time=1.0)
        self.f(arg, 'b', 'a', 0.9, 'd', 2.0)
        self.f(arg, 'a', 'c', 0.1, 'e', 2.0)
        self.f(arg, 'd', 'e', 0.7, 'f', 3.0)
        self.f(arg, 'f', 'd', 0.8, 'g', 4.0)
        self.f(arg, 'e', 'f', 0.2, 'h', 4.0)
        self.f(arg, 'b', 'g', 0.6, 'i', 5.0)
        self.f(arg, 'g', 'h', 0.5, 'j', 5.0)
        self.f(arg, 'c', 'h', 0.4, 'k', 5.0)
        arg.update_times()
        node_times = {u:arg.nodes.time[arg.node_ids[u]] for u in arg.node_ids}
        print(arg)
        arg.simplify(self.sample_input_ids)
        print(arg)
        new_node_times = {u:arg.nodes.time[arg.node_ids[u]] for u in arg.node_ids}
        for u in self.sample_input_ids:
            self.assertEqual(node_times[u], new_node_times[u])

    @unittest.skip
    def test_intermediate_simplify(self):
        # build initial tree sequence with just a, b, c
        nodes = six.StringIO("""\
        id      is_sample   population      time
        0       0           -1              1.00000000000000
        1       1           -1              0.00000000000000
        2       1           -1              0.00000000000000
        """)
        edges = six.StringIO("""\
        id      left            right           parent  children
        0       0.00000000      1.00000000      0       1,2
        """)
        init_ts = msprime.load_text(nodes=nodes, edges=edges)

        first_gen = {self.ids[k] : v for k, v in [('a', 0), ('b', 1), ('c', 2)]}
        arg = ftprime.ARGrecorder(ts=init_ts, node_ids=first_gen, time=1.0)
        self.f(arg, 'b', 'a', 0.9, 'd', 2.0)
        self.f(arg, 'a', 'c', 0.1, 'e', 2.0)
        self.f(arg, 'd', 'e', 0.7, 'f', 3.0)
        self.f(arg, 'f', 'd', 0.8, 'g', 4.0)
        # simplify
        print(arg)
        arg.simplify(samples=[self.ids[u] for u in ['b', 'c', 'e', 'f', 'g']])
        print(arg)
        self.f(arg, 'e', 'f', 0.2, 'h', 4.0)
        self.f(arg, 'b', 'g', 0.6, 'i', 5.0)
        self.f(arg, 'g', 'h', 0.5, 'j', 5.0)
        self.f(arg, 'c', 'h', 0.4, 'k', 5.0)
        print(arg)
        tss = arg.tree_sequence(self.sample_input_ids)
        self.check_trees(tss, self.true_tss)
