import ftprime
import msprime
import pytest
import six


def check_trees(tsa, tsb, npos=20):
    # check trees at a bunch of positions agree
    assert len(tsa.samples()) == len(tsb.samples())
    assert all([a == b for a, b in zip(tsa.samples(), tsb.samples())])
    check_positions = [k/npos for k in range(npos)]
    tsa_t = tsa.trees()
    tsb_t = tsb.trees()
    ta = next(tsa_t)
    tb = next(tsb_t)
    for pos in check_positions:
        while pos > ta.interval[1]:
            ta = next(tsa_t)
        while pos > tb.interval[1]:
            tb = next(tsb_t)
        check_tree(ta, tb, tsa.samples())
    return True

def check_tree(ta, tb, samples):
    for u in samples:
        for v in samples:
            assert ta.mrca(u,v) == tb.mrca(u,v)
    return True


# With `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:
# 
# 1. Begin with an individual `a` (and another anonymous one) at `t=0`.
# 2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
# 3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
# 4. `(d,e,0.7)->f` at `t=3`
# 5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
# 6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(c,h,0.4)->k` at `t=5`.
# 7. We sample `i`, `j` and `k`.
# 
 
# Here are the trees:
# ```
# t                  |              |              |             |             |             |             |             |             |            
#                                                                                                                                                   
# 0       --a--      |     --a--    |     --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--   
#        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
# 1     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b   |   c 
#       |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
# 2     | d   e |    |   | d   e |  |   | d   e |  |  | d   e |  |  | d   e    |  | d   e    |    d   e    |    d   e    |    d   e    |    d   e   
#       | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
# 3     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |    | f      |    | f     
#       | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
# 4     | g   h |    |   | g   h |  |   | g   h |  |  | g   h |  |  | g   h    |  | g   h    |    g   h    |    g   h    |    g   h    |    g   h   
#       |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
# 5     i   j   k    |   i   j   k  |   i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k 
#                                                                                                                                                   
#                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 
# ```
 
# and a labeling of the lineages
# ```
# t                  |              |              |             |             |             |             |             |             |            
#                                                                                                                                                   
# 0       --a--      |     --a--    |     --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--    |    --a--   
#        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
# 1     b   |   c    |   b   |   c  |   b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b       c  |  b   |   c 
#       |\ / \ /|    |   |\   \  |  |   |\     /|  |  |\     /|  |  |\     /   |  |\     /   |   \     /   |   \     /   |   \     /   |     /   /  
# 2     | d   e |    |   | d   e |  |   | d   e |  |  | d   e |  |  | d   e    |  | d   e    |    d   e    |    d   e    |    d   e    |    d   e   
#       | |\ /| |    |   |  \  | |  |   |  \  | |  |  |  \    |  |  |  \       |  |  \       |     \       |       /     |    |  /     |    |  /    
# 3     | | f | |    |   |   f | |  |   |   f | |  |  |   f   |  |  |   f      |  |   f      |      f      |      f      |    | f      |    | f     
#       | |/ \| |    |   |  /  | |  |   |  /  | |  |  |  / \  |  |  |  / \     |  |  / \     |     / \     |     / \     |    |  \     |    |  \    
# 4     | g   h |    |   | g   h |  |   | g   h |  |  | g   h |  |  | g   h    |  | g   h    |    g   h    |    g   h    |    g   h    |    g   h   
#       |/ \ / \|    |   |  \    |  |   |  \    |  |  |  \    |  |  |  \   \   |  |    / \   |   /   / \   |   /   / \   |   /   / \   |   /   / \  
# 5     i   j   k    |   i   j   k  |   i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k  |  i   j   k 
#                                                                                                                                                   
#                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 
# ```

def test_case():
    # build initial tree sequence with just a, b, c
    nodes = six.StringIO("""\
    id      is_sample   population      time
    0       0           -1              1.00000000000000
    1       1           -1              0.00000000000000
    2       1           -1              0.00000000000000
    """)
    edgesets = six.StringIO("""\
    id      left            right           parent  children
    0       0.00000000      1.00000000      0       1,2
    """)
    init_ts = msprime.load_text(nodes=nodes, edgesets=edgesets)

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
    edgesets = six.StringIO("""\
    id      left            right           parent  children
    0       0.40000000      0.50000000      7       10
    0       0.50000000      1.00000000      7       9,10
    0       0.00000000      0.50000000      6       9
    0       0.60000000      1.00000000      6       8
    0       0.00000000      0.20000000      5       6
    0       0.20000000      0.80000000      5       6,7
    0       0.80000000      1.00000000      5       7
    0       0.00000000      0.20000000      4       7
    0       0.70000000      1.00000000      4       5
    0       0.00000000      0.70000000      3       5
    0       0.80000000      1.00000000      3       6
    0       0.00000000      0.10000000      2       10
    0       0.10000000      0.40000000      2       4,10
    0       0.40000000      1.00000000      2       4
    0       0.00000000      0.60000000      1       3,8
    0       0.60000000      0.90000000      1       3
    0       0.00000000      0.10000000      0       1,2,4
    0       0.10000000      0.90000000      0       1,2
    0       0.90000000      1.00000000      0       1,2,3
    """)
    true_ts = msprime.load_text(nodes=nodes, edgesets=edgesets)

    big_trees = [
            { 'b':'a', 'c':'a', 'd':'b', 'e':'a', 'f':'d', 'g':'f', 'h':'e', 'i':'b', 'j':'g', 'k':'c' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'d', 'g':'f', 'h':'e', 'i':'b', 'j':'g', 'k':'c' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'d', 'g':'f', 'h':'f', 'i':'b', 'j':'g', 'k':'c' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'d', 'g':'f', 'h':'f', 'i':'b', 'j':'g', 'k':'h' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'d', 'g':'f', 'h':'f', 'i':'b', 'j':'h', 'k':'h' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'d', 'g':'f', 'h':'f', 'i':'g', 'j':'h', 'k':'h' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'e', 'g':'f', 'h':'f', 'i':'g', 'j':'h', 'k':'h' },
            { 'b':'a', 'c':'a', 'd':'b', 'e':'c', 'f':'e', 'g':'d', 'h':'f', 'i':'g', 'j':'h', 'k':'h' },
            { 'b':'a', 'c':'a', 'd':'a', 'e':'c', 'f':'e', 'g':'d', 'h':'f', 'i':'g', 'j':'h', 'k':'h' },
        ]
    end_time = 6.0
    ids = dict( [ (y,x) for x,y in enumerate(['a','b','c','d','e','f','g','h','i','j','k']) ] )
    true_times = [0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5]
    def f(lparent,rparent,breakpoint,child,btime):
        arg.add_individual(ids[child], btime)
        if breakpoint>0.0:
            arg.add_record(0.0, breakpoint, ids[lparent], (ids[child],))
        if breakpoint<1.0:
            arg.add_record(breakpoint, 1.0, ids[rparent], (ids[child],))

    first_gen = {ids[k] : v for k, v in [('a', 0), ('b', 1), ('c', 2)]}
    arg = ftprime.ARGrecorder(ts=init_ts, node_ids=first_gen, time=1.0)
    # 1. Begin with an individual `a` (and another anonymous one) at `t=0`.
    # taken care of in init_ts
    # arg.add_individual(ids['a'], 0.0)
    # # 2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
    # f('a', 'z', 1.0, 'b', 1.0)
    # f('a', 'z', 1.0, 'c', 1.0)
    # 3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
    f('b', 'a', 0.9, 'd', 2.0)
    f('a', 'c', 0.1, 'e', 2.0)
    # 4. `(d,e,0.7)->f` at `t=3`
    f('d', 'e', 0.7, 'f', 3.0)
    # 5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
    f('f', 'd', 0.8, 'g', 4.0)
    f('e', 'f', 0.2, 'h', 4.0)
    # 6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(c,h,0.4)->k` at `t=5`.
    f('b', 'g', 0.6, 'i', 5.0)
    f('g', 'h', 0.5, 'j', 5.0)
    f('c', 'h', 0.4, 'k', 5.0)
    # 7. We sample `i`, `j` and `k`.
    sample_ids = ('i', 'j', 'k')
    samples = [ids[x] for x in sample_ids]
    arg.mark_samples(samples=samples)
    arg.update_times()

    arg_ids = {k:arg.node_ids[ids[k]] for k in ids}
    assert arg.num_nodes == len(ids)
    assert arg.max_time == 5.0
    assert arg.sequence_length == 1.0
    for x in ids:
        assert arg.nodes.time[arg_ids[x]] == 5.0 - true_times[ids[x]]
        if x in sample_ids:
            assert arg.nodes.flags[arg_ids[x]] == msprime.NODE_IS_SAMPLE
        else:
            assert arg.nodes.flags[arg_ids[x]] == 0

    tss = arg.tree_sequence(samples)
    true_tss = true_ts.simplify()

    assert tss.num_nodes == true_tss.num_nodes

    for x in tss.dump_tables():
        print(x)

    for x in true_tss.dump_tables():
        print(x)

    check_trees(tss, true_tss)
