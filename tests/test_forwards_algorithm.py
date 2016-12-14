from ftprime import ARGrecorder
import msprime

 
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
    true_trees = [
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
    def f(lparent,rparent,breakpoint,child,btime):
        arg.add_individual(ids[child])
        times[child]=btime
        if breakpoint>0.0:
            arg.add_record(0.0,breakpoint,ids[lparent],(ids[child],),end_time-times[lparent],0)
        if breakpoint<1.0:
            arg.add_record(breakpoint,1.0,ids[rparent],(ids[child],),end_time-times[rparent],0)

    ids = dict( [ (y,3+x) for x,y in enumerate(['a','b','c','d','e','f','g','h','i','j','k']) ] )
    times = {}
    end_time = 6.0
    arg=ARGrecorder()
    # 1. Begin with an individual `a` (and another anonymous one) at `t=0`.
    arg.add_individual(ids['a'])
    times['a'] = 0.0
    # 2. `(a,?,1.0)->b` and `(a,?,1.0)->c` at `t=1`
    f('a','z',1.0,'b',1.0)
    f('a','z',1.0,'c',1.0)
    # 3. `(b,a,0.9)->d` and `(a,c,0.1)->e` and then `a` dies at `t=2`
    f('b','a',0.9,'d',2.0)
    f('a','c',0.1,'e',2.0)
    # 4. `(d,e,0.7)->f` at `t=3`
    f('d','e',0.7,'f',3.0)
    # 5. `(f,d,0.8)->g` and `(e,f,0.2)->h` at `t=4`.
    f('f','d',0.8,'g',4.0)
    f('e','f',0.2,'h',4.0)
    # 6. `(b,g,0.6)->i` and `(g,h,0.5)->j` and `(c,h,0.4)->k` at `t=5`.
    f('b','g',0.6,'i',5.0)
    f('g','h',0.5,'j',5.0)
    f('c','h',0.4,'k',5.0)
    # 7. We sample `i`, `j` and `k`.
    arg.add_samples(samples=[ids[x] for x in ('i','j','k')],
                    times=[end_time-times[x] for x in ('i','j','k')],
                    populations=[0 for x in ('i','j','k')] )
    ts=arg.tree_sequence()
    samples = [ (0,ids['i']), (1,ids['j']), (2,ids['k']) ]
    try:
        for x,y in zip(ts.trees(),true_trees):
            print(x)
            print(y)
            z = dict([ (ids[i],ids[y[i]]) for i in y.keys()] + samples)
            print(z)
            assert(all( [ x.get_parent(k)==z[k] for k in z.keys() ] ))
    except Exception as e:
        print('wrong tree!')
        print(e)


