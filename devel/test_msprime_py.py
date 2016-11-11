def trees(records):
    M = len(records)
    I = sorted(range(M), key=lambda j: (records[j].left, records[j].time))
    O = sorted(range(M), key=lambda j: (records[j].right, -records[j].time))
    pi = [-1 for j in range(max(r.node for r in records) + 1)]
    chi = [[] for j in range(max(r.node for r in records) + 1)]
    j = 0
    k = 0
    while j < M:
        x = records[I[j]].left
        while records[O[k]].right == x:
            h = O[k]
            print("\tout:", records[h])
            chi[records[h].node] = []
            for q in records[h].children:
                pi[q] = -1
            k += 1
        while j < M and records[I[j]].left == x:
            h = I[j]
            print("\tin:", records[h])
            chi[records[h].node] = records[h].children
            for q in records[h].children:
                pi[q] = records[h].node
            j += 1
        yield pi, chi

##
# Vanilla example #1
##

#                                                
# 1.0          6                                                                                                  
# 0.7         / \                                                                     5                            
#            /   \                                                                   / \                           
# 0.5       /     4                           4                                     /   4                                     
#          /     / \                         / \                                   /   / \                             
# 0.4     /     /   \                       /   3                                 /   /   \                             
#        /     /     \                     /   / \                               /   /     \                             
# 0.0   0     1       2                   1   0   2                             0   1       2                           
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)
#
#
# Note that 2 inherits from 3 on the whole chromosome but 3 is only specified on the middle interval.

records = [ 
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=3,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(1,  2),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(1,  3),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(1,  2),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=5,  children=(0,  4),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=6,  children=(0,  4),  time=1.0,  population=0)]

tr_it = trees(records)

next_st = next(tr_it)[0]
assert( all( x==y for x,y in zip(next_st,[6,4,4,-1,6,-1,-1]) ) )
next_st = next(tr_it)[0]
assert( all( x==y for x,y in zip(next_st,[3,4,3,4,-1,-1,-1]) ) )
next_st = next(tr_it)[0]
assert( all( x==y for x,y in zip(next_st,[5,4,4,-1,5,-1,-1]) ) )


##
# With partial dangling tip
##

#
# 1.0             7
# 0.7            / \                                                                     6
#               /   \                                                                   / \
# 0.5          /     5                           5                                     /   5
#             /     / \                         / \                                   /   / \
# 0.4        /     /   4                       /   4                                 /   /   4
#           /     /   / \                     /   / \                               /   /   / \
#          /     /   3   \                   /   /   \                             /   /   3   \
#         /     /         \                 /   /     \                           /   /         \
# 0.0    0     1           2               1   0       2                         0   1           2
#
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)


dangling_records = [
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(2,  3),  time=0.4,  population=0),  # left seg for 3
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(2,  3),  time=0.4,  population=0),  # right seg for 3
    msprime.CoalescenceRecord(left=0.0,  right=1.0,  node=5,  children=(1,  4),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=6,  children=(0,  5),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=7,  children=(0,  5),  time=1.0,  population=0),
    ]

dang_it = trees(dangling_records)

next_st = next(dang_it)[0]
assert( all( x==y for x,y in zip(next_st,[7,5,4,4,5,7,-1,-1]) ) )
next_st = next(dang_it)[0]
assert( all( x==y for x,y in zip(next_st,[4,5,4,-1,5,-1,-1,-1]) ) )
next_st = next(dang_it)[0]
assert( all( x==y for x,y in zip(next_st,[6,5,4,4,5,6,-1,-1]) ) )
