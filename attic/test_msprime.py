import msprime
import _msprime
from trees import trees

## In this file, testing of:
# Can have partial tips?  Yes, even if the largest ID is not an internal node, but can't plot them.
# Can have partial tips, if numbered increasing towards root(s)?  Yes.
# Can have unsampled tips, if numbered increasing towards root(s)?  Yes.
# Can have coalescence records with the same parent referring to different times?  No.
# Can have single-offspring coalescence records?  Yes.


# Requirements: (from msprime/lib/tree_sequence.c)
#   Input data must be time sorted.
#   Number of children must be at least 2.
#   Children are non-null and in ascending order
#   'left's must come before 'right's

# msprime.simulate(sample_size=None,
#   Ne=1,
#   length=None,
#   recombination_rate=None,
#   recombination_map=None,
#   mutation_rate=None,
#   population_configurations=None,
#   migration_matrix=None,
#   demographic_events=[],
#   samples=None,
#   random_seed=None,
#   num_replicates=None)

print("Simulate records.")

ts = msprime.simulate( sample_size=3, recombination_rate=1.0, random_seed=42 )

# >>> [x for x in ts.records()]
#[ CoalescenceRecord(left=0.1912159270586483,  right=0.8521429346530099,  node=3,  children=(0,  2),  time=0.40566159044942235,  population=0),
#  CoalescenceRecord(left=0.0,                 right=0.1912159270586483,  node=4,  children=(1,  2),  time=0.44077247376386364,  population=0),
#  CoalescenceRecord(left=0.1912159270586483,  right=0.8521429346530099,  node=4,  children=(1,  3),  time=0.44077247376386364,  population=0),
#  CoalescenceRecord(left=0.8521429346530099,  right=1.0,                 node=4,  children=(1,  2),  time=0.44077247376386364,  population=0),
#  CoalescenceRecord(left=0.8521429346530099,  right=1.0,                 node=5,  children=(0,  4),  time=0.7114579294844481,   population=0),
#  CoalescenceRecord(left=0.0,                 right=0.1912159270586483,  node=6,  children=(0,  4),  time=2.8161856234286375,   population=0)]

# >>> [ x.draw(path="tree_{}.svg".format(k)) for k,x in enumerate(ts.trees()) ]
# >>> [ x.get_interval() for x in ts.trees() ]
# Marginal trees are:
#
# 2.8          6
# 0.7         / \                                                                     5
#            /   \                                                                   / \
# 0.44      /     4                           4                                     /   4
#          /     / \                         / \                                   /   / \
# 0.4     /     /   \                       /   3                                 /   /   \
#        /     /     \                     /   / \                               /   /     \
# 0.0   0     1       2                   1   0   2                             0   1       2
#
# (0.0, 0.1912159270586483), (0.1912159270586483, 0.8521429346530099), (0.8521429346530099, 1.0)
#

###
print("Read records into msprime:")

records_0 = [
    msprime.CoalescenceRecord(left=0.1912159270586483,  right=0.8521429346530099,  node=3,  children=(0,  2),  time=0.40566159044942235,  population=0),
    msprime.CoalescenceRecord(left=0.0,                 right=0.1912159270586483,  node=4,  children=(1,  2),  time=0.44077247376386364,  population=0),
    msprime.CoalescenceRecord(left=0.1912159270586483,  right=0.8521429346530099,  node=4,  children=(1,  3),  time=0.44077247376386364,  population=0),
    msprime.CoalescenceRecord(left=0.8521429346530099,  right=1.0,                 node=4,  children=(1,  2),  time=0.44077247376386364,  population=0),
    msprime.CoalescenceRecord(left=0.8521429346530099,  right=1.0,                 node=5,  children=(0,  4),  time=0.7114579294844481,   population=0),
    msprime.CoalescenceRecord(left=0.0,                 right=0.1912159270586483,  node=6,  children=(0,  4),  time=2.8161856234286375,   population=0)]


ll_ts_0 = _msprime.TreeSequence()
ll_ts_0.load_records(records_0)
ts_0 = msprime.TreeSequence(ll_ts_0)

[ x==y for x,y in zip( ts.records(), ts_0.records() ) ]

[ x==y for x,y in zip( ts.trees(), ts_0.trees() ) ]

###
print("Adjust records")

# Round down the times:
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
# Note that one possibility is that 2 inherits from 3 on the whole chromosome but 3 is only specified on the middle interval.

###
print("Read records into msprime:")
print("--------------------------")

records = [
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=3,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(1,  2),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(1,  3),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(1,  2),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=5,  children=(0,  4),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=6,  children=(0,  4),  time=1.0,  population=0)]

true_trees = [ {0: 6, 1: 4, 2: 4, 3:-1, 4: 6, 5:-1, 6:-1},
               {0: 3, 1: 4, 2: 3, 3: 4, 4:-1, 5:-1, 6:-1},
               {0: 5, 1: 4, 2: 4, 3:-1, 4: 5, 5:-1, 6:-1} ]

ll_ts = _msprime.TreeSequence()
ll_ts.load_records(records)
ts = msprime.TreeSequence(ll_ts)

try:
    [ x.draw(path="tree_{}.svg".format(k)) for k,x in enumerate(ts.trees()) ]
    print('plotting works in c')
except:
    pass

try:
    for x,y in zip(ts.trees(),true_trees):
        assert(all( [ x.get_parent(k)==y[k] for k in y.keys() ] ))
        print(x)
        print(y)
except Exception as e:
    print('wrong tree!')
    print(e)


try:
    for x,(y,z) in zip(true_trees,trees(list(records))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python fails.')
    print(e)


########
print("Can have partial tips?")
print("--------------------------")
print("Yes, except plotting:")

# Yes, but plotting fails:

# Try inserting the invisible 3 in the left and right trees
# by adding a ghost offspring:
#
# 1.0          6
# 0.7         / \                                                                     5
#            /   \                                                                   / \
# 0.5       /     4                           4                                     /   4
#          /     / \                         / \                                   /   / \
# 0.4     /     /   3                       /   3                                 /   /   3
#        /     /   / \                     /   / \                               /   /   / \
#       /     /   7   \                   /   /   \                             /   /   7   \
#      /     /         \                 /   /     \                           /   /         \
# 0.0 0     1           2               1   0       2                         0   1           2
#
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)


records_1 = [
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=3,  children=(2,  7),  time=0.4,  population=0),  # left seg for 7
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=3,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=3,  children=(2,  7),  time=0.4,  population=0),  # right seg for 7
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(1,  3),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(1,  3),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(1,  3),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=5,  children=(0,  4),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=6,  children=(0,  4),  time=1.0,  population=0),
    ]

true_trees_1 = [ {0: 6, 1: 4, 2: 3, 3: 4, 4: 6, 5:-1, 6:-1, 7: 3},
                 {0: 3, 1: 4, 2: 3, 3: 4, 4:-1, 5:-1, 6:-1, 7:-1},
                 {0: 5, 1: 4, 2: 3, 3: 4, 4: 5, 5:-1, 6:-1, 7: 3}]

ll_ts_1 = _msprime.TreeSequence()
ll_ts_1.load_records(records_1)
ts_1 = msprime.TreeSequence(ll_ts_1)

try:
    [ x.draw(path="tree_{}.svg".format(k)) for k,x in enumerate(ts_1.trees()) ]
except Exception as e:
    print('plotting trees fails.')
    print(e)

try:
    for x,(y,z) in zip(true_trees_1,trees(list(records_1))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python fails.')
    print(e)

try:
    for x,y in zip(true_trees_1,ts_1.trees()):
        print(x)
        print(y)
        assert( all( [ x[k]==y.parent(k) for k in x.keys() ] ) )
        pass
except Exception as e:
    print('msprime does not agree.')
    print(e)




########
print("Can have partial tips, if numbered increasing towards root(s)?")
print("--------------------------")

# Yes!


# Same, but renumbered:
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


records_2 = [
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(2,  3),  time=0.4,  population=0),  # left seg for 3
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(2,  3),  time=0.4,  population=0),  # right seg for 3
    msprime.CoalescenceRecord(left=0.0,  right=1.0,  node=5,  children=(1,  4),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=6,  children=(0,  5),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=7,  children=(0,  5),  time=1.0,  population=0),
    ]

true_trees_2 = [ {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6:-1, 7:-1},
               {0: 4, 1: 5, 2: 4, 3:-1, 4: 5, 5:-1, 6:-1, 7:-1},
               {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6:-1, 7:-1}]

ll_ts_2 = _msprime.TreeSequence()
ll_ts_2.load_records(records_2)
ts_2 = msprime.TreeSequence(ll_ts_2)

try:
    for x,(y,z) in zip(true_trees_2,trees(list(records_2))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python gets wrong tree.')
    print(e)


try:
    [ x.draw(path="new_tree_{}.svg".format(k)) for k,x in enumerate(ts_2.trees()) ]
except Exception as e:
    print('fails in c with')
    print('  ', e)

try:
    for x,y in zip(true_trees_2,ts_2.trees()):
        print(x)
        print(y)
        assert( all( [ x[k]==y.parent(k) for k in x.keys() ] ) )
        pass
except Exception as e:
    print('msprime does not agree.')
    print(e)



########
print("Can have unsampled tips, if numbered increasing towards root(s)?")
print("--------------------------")

# Yes!


# Adding whole genome info for phantom
#
# 1.0             8
# 0.7            / \                                                                     7
#               /   \                                                                   / \
# 0.5          /     6                           6                                     /   6
#             /     / \                         / \                                   /   / \
# 0.4        /     /   5                       /   5                                 /   /   5
#           /     /   / \                     /   / \                               /   /   / \
# 0.2      /     /   <   \                   /   4   \                             /   /   <   \
#         /     /     \   \                 /   / \   \                           /   /     \   \
# 0.0    0     1       3   2               1   0   3   2                         0   1       3   2
#
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)


records_3 = [
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  3),  time=0.2,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=5,  children=(2,  3),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=5,  children=(2,  4),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=5,  children=(2,  3),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=1.0,  node=6,  children=(1,  5),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=7,  children=(0,  6),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=8,  children=(0,  6),  time=1.0,  population=0),
    ]

ll_ts_3 = _msprime.TreeSequence()
ll_ts_3.load_records(records_3)
ts_3 = msprime.TreeSequence(ll_ts_3)
for t in trees(list(ts_3.records())):
    print(t)
    pass

true_trees_3 = [ {0: 8, 1: 6, 2: 5, 3: 5, 4:-1, 5: 6, 6: 8, 7:-1, 8:-1},
               {0: 4, 1: 6, 2: 5, 3: 4, 4: 5, 5: 6, 6:-1, 7:-1, 8:-1},
               {0: 7, 1: 6, 2: 5, 3: 5, 4:-1, 5: 6, 6: 7, 7:-1, 8:-1}]

try:
    for x,y in zip(ts_3.trees(),true_trees_3):
        assert(all( [ x.get_parent(k)==y[k] for k in y.keys() ] ))
        print(x)
        print(y)
except Exception as e:
    print('wrong tree!')
    print(e)

## Works!

try:
    [ x.draw(path="new_tree_{}.svg".format(k)) for k,x in enumerate(ts_3.trees()) ]
except Exception as e:
    print('drawing trees fails')
    print(e)

try:
    for x,(y,z) in zip(true_trees_3,trees(list(records_3))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python gets wrong tree.')
    print(e)




########
print("Can have coalescence records with the same parent referring to different times?")
print("--------------------------")
print("No, as expected:")

# No.


#
# 1.0             7
# 0.7            / \                                                                     6
# 0.6           /   \   ....................    5     ..............................    / \
# 0.5          /     5                         / \                                     /   5
#             /     / \                       /   \                                   /   / \
# 0.4        /     /   4                     /     4                                 /   /   4
#           /     /   / \                   /     / \                               /   /   / \
#          /     /   3   \                 /     /   \                             /   /   3   \
#         /     /         \               /     /     \                           /   /         \
# 0.0    0     1           2             1     0       2                         0   1           2
#
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)


records_4 = [
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(2,  3),  time=0.4,  population=0),  # left seg for 3
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(2,  3),  time=0.4,  population=0),  # right seg for 3
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=5,  children=(1,  4),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=5,  children=(1,  4),  time=0.6,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=5,  children=(1,  4),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=6,  children=(0,  5),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=7,  children=(0,  5),  time=1.0,  population=0),
    ]

true_trees_4 = [ {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 7, 6:-1, 7:-1},
               {0: 4, 1: 5, 2: 4, 3:-1, 4: 5, 5:-1, 6:-1, 7:-1},
               {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6:-1, 7:-1}]


ll_ts_4 = _msprime.TreeSequence()
try:
    ll_ts_4.load_records(records_4)
    ts_4 = msprime.TreeSequence(ll_ts_4)
except Exception as e:
    print('loading trees fails')
    print(e)

# Unsuprisingly, fails at lib/tree_sequence.c :
#         } else if (self->trees.nodes.time[node] != records[j].time) {
# with:
#     Times associated with nodes not consistent between records

try:
    for x,(y,z) in zip(true_trees_4,trees(list(records_4))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python gets wrong tree.')
    print(e)




########
print("Can have single-offspring coalescence records?")
print("--------------------------")

# Yes!


#
# 1.0             7
# 0.7            / 6                                                                     6
#               /   \                                                                   / \
# 0.5          /     5                           5                                     /   5
#             /     / \                         / \                                   /   / \
# 0.4        /     /   4                       /   4                                 /   /   4
# 0.3       /     /   / \                     /   / \                               /   /   / \
#          /     /   3   \                   /   /   \                             /   /   3   \
#         /     /         \                 /   /     \                           /   /         \
# 0.0    0     1           2               1   0       2                         0   1           2
#
#          (0.0, 0.2),                   (0.2, 0.8),                             (0.8, 1.0)


records_5 = [
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(2,  3),  time=0.4,  population=0),  # left seg for 3
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  2),  time=0.4,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(2,  3),  time=0.4,  population=0),  # right seg for 3
    msprime.CoalescenceRecord(left=0.0,  right=1.0,  node=5,  children=(1,  4),  time=0.5,  population=0),
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=6,  children=(0,  5),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=6,  children=(5,   ),  time=0.7,  population=0),
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=7,  children=(0,  6),  time=1.0,  population=0),
    ]

true_trees_5 = [ {0: 7, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6: 7, 7:-1},
               {0: 4, 1: 5, 2: 4, 3:-1, 4: 5, 5:-1, 6:-1, 7:-1},
               {0: 6, 1: 5, 2: 4, 3: 4, 4: 5, 5: 6, 6:-1, 7:-1}]


ll_ts_5 = _msprime.TreeSequence()
try :
    ll_ts_5.load_records(records_5)
    ts_5 = msprime.TreeSequence(ll_ts_5)
except Exception as e:
    print(e)


try:
    for x,(y,z) in zip(true_trees_5,trees(list(records_5))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python gets wrong tree.')
    print(e)



try:
    for x,y in zip(ts_5.trees(),true_trees_5):
        assert(all( [ x.get_parent(k)==y[k] for k in y.keys() ] ))
        print(x)
        print(y)
except Exception as e:
    print('wrong tree!')
    print(e)


########
print("A more complex test with single offspring")
print("--------------------------")

# With `(i,j,x)->k` denoting that individual `k` inherits from `i` on `[0,x)` and from `j` on `[x,1)`:
# 1. Begin with an individual `3` (and another anonymous one) at `t=0`.
# 2. `(3,?,1.0)->4` and `(3,?,1.0)->5` at `t=1`
# 3. `(4,3,0.9)->6` and `(3,5,0.1)->7` and then `3` dies at `t=2`
# 4. `(6,7,0.7)->8` at `t=3`
# 5. `(8,6,0.8)->9` and `(7,8,0.2)->10` at `t=4`.
# 6. `(3,9,0.6)->0` and `(9,10,0.5)->1` and `(10,4,0.4)->2` at `t=5`.
# 7. We sample `0`, `1`, and `2`.
# Here are the trees:
# t                  |              |              |             |             |             |             |             |             |            
#                                                                                                                                                   
# 0       --3--      |     --3--    |     --3--    |    --3--    |    --3--    |    --3--    |    --3--    |    --3--    |    --3--    |    --3--   
#        /  |  \     |    /  |  \   |    /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /     \   |   /  |  \  
# 1     4   |   5    |   4   |   5  |   4       5  |  4       5  |  4       5  |  4       5  |  4       5  |  4       5  |  4       5  |  4   |   5 
#       |\ / \ /|    |   |\   \     |   |\     /   |  |\     /   |  |\     /|  |  |\     /|  |   \     /|  |   \     /|  |   \     /|  |     /   /| 
# 2     | 6   7 |    |   | 6   7    |   | 6   7    |  | 6   7    |  | 6   7 |  |  | 6   7 |  |    6   7 |  |    6   7 |  |    6   7 |  |    6   7 | 
#       | |\ /| |    |   |  \  |    |   |  \  |    |  |  \       |  |  \    |  |  |  \    |  |     \    |  |       /  |  |    |  /  |  |    |  /  | 
# 3     | | 8 | |    |   |   8 |    |   |   8 |    |  |   8      |  |   8   |  |  |   8   |  |      8   |  |      8   |  |    | 8   |  |    | 8   | 
#       | |/ \| |    |   |  /  |    |   |  /  |    |  |  / \     |  |  / \  |  |  |  / \  |  |     / \  |  |     / \  |  |    |  \  |  |    |  \  | 
# 4     | 9  10 |    |   | 9  10    |   | 9  10    |  | 9  10    |  | 9  10 |  |  | 9  10 |  |    9  10 |  |    9  10 |  |    9  10 |  |    9  10 | 
#       |/ \ / \|    |   |  \   \   |   |  \   \   |  |  \   \   |  |  \    |  |  |    /  |  |   /   /  |  |   /   /  |  |   /   /  |  |   /   /  | 
# 5     0   1   2    |   0   1   2  |   0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2  |  0   1   2 
#                                                                                                                                                   
#                    |   0.0 - 0.1  |   0.1 - 0.2  |  0.2 - 0.4  |  0.4 - 0.5  |  0.5 - 0.6  |  0.6 - 0.7  |  0.7 - 0.8  |  0.8 - 0.9  |  0.9 - 1.0 


true_trees_6 = [
        { 0:4, 1:9, 2:10, 3:-1, 4:3, 5:3, 6:4, 7:3, 8:6, 9:8, 10:7 },
        { 0:4, 1:9, 2:10, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:6, 9:8, 10:7 },
        { 0:4, 1:9, 2:10, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:6, 9:8, 10:8 },
        { 0:4, 1:9, 2:5, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:6, 9:8, 10:8 },
        { 0:4, 1:10, 2:5, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:6, 9:8, 10:8 },
        { 0:9, 1:10, 2:5, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:6, 9:8, 10:8 },
        { 0:9, 1:10, 2:5, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:7, 9:8, 10:8 },
        { 0:9, 1:10, 2:5, 3:-1, 4:3, 5:3, 6:4, 7:5, 8:7, 9:6, 10:8 },
        { 0:9, 1:10, 2:5, 3:-1, 4:3, 5:3, 6:3, 7:5, 8:7, 9:6, 10:8 }
    ]

records_6 = [ 
       msprime.CoalescenceRecord( left=0.5, right=1.0,  node=10, children=   (1,),  time=5.0-4.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.4,  node=10, children=   (2,),  time=5.0-4.0, population=0 ),
       msprime.CoalescenceRecord( left=0.6, right=1.0,  node=9,  children=   (0,),  time=5.0-4.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.5,  node=9,  children=   (1,),  time=5.0-4.0, population=0 ),
       msprime.CoalescenceRecord( left=0.8, right=1.0,  node=8,  children=  (10,),  time=5.0-3.0, population=0 ),
       msprime.CoalescenceRecord( left=0.2, right=0.8,  node=8,  children= (9,10),  time=5.0-3.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.2,  node=8,  children=   (9,),  time=5.0-3.0, population=0 ),
       msprime.CoalescenceRecord( left=0.7, right=1.0,  node=7,  children=   (8,),  time=5.0-2.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.2,  node=7,  children=  (10,),  time=5.0-2.0, population=0 ),
       msprime.CoalescenceRecord( left=0.8, right=1.0,  node=6,  children=   (9,),  time=5.0-2.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.7,  node=6,  children=   (8,),  time=5.0-2.0, population=0 ),
       msprime.CoalescenceRecord( left=0.4, right=1.0,  node=5,  children=  (2,7),  time=5.0-1.0, population=0 ),
       msprime.CoalescenceRecord( left=0.1, right=0.4,  node=5,  children=   (7,),  time=5.0-1.0, population=0 ),
       msprime.CoalescenceRecord( left=0.6, right=0.9,  node=4,  children=   (6,),  time=5.0-1.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.6,  node=4,  children=  (0,6),  time=5.0-1.0, population=0 ),
       msprime.CoalescenceRecord( left=0.9, right=1.0,  node=3,  children=(4,5,6),  time=5.0-0.0, population=0 ),
       msprime.CoalescenceRecord( left=0.1, right=0.9,  node=3,  children=  (4,5),  time=5.0-0.0, population=0 ),
       msprime.CoalescenceRecord( left=0.0, right=0.1,  node=3,  children=(4,5,7),  time=5.0-0.0, population=0 ),
       ]


ll_ts_6 = _msprime.TreeSequence()
try :
    ll_ts_6.load_records(records_6)
    ts_6 = msprime.TreeSequence(ll_ts_6)
except Exception as e:
    print(e)


try:
    for x,(y,z) in zip(true_trees_6,trees(list(records_6))):
        print(x)
        print(y)
        assert( all( [ x[k]==y[k] for k in range(len(y)) ] ) )
        pass
except Exception as e:
    print('python gets wrong tree.')
    print(e)



try:
    for x,y in zip(ts_6.trees(),true_trees_6):
        assert(all( [ x.get_parent(k)==y[k] for k in y.keys() ] ))
        print(x)
        print(y)
except Exception as e:
    print('wrong tree!')
    print(e)


