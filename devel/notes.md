Lets apply alg T to lines 166-174 of test, listed here


    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=4,  children=(2,  3),  time=0.4,  population=0),  # 1
    msprime.CoalescenceRecord(left=0.2,  right=0.8,  node=4,  children=(0,  2),  time=0.4,  population=0),  # 2
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=4,  children=(2,  3),  time=0.4,  population=0),  # 3
    msprime.CoalescenceRecord(left=0.0,  right=1.0,  node=5,  children=(1,  4),  time=0.5,  population=0),  # 4
    msprime.CoalescenceRecord(left=0.8,  right=1.0,  node=6,  children=(0,  5),  time=0.7,  population=0),  # 5
    msprime.CoalescenceRecord(left=0.0,  right=0.2,  node=7,  children=(0,  5),  time=1.0,  population=0),  # 6

The ordered index vectors would be (indexing starting at 1)

`I = [1, 4, 6, 2, 3, 5]` and `O = [6, 1, 2, 5, 4, 3]`.

1.  initialize by setting `pi_j = 0 for j in 1 to 7 + 1`. So `pi = 00000000` and `j =  1, k = 1, x= 0.0`. Also the number of records is `M = 6`
2. `Insert` for x = 0.0?
    -   j = 1; h=1, yes. pi =  00440000
    -   j = 2; h=4, yes. pi =  05445000
    -   j = 3; h=6, yes. pi =  75445700
    -   j = 4; h=2, no. x=0.0 =/= l_2 =0.2
3. `Visit` pi= 75445700, which matches A below although it's unclear 7 is the root. j = 4 < M = 6 so set x = 0.2
4. `Remove` for x = 0.2?
    -   k = 1; h=6, yes. pi = 05445000
    -   k = 2; h=1, yes. pi = 05005000
    -   k = 3; h=2, no. x=0.2 =/= r_2 = 0.8
5. `Insert` for x = 0.2?
    -   j = 4; h=2, yes. pi = 45405000
    -   j = 5; h=3, no. x=0.2 =/= l_3 = 0.8
6. `Visit` pi= 4540500 which matches B although it's unclear that 5 is the root. j= 5 < M=6 so set x = 0.8
7. `Remove` for x = 0.8?
    -   k = 3; h=2, yes. pi = 05005000
    -   k = 4; h=5, no. x=0.8 =/= r_5 = 1.0
8. `Insert` for x = 0.8?
    -   j = 5; h=3, yes. pi = 05445000
    -   j = 6; h=5, yes. pi = 65445600
    -   j = 7; j > M = 6 so goto visit
9. `Visit` pi= 65445600 which matches C although it's unclear that 6 is the root. j= 7 > M so terminate.

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
    #          A(0.0, 0.2),                   B(0.2, 0.8),                             C(0.8, 1.0)


# output from `trees.py`

lists the records in and out and prints each tree visited (in two ways):

    in: CoalescenceRecord(left=0.0, right=0.2, node=4, children=(2, 3), time=0.4, population=0)
    in: CoalescenceRecord(left=0.0, right=1.0, node=5, children=(1, 4), time=0.5, population=0)
    in: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(0, 5), time=1.0, population=0)
    ([7, 5, 4, 4, 5, 7, -1, -1], [[], [], [], [], (2, 3), (1, 4), [], (0, 5)])
    out: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(0, 5), time=1.0, population=0)
    out: CoalescenceRecord(left=0.0, right=0.2, node=4, children=(2, 3), time=0.4, population=0)
    in: CoalescenceRecord(left=0.2, right=0.8, node=4, children=(0, 2), time=0.4, population=0)
    ([4, 5, 4, -1, 5, -1, -1, -1], [[], [], [], [], (0, 2), (1, 4), [], []])
    out: CoalescenceRecord(left=0.2, right=0.8, node=4, children=(0, 2), time=0.4, population=0)
    in: CoalescenceRecord(left=0.8, right=1.0, node=4, children=(2, 3), time=0.4, population=0)
    in: CoalescenceRecord(left=0.8, right=1.0, node=6, children=(0, 5), time=0.7, population=0)
    ([6, 5, 4, 4, 5, 6, -1, -1], [[], [], [], [], (2, 3), (1, 4), (0, 5), []])

# The error in c

Per Peter's sleuthing, the low level code throws an error with these records
that relates to the in_count vs out_count

    #     if (first_tree) {
    #     } else {
    #         if (in_count != out_count) {
    #             ret = MSP_ERR_BAD_COALESCENCE_RECORDS_12;
    #             goto out;
    #         }

The [in_count](https://github.com/ashander/msprime/blob/acd9e3aad2fc828115dd2e6ba5a28b0b097abad5/lib/tree_sequence.c#L2183) and
[out_count](https://github.com/ashander/msprime/blob/acd9e3aad2fc828115dd2e6ba5a28b0b097abad5/lib/tree_sequence.c#L2158) 
both accumulate the number of _children_ (less 1 for some reason) in or out when adding or removing a set of records. 

The above snippet appears just before 'returning' the tree and seems to be
a check that, at any given point after the first tree, the cumulative number of
children (beyond the first child) added or removed is
equal.
For binary trees (like we have) this is just that the number of records in and
out match.

This clearly isn't the case for the example worked through in this document.

But! It is the case both for the first simple example in `test_msprime.py` and
the final example that 'works' (but maybe has incorrect records)
