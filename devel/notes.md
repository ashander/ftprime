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
