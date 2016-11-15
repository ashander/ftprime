from ftprime.alg import main
import tempfile
import msprime
import _msprime


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


def parent_dict(pi):
    parent_dict = {}
    for child, parent in enumerate(pi):
        parent_dict[child] = parent
    return parent_dict


population = main()
#  what main() does
#    p = Population(0.0, state={'a': 1, '?': -1})
#    births = [Birth('a', '?', 1.0, 'b'), ]
#    deaths = ['?',]
#    p._onegen(births, deaths)
#
#    births = [Birth('a', 'b', 0.5, 'c'), ]
#    p._onegen(births, deaths)
#
#    births = [Birth('a', 'c', 0.2, 'd'),
#              Birth('c', 'b', 0.6, 'e'), ]
#    p._onegen(births, deaths)
#
# output
#Births [Birth(left='a', right='?', x=1.0, offspring='b')]
#{'b': 3, 'a': 2, '?': -1}
#Births [Birth(left='a', right='b', x=0.5, offspring='c')]
#{'b': 5, 'a': 4, 'c': 6}
#Births [Birth(left='a', right='c', x=0.2, offspring='d'), Birth(left='c', right='b', x=0.6, offspring='e')]
#{'b': 5, 'a': 7, 'd': 9, 'c': 8}
#{'b': 10, 'a': 7, 'e': 11, 'd': 9, 'c': 8}
#--------------------- fin ------------
#CoalescenceRecord(left=0.0, right=1.0, node=1, children=(2, 3), time=1.0, population=0)
#CoalescenceRecord(left=0.0, right=0.5, node=2, children=(4, 6), time=2.0, population=0)
#CoalescenceRecord(left=0.5, right=1.0, node=3, children=(5, 6), time=2.0, population=0)
#CoalescenceRecord(left=0.0, right=0.2, node=4, children=(7, 9), time=3.0, population=0)
#CoalescenceRecord(left=0.6, right=1.0, node=5, children=(10, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.0, right=0.2, node=6, children=(8, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.2, right=0.5, node=6, children=(8, 9, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.5, right=0.6, node=6, children=(8, 9, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.6, right=1.0, node=6, children=(8, 9, 11), time=3.0, population=0)

mapping = population.renumber()
print("\n  ------ feeding a renumbered version into the msprime alg T ------")
print("      the mapping from NEW to OLD:")
print(mapping)
print("\n state:")
print(population.state)
print("(in old numbers:", [(k, mapping[v]) for k, v in population.state.items()])
print("\n")
pytrees = []
for pi, chi in trees(list(population)):
    print("  trees: ", pi, chi)
    pytrees.append(pi)
    t_1_renum = [[mapping[j] for j in els] if len(els) > 0 else [] for els in chi]
    t_1_renum.reverse()
    print("(in old numbers:", t_1_renum)
#
#  ------ feeding a renumbered version into the msprime alg T ------
#      the mapping from NEW to OLD:
#{0: 10, 1: 9, 2: 8, 3: 7, 4: 6, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1}
#
# state:
#{'b': 1, 'a': 4, 'e': 0, 'd': 2, 'c': 3}
#(in old numbers: [('b', 9), ('a', 6), ('e', 10), ('d', 8), ('c', 7)]
#
#
#	in: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(2, 4), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.2, node=5, children=(0, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.5, node=9, children=(5, 7), time=1.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=1.0, node=10, children=(8, 9), time=2.0, population=0)
#  trees:  [5, -1, 7, 5, 7, 9, -1, 9, 10, 10, -1] [[], [], [], [], [], (0, 3), [], (2, 4), [], (5, 7), (8, 9)]
#(in old numbers: [[2, 1], [5, 3], [], [8, 6], [], [10, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(2, 4), time=0.0, population=0)
#	out: CoalescenceRecord(left=0.0, right=0.2, node=5, children=(0, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.2, right=0.5, node=5, children=(0, 2, 3), time=0.0, population=0)
#  trees:  [5, -1, 5, 5, -1, 9, -1, 9, 10, 10, -1] [[], [], [], [], [], (0, 2, 3), [], [], [], (5, 7), (8, 9)]
#(in old numbers: [[2, 1], [5, 3], [], [], [], [10, 8, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.0, right=0.5, node=9, children=(5, 7), time=1.0, population=0)
#	out: CoalescenceRecord(left=0.2, right=0.5, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=0.6, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=1.0, node=8, children=(5, 6), time=1.0, population=0)
#  trees:  [5, -1, 5, 5, -1, 8, 8, -1, 10, 10, -1] [[], [], [], [], [], (0, 2, 3), [], [], (5, 6), [], (8, 9)]
#(in old numbers: [[2, 1], [], [5, 4], [], [], [10, 8, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.5, right=0.6, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=6, children=(0, 1), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=5, children=(0, 2, 3), time=0.0, population=0)
#  trees:  [5, 6, 5, 5, -1, 8, 8, -1, 10, 10, -1] [[], [], [], [], [], (0, 2, 3), (0, 1), [], (5, 6), [], (8, 9)]
#(in old numbers: [[2, 1], [], [5, 4], [], [10, 9], [10, 8, 7], [], [], [], [], []]
# According to my renumbering, we get
#
#3  |       a       |       1             1           1              1        |                    |
#   |      / \      |      / \           / \         / \            / \       |                    |
#   |     a   b     |     2   3         2   3       2   3          2   3      |                    |
#   |    / \ / \    |    / \           / \             / \            / \     |                    |
#   |   a   c   b   |   4   6         4   6           6   5          6   5    |    7:a             |
#   |  / \ /|\ / \  |  /|   |\           /|\         /|\            /|   |\   |    9:d             |
#   | a   d | e   b | 7 9   8 11        8 9 11      8 9 11         8 9  11 10 |                    |
#        /  |  \    |                                                         |    8:c             |
#   |   d   c   e   |                                                         |   11:e             |
#   |               |                                                         |   10:b             |
#   |               |  [0.0,0.2)      [0.2,0.5)    [0.5,0.6)      [0.6,1.0)   |                    |
#   |               |                                                         |                    |
def test_load_directly():
    population = main()
    mapping = population.renumber()
    my_ll_ts = _msprime.TreeSequence()
    my_ll_ts.load_records(list(population.records))
    ts = msprime.TreeSequence(my_ll_ts)
    return ts

def test_load_from_file():
    population = main()
    mapping = population.renumber()
    with tempfile.NamedTemporaryFile(mode='w') as f:
        population.write_records(f)
        f.flush()
        msprime.load_txt(f.name)

def test_records_match():
    ts = test_load_directly()
    my_ll_ts = ts.ll_tree_sequence
    for x, y in zip(list(population), ts.records()):
        # records match
        assert x == y
    return my_ll_ts

def test_records_can_iterate():
    ts = test_load_directly()
    for t in ts.tress():
        pass

def test_compare_trees():
    ts = test_load_directly()
    print("# this is what we get from tree iteration")
    for t, pyt in zip(ts.trees(), pytrees):
        print("#a tree")
        t  = t.parent_dict
        p  = parent_dict(pyt)
        print("#  msprime", t)
        print("#  trees.py", p)
        print("#   differ at")
        clashing_children = set(t.keys()).symmetric_difference(p.keys())
        for k in clashing_children:
            print("#    child:", k)
            try:
                v = t[k]
                print("#      <(msprime)<", v)
            except:
                print("#      <(msprime)<", "no entry")
            try:
                v = p[k]
                print("#      >(treespy)>", v)
            except:
                print("#      >(msprime)>", "no entry")

# this is what we get from tree iteration
#a tree
#  msprime {0: 5, 1: -1, 2: 7, 3: 5, 4: 7, 5: 9, 7: 9, 9: 10, 10: -1}
#  trees.py {0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 7: -1, 8: 10, 9: 10, 10: -1}
#   differ at
#    child: 6
#      <(msprime)< no entry
#      >(treespy)> 8
#    child: 8
#      <(msprime)< no entry
#      >(treespy)> 10
#a tree
#  msprime {0: 5, 1: -1, 2: 5, 3: 5, 4: -1, 5: 9, 9: 10, 10: -1}
#  trees.py {0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 7: -1, 8: 10, 9: 10, 10: -1}
#   differ at
#    child: 6
#      <(msprime)< no entry
#      >(treespy)> 8
#    child: 7
#      <(msprime)< no entry
#      >(treespy)> -1
#    child: 8
#      <(msprime)< no entry
#      >(treespy)> 10
#a tree
#  msprime {0: 5, 1: -1, 2: 5, 3: 5, 4: -1, 5: 8, 8: 10, 10: -1}
#  trees.py {0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 7: -1, 8: 10, 9: 10, 10: -1}
#   differ at
#    child: 6
#      <(msprime)< no entry
#      >(treespy)> 8
#    child: 7
#      <(msprime)< no entry
#      >(treespy)> -1
#    child: 9
#      <(msprime)< no entry
#      >(treespy)> 10
#a tree
#  msprime {0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 8: 10, 10: -1}
#  trees.py {0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 7: -1, 8: 10, 9: 10, 10: -1}
#   differ at
#    child: 7
#      <(msprime)< no entry
#      >(treespy)> -1
#    child: 9
#      <(msprime)< no entry
#      >(treespy)> 10
