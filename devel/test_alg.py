from ftprime.alg import main
import msprime
import _msprime
from trees import trees

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
#{'?': 0, 'a': 2, 'b': 3}
#Births [Birth(left='a', right='b', x=0.5, offspring='c')]
#{'a': 4, 'c': 6, 'b': 5}
#Births [Birth(left='a', right='c', x=0.2, offspring='d'), Birth(left='c', right='b', x=0.6, offspring='e')]
#{'a': 7, 'c': 8, 'd': 9, 'b': 5}
#{'e': 11, 'b': 10, 'a': 7, 'c': 8, 'd': 9}
#--------------------- fin ------------
#CoalescenceRecord(left=0.0, right=1.0, node=1, children=(2, 3), time=1.0, population=0)
#CoalescenceRecord(left=0.0, right=0.5, node=2, children=(4, 6), time=2.0, population=0)
#CoalescenceRecord(left=0.5, right=1.0, node=3, children=(5, 6), time=2.0, population=0)
#CoalescenceRecord(left=0.0, right=0.2, node=4, children=(7, 9), time=3.0, population=0)
#CoalescenceRecord(left=0.6, right=1.0, node=5, children=(10, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.0, right=0.2, node=6, children=(8, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.2, right=0.5, node=6, children=(8, 9, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.5, right=0.6, node=6, children=(8, 9, 11), time=3.0, population=0)
#CoalescenceRecord(left=0.6, right=1.0, node=6, children=(8, 9), time=3.0, population=0)


mapping = population.renumber()
print("\n  ------ feeding a renumbered version into the msprime alg T ------")
print("      the mapping from NEW to OLD:")
print(mapping)
print("\n state:")
print(population.state)
print("(in old numbers:", [(k, mapping[v]) for k, v in population.state.items()])
print("\n")
for t in trees(list(population)):
    print("  trees: ", t)
    t_1_renum = [[mapping[j] for j in els] if len(els) > 0 else [] for els in t[1]]
    t_1_renum.reverse()
    print("(in old numbers:", t_1_renum)
    pass
#
#  ------ feeding a renumbered version into the msprime alg T ------
#      the mapping from NEW to OLD:
#{0: 10, 1: 9, 2: 8, 3: 7, 4: 6, 5: 5, 6: 4, 7: 3, 8: 2, 9: 1}
#
# state:
#{'a': 4, 'b': 1, 'e': 0, 'd': 2, 'c': 3}
#(in old numbers: [('a', 6), ('b', 9), ('e', 10), ('d', 8), ('c', 7)]
#
#
#	in: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(2, 4), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.2, node=5, children=(0, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.5, node=9, children=(5, 7), time=1.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=1.0, node=10, children=(8, 9), time=2.0, population=0)
#  trees:  ([5, -1, 7, 5, 7, 9, -1, 9, 10, 10, -1], [[], [], [], [], [], (0, 3), [], (2, 4), [], (5, 7), (8, 9)])
#(in old numbers: [[2, 1], [5, 3], [], [8, 6], [], [10, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.0, right=0.2, node=7, children=(2, 4), time=0.0, population=0)
#	out: CoalescenceRecord(left=0.0, right=0.2, node=5, children=(0, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.2, right=0.5, node=5, children=(0, 2, 3), time=0.0, population=0)
#  trees:  ([5, -1, 5, 5, -1, 9, -1, 9, 10, 10, -1], [[], [], [], [], [], (0, 2, 3), [], [], [], (5, 7), (8, 9)])
#(in old numbers: [[2, 1], [5, 3], [], [], [], [10, 8, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.0, right=0.5, node=9, children=(5, 7), time=1.0, population=0)
#	out: CoalescenceRecord(left=0.2, right=0.5, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=0.6, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=1.0, node=8, children=(5, 6), time=1.0, population=0)
#  trees:  ([5, -1, 5, 5, -1, 8, 8, -1, 10, 10, -1], [[], [], [], [], [], (0, 2, 3), [], [], (5, 6), [], (8, 9)])
#(in old numbers: [[2, 1], [], [5, 4], [], [], [10, 8, 7], [], [], [], [], []]
#	out: CoalescenceRecord(left=0.5, right=0.6, node=5, children=(0, 2, 3), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=6, children=(0, 1), time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=5, children=(2, 3), time=0.0, population=0)
#  trees:  ([6, 6, 5, 5, -1, 8, 8, -1, 10, 10, -1], [[], [], [], [], [], (2, 3), (0, 1), [], (5, 6), [], (8, 9)])
#(in old numbers: [[2, 1], [], [5, 4], [], [10, 9], [8, 7], [], [], [], [], []]
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
with open('tmp.tsv', 'w') as f:
    population.write_records(f)

ts = msprime.load_txt('tmp.tsv')
my_ll_ts = ts.ll_tree_sequence

for t in ts.trees():
    print(t)
# this is what we get from tree iteration
#{0: 5, 1: -1, 2: 7, 3: 5, 4: 7, 5: 9, 7: 9, 9: 10, 10: -1}
#{0: 5, 1: -1, 2: 5, 3: 5, 4: -1, 5: 9, 9: 10, 10: -1}
#{0: 5, 1: -1, 2: 5, 3: 5, 4: -1, 5: 8, 8: 10, 10: -1}
#{0: 6, 1: 6, 2: 5, 3: 5, 4: -1, 5: 8, 6: 8, 8: 10, 10: -1}
#
