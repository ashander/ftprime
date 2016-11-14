from ftprime.alg import main
from trees import trees

population = main()
#  what main() does
#    p = Population(0.0, state={'a': 1, '?': 0})
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
print("      the mapping:")
print(mapping)
print("\n state:")
print(population.state)
print("\n")
for t in trees(list(population)):
    print(t)
    pass

#  ------ feeding a renumbered version into the msprime alg T ------
#      the mapping:
#{0: 12, 1: 11, 2: 10, 3: 9, 4: 8, 5: 7, 6: 6, 7: 5, 8: 4, 9: 3, 10: 2, 11: 1}
#
# state:
#{'a': 5, 'c': 4, 'b': 2, 'e': 1, 'd': 3}
#
#
#	in: CoalescenceRecord(left=0.0, right=0.2, node=8, children=[3, 5], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.2, node=6, children=[1, 4], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=0.5, node=10, children=[6, 8], time=1.0, population=0)
#	in: CoalescenceRecord(left=0.0, right=1.0, node=11, children=[9, 10], time=2.0, population=0)
#([-1, 6, -1, 8, 6, 8, 10, -1, 10, 11, 11, -1], [[], [], [], [], [], [], [1, 4], [], [3, 5], [], [6, 8], [9, 10]])
#	out: CoalescenceRecord(left=0.0, right=0.2, node=8, children=[3, 5], time=0.0, population=0)
#	out: CoalescenceRecord(left=0.0, right=0.2, node=6, children=[1, 4], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.2, right=0.5, node=6, children=[1, 3, 4], time=0.0, population=0)
#([-1, 6, -1, 6, 6, -1, 10, -1, 10, 11, 11, -1], [[], [], [], [], [], [], [1, 3, 4], [], [], [], [6, 8], [9, 10]])
#	out: CoalescenceRecord(left=0.0, right=0.5, node=10, children=[6, 8], time=1.0, population=0)
#	out: CoalescenceRecord(left=0.2, right=0.5, node=6, children=[1, 3, 4], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=0.6, node=6, children=[1, 3, 4], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.5, right=1.0, node=9, children=[6, 7], time=1.0, population=0)
#([-1, 6, -1, 6, 6, -1, 9, 9, -1, 11, 11, -1], [[], [], [], [], [], [], [1, 3, 4], [], [], [6, 7], [], [9, 10]])
#	out: CoalescenceRecord(left=0.5, right=0.6, node=6, children=[1, 3, 4], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=7, children=[1, 2], time=0.0, population=0)
#	in: CoalescenceRecord(left=0.6, right=1.0, node=6, children=[3, 4], time=0.0, population=0)
#([-1, 7, 7, 6, 6, -1, 9, 9, -1, 11, 11, -1], [[], [], [], [], [], [], [3, 4], [1, 2], [], [6, 7], [], [9, 10]])
