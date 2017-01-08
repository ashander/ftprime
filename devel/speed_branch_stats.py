
import msprime
from timeit import default_timer as timer

ts = msprime.simulate(40, recombination_rate=1000)

A = [ [a] for a in ts.samples() ]
n = len(A)
def f(x):
    return [ float(x[j]*(1-x[k]) + (1-x[j])*x[k]) for j in range(n) for k in range(n) if j<k ] 

start = timer()

brlen = [ msprime.branch_length_diversity(ts,A[j],A[k]) for j in range(len(A)) for k in range(len(A)) if j<k ]

mid1 = timer()

bsv = msprime.branch_stats_vector(ts,A,f)

mid2 = timer()

ni = msprime.branch_stats_vector_node_iter(ts,A,f,method='length')

end = timer()

print("Should be zero:",max([abs(a-b) for a,b in zip(brlen,bsv)]))
print("Should be zero:",max([abs(a-b) for a,b in zip(brlen,ni)]))

print("branch_length_diversity():", mid1-start)
print("branch_stats_vector()", mid2-mid1)
print("branch_stats_vector_node_iter()", end-mid2)
