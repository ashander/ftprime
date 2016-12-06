from ftprime import merge_records
import msprime
from wf import wf

# test 1

new_rec=msprime.CoalescenceRecord(left=0.8, right=1.0, node=18, children=(22,), time=2, population=0)
old_recs=[
        msprime.CoalescenceRecord(left=0.0, right=0.6, node=18, children=(19,), time=2, population=0),
        msprime.CoalescenceRecord(left=0.6, right=1.0, node=18, children=(19, 20), time=2, population=0)]
merge_records(new_rec,old_recs)

old_recs

right_answer=[
    msprime.CoalescenceRecord(left=0.0, right=0.6, node=18, children=(19,), time=2, population=0),
    msprime.CoalescenceRecord(left=0.6, right=0.8, node=18, children=(19,20), time=2, population=0),
    msprime.CoalescenceRecord(left=0.8, right=1.0, node=18, children=(19, 20, 22), time=2, population=0),
    ]

if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
    raise ValueError('whoops')

# test 2

new_rec = msprime.CoalescenceRecord(left=0.5, right=1.0, node=4, children=(13,), time=6, population=0)
old_recs = [ msprime.CoalescenceRecord(left=0.0, right=0.9, node=4, children=(7, 12), time=6, population=0),
             msprime.CoalescenceRecord(left=0.9, right=1.0, node=4, children=(7,), time=6, population=0) ]
merge_records(new_rec,old_recs)

right_answer=[
    msprime.CoalescenceRecord(left=0.0, right=0.5, node=4, children=(7, 12), time=6, population=0),
    msprime.CoalescenceRecord(left=0.5, right=0.9, node=4, children=(7, 12, 13), time=6, population=0),
    msprime.CoalescenceRecord(left=0.9, right=1.0, node=4, children=(7, 13), time=6, population=0)
    ]

if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
    raise ValueError('whoops')

# test 3

new_rec = msprime.CoalescenceRecord(left=0.1, right=1.0, node=9, children=(11,), time=4, population=0)
old_recs = [msprime.CoalescenceRecord(left=0.0, right=0.1, node=9, children=(11,), time=4, population=0)]
merge_records(new_rec,old_recs)

right_answer=[
    msprime.CoalescenceRecord(left=0.0, right=0.1, node=9, children=(11,), time=4, population=0),
    msprime.CoalescenceRecord(left=0.1, right=1.0, node=9, children=(11,), time=4, population=0)
    ]

if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
    raise ValueError('whoops')


#### test simulation runs

records = wf(N=5,ngens=5,nsamples=5,survival=0.5)

for x in records.dump_records():
    print(x)

ts = records.tree_sequence()

for t in ts.trees():
    print(t)

print("Mean pairwise diversity:",ts.get_pairwise_diversity())
print("(should be zero)")

