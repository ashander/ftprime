import msprime
import _msprime
from trees import trees
from wf import wf, merge_records

def build_tree_sequence(records, mutations=[]):
    ts = _msprime.TreeSequence()
    ts.load_records(records)
    ts.set_mutations(mutations)
    return msprime.TreeSequence(ts)

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


records = wf(N=5,ngens=5,nsamples=2,survival=0.5)

for x in records:
    print(x)

ts = build_tree_sequence(records)

for t in ts.trees():
    print(t)
