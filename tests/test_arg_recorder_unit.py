from ftprime import merge_records
import msprime


def test_1():

    new_rec=msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(22,))
    old_recs=[
            msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19,)),
            msprime.Edgeset(left=0.6, right=1.0, parent=18, children=(19, 20))]
    merge_records(new_rec,old_recs)

    old_recs

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19,)),
        msprime.Edgeset(left=0.6, right=0.8, parent=18, children=(19,20)),
        msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(19, 20, 22)),
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')

def test_2():

    new_rec = msprime.Edgeset(left=0.5, right=1.0, parent=4, children=(13,))
    old_recs = [ msprime.Edgeset(left=0.0, right=0.9, parent=4, children=(7, 12)),
                 msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7,)) ]
    merge_records(new_rec,old_recs)

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.5, parent=4, children=(7, 12)),
        msprime.Edgeset(left=0.5, right=0.9, parent=4, children=(7, 12, 13)),
        msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7, 13))
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')

def test_3():

    new_rec = msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,))
    old_recs = [msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,))]
    merge_records(new_rec,old_recs)

    right_answer=[
        msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,)),
        msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,))
        ]

    if not all([ a==b for a,b in zip(right_answer,old_recs) ]):
        raise ValueError('whoops')
