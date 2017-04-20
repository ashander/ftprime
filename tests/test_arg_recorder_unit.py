import pytest
import msprime

from ftprime import merge_records


@pytest.mark.parametrize(('new_rec', 'old_recs', 'right_answer'), [
    (msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(22, )),
     [
         msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19, )),
         msprime.Edgeset(left=0.6, right=1.0, parent=18, children=(19, 20))],
     [
         msprime.Edgeset(left=0.0, right=0.6, parent=18, children=(19, )),
         msprime.Edgeset(left=0.6, right=0.8, parent=18, children=(19, 20)),
         msprime.Edgeset(left=0.8, right=1.0, parent=18, children=(19, 20, 22)),
     ]),
    (msprime.Edgeset(left=0.5, right=1.0, parent=4, children=(13, )),
     [
         msprime.Edgeset(left=0.0, right=0.9, parent=4, children=(7, 12)),
         msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7, ))
      ],
     [
         msprime.Edgeset(left=0.0, right=0.5, parent=4, children=(7, 12)),
         msprime.Edgeset(left=0.5, right=0.9, parent=4, children=(7, 12, 13)),
         msprime.Edgeset(left=0.9, right=1.0, parent=4, children=(7, 13))
     ]),
    (msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,)),
     [
         msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,)),
     ],
     [
         msprime.Edgeset(left=0.0, right=0.1, parent=9, children=(11,)),
         msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,))
     ]),
])
def test_merge_records(new_rec, old_recs, right_answer):
    merge_records(new_rec, old_recs)
    if not all([a == b for a, b in zip(right_answer, old_recs)]):
        raise ValueError('whoops')
