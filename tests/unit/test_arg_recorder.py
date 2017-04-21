import pytest
import six
import msprime

from ftprime import merge_records, ARGrecorder


def num_lines(string):
    return len(str(string).rstrip().split('\n'))


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


@pytest.mark.parametrize(('new_rec', 'old_recs', 'exception'), [
    (msprime.Edgeset(left=0.1, right=1.0, parent=9, children=(11,)),
     [
         msprime.Edgeset(left=0.0, right=0.1, parent=10, children=(11,)),
     ],
     ValueError
     ),
])
def test_merge_bad_records(new_rec, old_recs, exception):
    with pytest.raises(exception):
        merge_records(new_rec, old_recs)


class TestARGrecorderEmpty():
    def setup_method(self):
        '''setup for all tests in class'''
        self.arg = ARGrecorder()

    def test_init(self):
        assert self.arg.num_nodes == 0

    def test_str(self):
        assert num_lines(self.arg) == 3

    def test_add_individual(self):
        with pytest.raises(ValueError):
            self.arg.add_individual('a', time=0)  # name is converted to an int
        self.arg.add_individual('0', time=0.0)  # so some strings work...
        assert self.arg.num_nodes == 1
        self.arg.add_individual(0, time=1.0)  # won't add with same name
        assert self.arg.num_nodes == 1
        self.arg.add_individual(100, time=1.0)  # num nodes is 1 + int(name)
        assert self.arg.num_nodes == 101

    def test_add_record(self):
        with pytest.raises(ValueError):
            self.arg.add_record(left=0.0, right=1.0, parent=2, children=(0, 1))

    def test_edgeset(self):
        with pytest.raises(StopIteration):
            next(self.arg.edgesets())  # empty to start


class TestARGrecorderFilled():
    def setup_method(self):
        '''setup for all tests in class'''
        self.parent1 = 2
        self.time1 = 1.0
        self.parent2 = 5
        self.time2 = 2.0
        self.arg = ARGrecorder()
        self.arg.add_individual(self.parent1, time=self.time1)
        self.arg.add_record(left=0.0, right=1.0, parent=2, children=(0, 1))
        self.arg.add_individual(self.parent2, time=self.time2)

    def test_add_record(self):
        assert len(self.arg) == 2
        # adding more children to the one parent does not increase length of
        # the ARGrecorder
        self.arg.add_record(left=0.0, right=1.0, parent=2, children=(3, 4))
        assert len(self.arg) == 2
        # adding an new parent does increase the length
        self.arg.add_individual(5, time=10.0)
        self.arg.add_record(left=0.0, right=1.0, parent=5, children=(3, 4))
        assert len(self.arg) == 2

    def test_edgeset(self):
        self.arg.add_individual(5, time=0.0)
        with pytest.raises(StopIteration):
            e_iter = self.arg.edgesets()
            next(e_iter)
            next(e_iter)  # no edgesets for added parent

        # add a new parent and see the edgesets come out in reverse order
        self.arg.add_record(left=0.0, right=1.0, parent=5, children=(3, 4))
        e_iter = self.arg.edgesets()
        e1 = next(e_iter)
        e2 = next(e_iter)
        assert e1.parent == 5
        assert e2.parent == self.parent1

    def test_edgeset_table(self):
        self.arg.add_individual(5, time=0.0)
        self.arg.add_record(left=0.0, right=1.0, parent=5, children=(3, 4))
        et = self.arg.edgeset_table()
        assert num_lines(et) == 3

    def test_node_table(self):
        nt = self.arg.node_table()
        assert num_lines(nt) == self.arg.num_nodes + 1
        # one line for each node and one for header

    def test_sample_table(self):
        out = six.StringIO()
        self.arg.dump_sample_table(out)
        print(out.getvalue())
        assert num_lines(out.getvalue()) == 1  # no samples


class TestARGrecorderGoodSampledPopulation():
    def setup_method(self):
        '''setup for all tests in class'''
        self.arg = ARGrecorder()
        self.arg.add_individual(2, time=3.0)
        self.arg.add_individual(3, time=3.0)
        self.arg.add_individual(4, time=1.0)
        self.arg.add_record(left=0.5, right=1.0, parent=3, children=(4, ))
        self.arg.add_record(left=0.0, right=0.5, parent=2, children=(4, ))
        self.samples = [3, 4]
        self.arg.add_samples(samples=self.samples, length=1.0)
        # trees in format of msprime.SparseTree.parent_dict:
        # {node: parent}
        # after sampling 3, 4 nodes 0, 1 have parents 3, 4
        # over the whole range; 0.5 is a breakpoint for internal node
        # switching from 2 to 3
        self.tree_1 = {0: 3, 1: 4, 4: 2}
        self.tree_2 = {0: 3, 1: 4, 4: 3}

    def test_add_samples(self):
        out = six.StringIO()
        self.arg.dump_sample_table(out)
        print('\nSelf:\n')
        print(self.arg)
        print('\nSamples:\n')
        print(out.getvalue())
        assert num_lines(out.getvalue()) == len(self.samples) + 1  # for header

    def test_tree_sequence(self):
        ts = self.arg.tree_sequence()
        trees = [t.parent_dict for t in ts.trees()]
        assert self.tree_2 in trees
        assert self.tree_1 in trees
