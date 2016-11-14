from collections import namedtuple
from itertools import count
from msprime import CoalescenceRecord
from sortedcontainers import SortedSet, SortedList
from .merge import merge_records

Birth = namedtuple('Birth', 'left right x offspring')
Birth.__doc__ = ','.join(['left (Individual)',
                          'right (Individual)',
                          'x (float)',
                          'offspring (Individual)'])

class NodeItems(dict):
    '''
    a list of items keyed by parent node

    a dict but with a list of nodes as the value
    '''

    def __setitem__(self, key, value):
        '''
        keep list of values for each key, add to it with each set
        '''
        self.setdefault(key, []).append(value)

    def update(self, *E, **F):
        raise NotImplementedError

    def get(self, k):
        raise NotImplementedError


def main():
    p = Population(0.0, state={'a': 1, '?': 0})
    births = [Birth('a', '?', 1.0, 'b'), ]
    deaths = ['?',]
    p._onegen(births, deaths)

    births = [Birth('a', 'b', 0.5, 'c'), ]
    p._onegen(births, deaths)

    births = [Birth('a', 'c', 0.2, 'd'),
              Birth('c', 'b', 0.6, 'e'), ]
    p._onegen(births, deaths)

    print('--------------------- fin ------------')
    for r in p: print(r)

    return p


class Population(object):

    def __init__(self, time=0.0, state=dict()):
        '''A pop that spits out records

            state (dict): mapping individuals to labels
            time (float): the present time
        '''
        self._time = time
        self._state = state
        self.left_end = 0.0
        self.right_end = 1.0
        self._records = []
        self._breakpoints = SortedSet()
        # next_label (callable): to generate labels (integers)
        next_label = 1 + max(v for v in state.values())
        self._labels = count(next_label)
        self.add_breakpoint(0.0)
        self.add_breakpoint(1.0)
        self._final = False

    @property
    def state(self):
        return self._state

    @property
    def time(self):
        return self._time

    @property
    def next_label(self):
        return next(self._labels)

    @property
    def breakpoints(self):
        return self._breakpoints

    def add_breakpoint(self, breakpoint):
        self._breakpoints.add(breakpoint)

    def __iter__(self):
        return iter(self._records)

    def _onegen(self, births, deaths=None):
        '''
        args:
            births (list): of Birth(left, right, x, j) left, right, j are individuals
            deaths: (list) of individuals
        '''
        old_labels = dict()  # map individuals to their old label
        records = NodeItems()
        print("Births", births)
        self._time += 1.0
        for i, j, x, k in births:
            if x == self.right_end:
                # i only
                old_i, new_i, old_labels = self._update_parent_label(i, old_labels)
            elif x == self.left_end:
                # j only
                old_j, new_j, old_labels = self._update_parent_label(j, old_labels)
                k_lab = self.next_label
            else:
                # both
                old_i, new_i, old_labels = self._update_parent_label(i, old_labels)
                old_j, new_j, old_labels = self._update_parent_label(j, old_labels)
            k_lab = self.next_label
            self.state[k] = k_lab
            print(self.state)
            self.add_breakpoint(x)
            bp = self._breakpoints
            x_idx = bp.bisect_left(x)
            for i in range(1, len(bp)):
                if i <= x_idx:
                    records[old_i] = CoalescenceRecord(bp[i-1], bp[i], old_i,
                                                       sorted([new_i, k_lab]),
                                                       self.time, population=0)
                elif i > x_idx:
                    records[old_j] = CoalescenceRecord(bp[i-1], bp[i], old_j,
                                                       sorted([new_j, k_lab]),
                                                       self.time, population=0)
        if deaths is not None:
            for ind in deaths:
                try:
                    self.state.pop(ind)
                except:
                    pass
        for recs in records.values():
            comp,  inc = merge_records(recs, debug=False)
            assert len([i for i in inc]) == 0
            self._records.extend(comp)

    def _update_parent_label(self, ind, old_labels):
        try:
            old_lab = old_labels[ind]
            new_lab = self.state[ind]
        except KeyError:
            old_lab = self.state[ind]
            old_labels[ind] = old_lab
            new_lab = self.next_label
            self.state[ind] = new_lab
        return old_lab, new_lab, old_labels

    def renumber(self):
        self._final = True
        maxn = self._maxnode()
        maxt = self.time
        self._records = self._finalize(self, maxt=maxt, maxn=maxn)
        for ind, label in self.state.items():
            self.state[ind] = maxn - label

    def _maxnode(self):
        ''' maximum node number

        such that renumbering using this results in nodes with numbers > 0'''
        if not self._final:
            raise ValueError("should not be called outside renumbering")
        return self.next_label

    def _finalize(self, iterable, maxt, maxn):
        recs = SortedList([], key=lambda rec: rec.time)
        for r in iterable:
            time = maxt - r.time
            node = maxn - r.node
            renumbered_children = (maxn - c for c in r.children)
            recs.add(CoalescenceRecord(
                    time=time,
                    node=node,
                    children=sorted(renumbered_children),
                    population=r.population,
                    left=r.left,
                    right=r.right,
                ))
        return recs

if __name__ == "__main__":
    main()
