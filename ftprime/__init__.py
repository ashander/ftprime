from .chromosome import Chromosome, Segment
from numpy.random import randint, random, poisson, permutation
from itertools import count
from .merge import merge_records, SortedList, NodeItems, CoalescenceRecord


class Individual(object):
    def __init__(self, chromosomes: tuple, age: int=0):
        '''chromosome ids only are stored'''
        try:
            try:
                it = iter(chromosomes)
            except Exception as e:
                raise TypeError(e.value.message)
            self._cs = tuple([next(it) for i in range(2)])
            self._cs[0].identity
            self._cs[1].identity
        except Exception:
            raise ValueError("Can't initialize with with non-Chromosomes" +
                             " -- must have identity set!")

        self._sc = list(self._cs)
        self._sc.reverse()
        self._a = age

    def age(self) -> int:
        return self._a

    def __iter__(self) -> tuple:
        '''the chromosomes

        return:
            chromosomes: tuple of Chromosomes
        '''
        return iter(self._cs)

    def __str__(self):
        return "Individual(" + \
            ",".join((str(c) for c in iter(self))) + \
            ")"

    def __repr__(self):
        return str(self)

#    @profile
    def raw_gametes(self, number: int):
        ''' iterator over gametes

        args:
            number: integer number of gametes to produce

        return:
            bp (float)
            lp_id (int)
            rp_id (int)
        '''
        parents_and_breakpoints = ((bp,
                                    self._cs[ind].identity,
                                    self._sc[ind].identity,
                                    )
                                   for ind in randint(2, size=number)
                                   for bp in random(size=number))
        for r in range(number):
            yield next(parents_and_breakpoints)

    def gametes(self, number: int, tagger=count()):
        for other_args, idnum in zip(self.raw_gametes(number), tagger):
            yield Chromosome(idnum, *other_args)


class Population(object):
    def __init__(self, size: int=10, verbose=False):
        '''make  some folks

        size (int): how many? can be passed a float but rounds down
        '''
        if verbose:
            if type(size) is not int:
                raise Warning("size", size, "coerced to integer")
        self._time = 0.0
        self._total_size = self._size = int(size)
        self._chromosome_ids = count()
        self._final = False
        # decreasing 'infinite' count(__import__('sys').maxsize, -1)
        self._records = SortedList([], key=lambda rec: rec.time)
        self._nc = NodeItems()
        initial_chroms = (Chromosome(id)
                          for ct, id in zip(range(size * 2),
                                            self._chromosome_ids))
        initial_inds = (Individual(chromosomes=(c1, c2)) for c1, c2 in
                        zip(*[initial_chroms] * 2))   # takes 2
        self._individuals = [i for i in initial_inds]  # self._store_inds()

    def __str__(self):
        return "size: {} & time: {}\n".format(self.size(), self.time()) +\
            "\n".join((str(r) for r in self._records))

    def __repr__(self):
        return str(self)

    def time(self):
        return self._time

    def _updatetime(self, increment=0.5):
        self._time += increment

    def __iter__(self):
        return iter(self._individuals)

    def chromosomes(self):
        return (c for i in iter(self) for c in i)

    def individuals_randomized(self):
        return iter(permutation(self._individuals))

    def size(self):
        return self._size

    def _maxtime(self):
        ''' the last time that is stored in a record'''
        if not self._final:
            raise ValueError("should not be called outside finalize")
        last_rec = self._records.pop()
        self._records.append(last_rec)
        return min(last_rec.time, self._time)

    def _maxnode(self):
        if not self._final:
            raise ValueError("should not be called outside finalize")
        return next(self._chromosome_ids) - 1

    def generation(self, average_offspring: float):
        ''' make one generation and store the resulting segments,

        then merge these to complete records information'''
        new_inds = self._store_inds(self._next_generation(average_offspring))
        self._updatetime()
        self._individuals = new_inds  # non-overlapping gens
        self._total_size = self.size() + len(self._individuals)
        self._size = len(self._individuals)
        self._merge_segs_to_records()

    def _store_inds(self, iterable_of_individuals):
        new_inds = []
        for ind in iterable_of_individuals:
            new_inds.append(ind)
            for c in ind:
                if c.breakpoint is None:
                    s = Segment(left=c.left_end,
                                right=c.right_end,
                                node=c.lparent,
                                children=(c.identity,),
                                time=self.time())
                    self._nc[s.node] = s
                else:
                    sl = Segment(left=c.left_end,
                                 right=c.breakpoint,
                                 node=c.lparent,
                                 children=(c.identity,),
                                 time=self.time())
                    sr = Segment(left=c.breakpoint,
                                 right=c.right_end,
                                 node=c.rparent,
                                 children=(c.identity,),
                                 time=self.time())
                    self._nc[sl.node] = sl
                    self._nc[sr.node] = sr
        return new_inds

    def records(self):
        return iter(self._records)

    def unmerged_records(self):
        yield from self._nc.items()

    def finalize(self):
        self._final = True
        maxn = self._maxnode()
        maxt = self._maxtime()
        self._records = self._finalize(self.records(), maxt=maxt, maxn=maxn)
#        for k, v in self.unmerged_records():
#            for r in self._finalize(v, maxt=maxt, maxn=maxn):
#                self._nc[k] = r

    def _finalize(self, iterable, maxt, maxn):
        recs = SortedList([], key=lambda rec: rec.time)
        for r in iterable:
            time = maxt - r.time
            node = maxn - r.node
            recs.add(CoalescenceRecord(
                    time=time,
                    node=node,
                    children=[maxn - c for c in r.children],
                    population=r.population,
                    left=r.left,
                    right=r.right,
                ))
        return recs

    def _merge_singletons(self):
        unmerged = (s for s in self._unmerged_children(segments))
        children = (c for s in segments for c in s.children)
        for c in children:
            try:
                for s in self._nc.pop(c):
                    yield s
            except:
                pass

    def _merge_segs_to_records(self, debug=False):
        for k in self._nc.keys():
            segments = self._nc.pop(k)
            comp, incomp = merge_records(segments)
            for c in comp:
                self._records.add(c)

            for i in incomp:
                self._nc[k] = i

    # helpers to advance generation
    def _next_generation(self, average_offspring):
        for childs in self._gametes_from_pairs(average_offspring):
            for cc in childs:
                yield Individual(cc)

    def _gametes_from_pairs(self, average_offspring):
        for i1, i2, on in self._mating_pairs(average_offspring):
            yield zip(i1.gametes(on, tagger=self._chromosome_ids),
                      i2.gametes(on, tagger=self._chromosome_ids))

    def _mating_pairs(self, average_offspring):
        offspring_it = iter(poisson(average_offspring, self.size()))
        for i1, i2, on in zip(*[self.individuals_randomized()] * 2, offspring_it):
            yield i1, i2, on

    def write_records(self, output, header=True, precision=6):
        """
        Writes the records for this tree sequence to the specified file in a
        tab-separated format. If ``header`` is True, the first line of this
        file contains the names of the columns, i.e., ``left``, ``right``,
        ``node``, ``children``, ``time`` and ``population``. After the
        optional header, the records are written to the file in
        tab-separated form in order of non-decreasing time. The ``left``,
        ``right`` and ``time`` fields are base 10 floating point values
        printed to the specified ``precision``. The ``node`` and
        ``population`` fields are base 10 integers. The ``children`` column
        is a comma-separated list of base 10 integers, which must contain at
        least two values.

        Example usage:

        >>> with open("records.txt", "w") as records_file:
        >>>     tree_sequence.write_records(records_file)

        :param File output: The file-like object to write the tab separated
        output.
        :param bool header: If True, write a header describing the column
        names in the output.
        :param int precision: The number of decimal places to print out for
        floating point columns.

        License:
            this function is Copyright (C) 2015 Jerome Kelleher and reused
            under the GPLv3
        """
        if not self._final:
            raise ValueError("need to finalize records first")
        if header:
            print(
                "left", "right", "node", "children",
                "time", "population", sep="\t", file=output)

        for record in self.records():
            children = ",".join(str(c) for c in record.children)
            row = (
                "{left:.{precision}f}\t"
                "{right:.{precision}f}\t"
                "{node}\t"
                "{children}\t"
                "{time:.{precision}f}\t"
                "{population}\t").format(
                    precision=precision,
                    left=record.left, right=record.right,
                    node=record.node, children=children,
                    time=record.time, population=record.population)
            print(row, file=output)

if __name__ == "__main__":
    size = 5
    p = Population(size=size)
    for g in p.generation(average_offspring=2.5):
        print(g)

