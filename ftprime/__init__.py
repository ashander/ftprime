from .alg import Population, main
from .chromosome import Chromosome, Segment
from numpy.random import randint, random
from itertools import count


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

    @property
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
