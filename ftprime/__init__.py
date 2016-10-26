from numpy.random import randint, random, poisson
from collections import namedtuple
from itertools import count


Chromosome = namedtuple('Chromosome',
                        'identity breakpoint lparent rparent left_end right_end')
default_chromosome_values = (-1, -2, 0.0, 1.0)
Chromosome.__new__.__defaults__ = default_chromosome_values
Chromosome.__doc__ = \
    '''a chromosome is produced by recombination between two parents

        Args:
            identity: integer ID for this chromosome
            breakpoint: separates the right from left parts
            lparent: integer ID for the parent of the left part
            rparent: integer ID for the parent of the left part

        Optional args:
            left_end
            right_end -- NOTE no error checking to make sure bp is between these
    '''


def chromosome_str(self):
    if self.breakpoint is not None:
        bp = "{breakpoint:.{bpdigits}f})(".format(breakpoint=self.breakpoint,
                                                 bpdigits=2)
        info = "{lp:<4}~ {breakpoint} ~{rp:>4}".format(lp=self.lparent,
                                                       rp=self.rparent,
                                                       breakpoint=bp,)
    elif self.breakpoint is None:
        if self.rparent != self.lparent:
            raise ValueError("R parent and L parent must be equal if None break")
        info = "     {}       ".format(self.rparent)
    else:
        raise TypeError("Illegal type of breakpoint")

    output = "{id:8}: {le:.{epdigits}f}.[" + info + "]{re:.{epdigits}f}."
    return output.format(
        id=self.identity,
        le=self.left_end,
        re=self.right_end,
        epdigits=0,
    )

Chromosome.__str__ = chromosome_str
Chromosome.__repr__ = Chromosome.__str__


class Individual(object):
    def __init__(self, chromosomes: tuple, age: int=0):
        '''chromosome ids only are stored'''
        try:
            try:
                it = iter(chromosomes)
            except Exception as e:
                raise TypeError(e.value.message)
            self._cs = tuple([next(it).identity for i in range(2)])
        except Exception:
            raise ValueError("Can't initialize with with non-Chromosomes" +
                             " -- must have identity set!")

        self._sc = list(self._cs)
        self._sc.reverse()
        self._a = age

    def age(self) -> int:
        return self._a

    def chromosomes(self) -> tuple:
        '''the chromosomes

        return:
            chromosomes: tuple of Chromosomes
        '''
        return self._cs

    def __str__(self):
        return "Individual(" + \
            ", ".join((str(c) for c in self.chromosomes())) + \
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
        parents_and_breakpoints = ((bp, self._cs[ind], self._sc[ind])
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

        self._total_size = self._size = int(size)
        self._counter = count()
        self._chromosomes = dict()
        initial_chroms = ((i, Chromosome(i, bp))
                          for i, bp in zip(self._counter, random(size=size * 2)))
        self._chromosomes.update(initial_chroms)
        self._individuals = [Individual(chromosomes=(c1, c2)) for c1, c2 in
                             zip(*[self.chromosomes()] * 2)]  # takes 2

    def chromosomes(self):
        return iter(self._chromosomes.values())

    def individuals(self):
        return iter(self._individuals)

    def size(self):
        return self._size

    def _mating_pairs(self, average_offspring):
        offspring_it = iter(poisson(average_offspring, self.size()))
        for i1, i2, on in zip(*[self.individuals()] * 2, offspring_it):
            yield i1, i2, on

    def _gametes_from_pairs(self, average_offspring):
        for i1, i2, on in self._mating_pairs(average_offspring):
            yield zip(i1.gametes(on, tagger=self._counter),
                      i2.gametes(on, tagger=self._counter))

    def generation(self, average_offspring=1.5):
        for childs in self._gametes_from_pairs(average_offspring):
            for cc in childs:
                yield Individual(cc)

if __name__ == "__main__":
    size = 5
    p = Population(size=size)
    p.generation(average_offspring=2.5)
