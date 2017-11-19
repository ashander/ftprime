import msprime
from ftprime import ARGrecorder
from itertools import count
import random

from .breakpoints import (
    random_breakpoint,
    random_mutations,
)


"""
Same as wf but using __call__() interface..
"""

def wf_call(N, ngens, nsamples, survival=0.0, mutation_rate=0.0, simplify_interval=10,
       debug=False, seed=None) :
    '''
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population of size N
    for ngens generations, in which each individual survives with probability
    survival and only those who die are replaced.  The chromosome is 1.0
    Morgans long, and the mutation rate is in units of
    mutations/Morgan/generation.

    NOTE: mutation rate not implemented!

    Outputs an ARGrecorder object for the simulation.  In the final generation,
    a random set of individuals are chosen to be samples.
    '''
    if seed is not None:
        random.seed(seed)
    labels = count(0, 1)
    pop = [next(labels) for k in range(N)]
    # initial population
    init_ts = msprime.simulate(N, recombination_rate=1.0)
    init_samples = init_ts.samples()
    records = ARGrecorder(ts=init_ts,
                          node_ids={k:init_samples[k] for k in range(N)})

    for t in range(1, 1+ngens) :
        if debug:
            print("t:", t)
            print("pop:", pop)
            print(records)

        if (t % simplify_interval) == 0:
            records.simplify(pop)

        dead = [(random.random() > survival) for k in pop]
        # this is: offspring ID, lparent, rparent, breakpoint
        new_inds = [(next(labels), random.choice(pop), random.choice(pop),
                     random_breakpoint(), random_mutations(mutation_rate))
                    for k in range(sum(dead))]
        if debug:
            print("Replacing", sum(dead), "individuals.")
        js, parents, times, pops, offsprings, lefts, rights = \
            zip(*record_info_iter(new_inds, debug, time=t, dead=dead))
        for j, offspring in zip(js, offsprings):
            pop[j] = offspring
        records(parents=parents, times=times,
                populations=pops, childrens=((o, ) for o in offsprings),
                lefts=lefts, rights=rights)

    if debug:
        print("Done, now sampling.")
        print(records)
        print("pop:", pop)

    # restrict to a random subsample
    samples = random.sample(pop, nsamples)
    records.simplify(samples)

    if debug:
        print("Done.")
        print("samples:", samples)
        print(records)

    return records


def record_info_iter(new_inds, debug, time, dead):
    '''
    yields:
        pop_idx, parents, times, populations, childrens, lefts, rights
    '''
    j = 0
    for offspring, lparent, rparent, bp, muts in new_inds :
        while not dead[j] :
            j+=1
        if debug:
            print("--->", offspring, lparent, rparent, bp)
        if bp > 0.0 :
            yield j, lparent, time, 0, offspring, 0.0, bp
        if bp < 1.0 :
            yield j, rparent, time, 0, offspring, bp, 1.0

# TODO add something like above for muts
# for mut in muts:
#     records.add_mutation(position=mut, node=input_id,
#                          derived_state=1, ancestral_state=0)
