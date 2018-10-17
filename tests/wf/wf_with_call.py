import msprime
from ftprime import ARGrecorder
from itertools import count
import numpy as np

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

    Outputs an ARGrecorder object for the simulation.  In the final generation,
    a random set of individuals are chosen to be samples.
    '''
    if seed is not None:
        np.random.seed(seed)
    else:
        seed = np.random.randint(1e8)
    labels = count(0, 1)
    pop = [next(labels) for k in range(N)]
    # initial population
    init_ts = msprime.simulate(N, recombination_rate=1.0, random_seed=seed)
    init_samples = init_ts.samples()
    records = ARGrecorder(ts=init_ts, node_ids={k:init_samples[k] for k in range(N)})

    for t in range(1, 1+ngens) :
        if debug:
            print("t:", t)
            print("pop:", pop)
            print(records)

        if (t % simplify_interval) == 0:
            records.simplify(pop)

        dead = [(np.random.uniform() > survival) for k in pop]
        # this is: offspring ID, lparent, rparent, breakpoint
        new_inds = [(next(labels), np.random.choice(pop), np.random.choice(pop), 
                     random_breakpoint(), random_mutations(mutation_rate))
                    for k in range(sum(dead))]
        j=0
        if debug:
            print("Replacing", sum(dead), "individuals.")
        for offspring, lparent, rparent, bp, muts in new_inds :
            if debug:
                print("--->", offspring, lparent, rparent, bp)
            while not dead[j] :
                j+=1
            pop[j]=offspring
            j+=1
            if bp > 0.0 :
                records(parent=lparent, time=t, population=0, child=offspring,
                        left=0.0, right=bp)
            if bp < 1.0 :
                records(parent=rparent, time=t, population=0, child=offspring,
                        left=bp, right=1.0)
            for mut in muts:
                records.add_mutation(position=mut, node=offspring, 
                                     derived_state=b'1', ancestral_state=b'0')

    if debug:
        print("Done, now sampling.")
        print(records)
        print("pop:", pop)

    # restrict to a random subsample
    samples = np.random.choice(pop, nsamples, replace=False)
    records.simplify(samples)

    if debug:
        print("Done.")
        print("samples:", samples)
        print(records)

    return records

