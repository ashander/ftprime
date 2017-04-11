from ftprime import ARGrecorder
from itertools import count
import random

def random_breakpoint() :
    return min(1.0,max(0.0, 2*random.random()-0.5))

def wf(N, ngens, nsamples, survival=0.0, debug=False) :
    '''
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population
    of size N for ngens generations,
    in which each individual survives with probability survival
    and only those who die are replaced.

    Keeps track of a dict whose keys are individuals
    and whose entries are Edgesets.

    Outputs an ARGrecorder object for the simulation.
    In the final generation, a random set of nsamples individuals are labeled 1...nsamples.

    '''
    labels = count(nsamples,1)
    pop = [ next(labels) for k in range(N) ]
    # insert individuals into records in order of birth to ensure we output records in time-sorted order
    records = ARGrecorder()
    for x in pop:
        records.add_individual(x,ngens+1)

    for t in range(ngens) :
        if debug:
            print("t:",t)
            print("pop:", pop)
        dead = [ (random.random() > survival) for k in pop ]

        # this is: offspring ID, lparent, rparent, breakpoint
        new_inds = [ 
                (next(labels),random.choice(pop),random.choice(pop),random_breakpoint()) 
                for k in range(sum(dead)) ]
        j=0
        for offspring,lparent,rparent,bp in new_inds :
            if debug:
                print(offspring,lparent,rparent,bp)
            while not dead[j] :
                j+=1
            pop[j]=offspring
            j+=1
            records.add_individual(offspring,ngens-t)
            if bp > 0.0 :
                records.add_record(left=0.0, right=bp, parent=lparent, children=(offspring,))
            if bp < 1.0 :
                records.add_record( left=bp, right=1.0, parent=rparent, children=(offspring,))
        # print(records)

    # add phony records that stand in for sampling
    samples=random.sample(pop,nsamples)
    records.add_samples(samples=samples,length=1.0)

    return records
