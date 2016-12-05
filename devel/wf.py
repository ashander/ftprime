from pedrecorder import PedigreeRecorder
from itertools import count
import random

def random_breakpoint() :
    return min(1.0,max(0.0, 2*random.random()-0.5)) 

def wf(N,ngens,nsamples,survival=0.0) :
    '''
    SIMPLE simulation of a bisexual, haploid Wright-Fisher population
    of size N for ngens generations,
    in which each individual survives with probability survival
    and only those who die are replaced.

    Keeps track of a dict whose keys are individuals
    and whose entries are CoalescenceRecords.

    Outputs a list of CoalescenceRecords.
       CoalescenceRecord(left, right, node, children, time-from-end, population)
    in nondecreasing order by time-from-end;
    in the final generation, a random set of nsamples individuals are labeled 1...nsamples.

    '''
    labels = count(nsamples,1)
    pop = [ next(labels) for k in range(N) ]
    time = dict( (x,ngens+1) for x in pop )
    # insert individuals into records in order of birth to ensure we output records in time-sorted order
    records = PedigreeRecorder((x,[]) for x in pop)

    for t in range(ngens) :
        print("t:",t)
        print("pop:", pop)
        dead = [ (random.random() > survival) for k in pop ]

        # this is: offspring ID, lparent, rparent, breakpoint
        new_inds = [ (next(labels),random.choice(pop),random.choice(pop),random_breakpoint()) for k in range(sum(dead)) ]
        j=0
        for offspring,lparent,rparent,bp in new_inds :
            print(offspring,lparent,rparent,bp)
            while not dead[j] :
                j+=1
            pop[j]=offspring
            j+=1
            time[offspring]=ngens-t
            records[offspring]=[]
            if bp > 0.0 :
                records.add_record(left=0.0, right=bp, parent=lparent, children=(offspring,), time=time[lparent], population=0)
            if bp < 1.0 :
                records.add_record( left=bp, right=1.0, parent=rparent, children=(offspring,), time=time[rparent], population=0)
        # print(records)

    # add phony records that stand in for sampling
    samples=random.sample(pop,nsamples)
    for k,parent in enumerate(samples):
        records.add_record(left=0.0, right=1.0, parent=parent, children=(k,), time=time[parent], population=0)

    return [ r for a in reversed(records.values()) for r in a ]
