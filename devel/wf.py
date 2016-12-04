from msprime import CoalescenceRecord
from itertools import count
from collections import OrderedDict
import random

def random_breakpoint() :
    return min(1.0,max(0.0, 2*random.random()-0.5)) 

def merge_records(new,existing) :
    #  Incorporate a new record (l,r,x,c,t[x])
    #    into a list of existing ones (a,b,x,C,t[x]) sorted on left endpoint.
    #  Keeping them in sorted order simplifies the procedure
    #    (makes it so we don't have to split the new record).
    k=0
    cur_left=new.left
    print("MR: -----")
    print("adding", new)
    print("    to", existing)
    while (k<len(existing)) and (cur_left<new.right):
        left,right,node,children,time,population=existing[k]
        print("k:",k)
        print("existing:",existing[k])
        print("cur_left:",cur_left)
        if new.node != node:
            raise ValueError("Trying to merge records with different parents.")
        if new.time != time:
            raise ValueError("Trying to merge records with different times.")
        if right<=cur_left:
            # no overlap
            print("no overlap")
            k+=1
            continue
        if cur_left < left:
            print("dangling left")
            existing.insert( k, CoalescenceRecord(
                left=cur_left, right=min(new.right,left), node=node, children=new.children, time=time, population=new.population) )
            cur_left=min(new.right,left)
            k+=1
            continue
        combined_children=tuple(sorted(children+new.children))
        combined_rec=CoalescenceRecord(
            left=cur_left, right=min(new.right,right), node=new.node, children=combined_children, time=new.time, population=new.population)
        if cur_left == left:
            print("equal left")
            if new.right < right:
                print("overlap right")
                mod_rec=CoalescenceRecord(
                        left=new.right, right=right, node=node, children=children, time=time, population=population )
                existing[k]=mod_rec
                k+=1
                existing.insert(k,combined_rec)
                k+=1
            else:
                print("dangling right")
                existing[k]=combined_rec
                k+=1
        else:
            # here we know that left < cur_left < right
            print("overlap left")
            mod_rec=CoalescenceRecord(
                    left=left, right=cur_left, node=node, children=children, time=time, population=population )
            existing[k]=mod_rec
            k+=1
            existing.insert(k,combined_rec)
            k+=1
            if new.right < right:
                print("overlap right")
                existing.insert(k,CoalescenceRecord(
                    left=new.right, right=right, node=node, children=children, time=time, population=new.population))
                k+=1
        cur_left=min(new.right,right)
    # add whatever's left at the end
    if cur_left < new.right:
        existing.insert(k,CoalescenceRecord(
            left=cur_left, right=new.right, node=new.node, children=new.children, time=new.time, population=new.population) )
    print("getting")
    for x in existing:
        print("   ", x)
    return None


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
    records = OrderedDict((x,[]) for x in pop)

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
                merge_records( CoalescenceRecord(left=0.0, right=bp, node=lparent, children=(offspring,), time=time[lparent], population=0), records[lparent] )
            if bp < 1.0 :
                merge_records( CoalescenceRecord(left=bp, right=1.0, node=rparent, children=(offspring,), time=time[rparent], population=0), records[rparent] )
        # print(records)

    # add phony records that stand in for sampling
    samples=random.sample(pop,nsamples)
    for k,parent in enumerate(samples):
        new_rec=CoalescenceRecord(left=0.0, right=1.0, node=parent, children=(k,), time=time[parent], population=0)
        if not parent in records:
            records[parent]=[new_rec]
        else:
            merge_records(new_rec,records[parent])

    return [ r for a in reversed(records.values()) for r in a ]
