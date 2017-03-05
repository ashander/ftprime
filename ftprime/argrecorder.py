import msprime
import _msprime
from collections import OrderedDict

class ARGrecorder(OrderedDict):
    '''
    Keys are individuals, and values are tuples, whose first entry is the birth time of the parent,
    and the second entry is a ordered list of nonoverlapping CoalescenceRecords.
    This inherits from OrderedDict so that if individuals are added by birth order
    then records can easily be output ordered by time.
    
    Note that the 'time' fields must all be *strictly* greater than 0.
    '''

    def add_individual(self,name,time):
        '''Add a new individual.
        We need to add individuals when they are *born*,
        rather than the first time they reproduce, to ensure
        that records are output in order by birth time of the parent.
        '''
        if not name in self:
            self[name]=(time,[])

    def add_record(self,left,right,parent,children,population=0):
        '''
        Add records corresponding to a reproduction event in which children (a
        tuple of IDs) inherit from parent (a single ID) on the interval
        [left,right). 
        '''
        if not parent in self.keys():
            raise ValueError("Parent "+str(parent)+"'s birth time has not been recorded with .add_individual().")
        time=self[parent][0]
        new_rec=msprime.CoalescenceRecord(
                left=left,
                right=right,
                node=parent,
                children=children,
                time=time,
                population=population)
        merge_records(new_rec,self[parent][1])

    def dump_records(self):
        '''
        Returns an iterator of CoalescenceRecords, sorted in reverse order by
        the order the parents were input; so if parents are input in order
        by birth time, the output will be in the time order required by
        msprime.
        '''
        for y in reversed(self.values()):
            for z in y[1]:
                yield z

    def tree_sequence(self,samples,mutations=None):
        '''
        Produce a tree sequence from the ARG.
        '''
        ts = _msprime.TreeSequence()
        ts.load_records(
                records=list(self.dump_records()))
        if mutations is not None:
            ts.set_mutations(mutations)
        return msprime.TreeSequence(ts)

    def add_samples(self,samples,length,populations=None):
        '''
        Add phony records that stand in for sampling the IDs in `samples`,
        whose populations are given in `populations`, on a chromosome of
        total length `length`.
        '''
        if populations is None:
            populations = [0 for x in samples]
        for k,parent in enumerate(samples):
            self.add_record(left=0.0, right=length, parent=parent, 
                        children=(k,), population=populations[k])


def merge_records(new,existing) :
    '''
    Incorporate a new record (l,r,x,c,t[x])
      into a list of existing ones (a,b,x,C,t[x]) sorted on left endpoint.
    Keeping them in sorted order simplifies the procedure
      (makes it so we don't have to split the new record).
    '''
    k=0
    cur_left=new.left
    # print("MR: -----")
    # print("adding", new)
    # print("    to", existing)
    while (k<len(existing)) and (cur_left<new.right):
        left,right,node,children,time,population=existing[k]
        # print("k:",k)
        # print("existing:",existing[k])
        # print("cur_left:",cur_left)
        if new.node != node:
            raise ValueError("Trying to merge records with different parents.")
        if new.time != time:
            raise ValueError("Trying to merge records with different times.")
        if right <= cur_left:
            # no overlap
            # print("no overlap")
            k+=1
            continue
        if cur_left < left:
            # print("dangling left")
            existing.insert(k, msprime.CoalescenceRecord(
                left=cur_left, 
                right=min(new.right,left), 
                node=node, 
                children=new.children, 
                time=time, 
                population=new.population))
            cur_left=min(new.right,left)
            k+=1
            continue
        combined_children=tuple(sorted(children+new.children))
        combined_rec=msprime.CoalescenceRecord(
            left=cur_left, 
            right=min(new.right,right), 
            node=new.node, 
            children=combined_children, 
            time=new.time, 
            population=new.population)
        if cur_left == left:
            # print("equal left")
            if new.right < right:
                # print("overlap right")
                mod_rec=msprime.CoalescenceRecord(
                        left=new.right, 
                        right=right, 
                        node=node, 
                        children=children, 
                        time=time, 
                        population=population)
                existing[k]=combined_rec
                k+=1
                existing.insert(k,mod_rec)
                k+=1
            else:
                # print("dangling right")
                existing[k]=combined_rec
                k+=1
        else:
            # here we know that left < cur_left < right
            # print("overlap left")
            mod_rec=msprime.CoalescenceRecord(
                    left=left, 
                    right=cur_left, 
                    node=node, 
                    children=children, 
                    time=time, 
                    population=population)
            existing[k]=mod_rec
            k+=1
            existing.insert(k,combined_rec)
            k+=1
            if new.right < right:
                # print("overlap right")
                existing.insert(k,msprime.CoalescenceRecord(
                    left=new.right, 
                    right=right, 
                    node=node, 
                    children=children, 
                    time=time, 
                    population=new.population))
                k+=1
        cur_left=min(new.right,right)
    # add whatever's left at the end
    if cur_left < new.right:
        existing.insert(k,msprime.CoalescenceRecord(
            left=cur_left, 
            right=new.right, 
            node=new.node, 
            children=new.children, 
            time=new.time, 
            population=new.population))
    # print("getting")
    # for x in existing:
    #     print("   ", x)
    return None

