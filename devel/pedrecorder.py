import msprime
import _msprime
from collections import OrderedDict

class PedigreeRecorder(OrderedDict):
    '''
    Keys are individuals whose values are ordered lists of nonoverlapping CoalescenceRecords.
    This inherits from OrderedDict so that if individuals are added by birth order
    then records can easily be output ordered by time.
    '''

    def add_individual(self,name):
        self[name]=[]

    def add_record(self,left,right,parent,children,time,population):
        '''
        Add records corresponding to a reproduction event
        in which children inherits from parent on the interval [left,right). 
        '''
        if not parent in self:
            self[parent]=[]
        new_rec=msprime.CoalescenceRecord(left=left,right=right,node=parent,children=children,time=time,population=population)
        merge_records(new_rec,self[parent])

    def dump_records(self):
        return flatten_once(reversed(self.values()))

    def tree_sequence(self,mutations=None):
        ts = _msprime.TreeSequence()
        ts.load_records(list(self.dump_records()))
        if mutations is not None:
            ts.set_mutations(mutations)
        return msprime.TreeSequence(ts)

    def add_samples(self,samples,times,populations):
        # add phony records that stand in for sampling the IDs in samples
        # with corresponding times
        for k,parent in enumerate(samples):
            self.add_record(left=0.0, right=1.0, parent=parent, children=(k,), time=times[k], population=populations[k])


def flatten_once(x):
    '''
    Returns an iterator over the elements of a once-nested list.
    '''
    for y in x:
        for z in y:
            yield z


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
        if right<=cur_left:
            # no overlap
            # print("no overlap")
            k+=1
            continue
        if cur_left < left:
            # print("dangling left")
            existing.insert( k, msprime.CoalescenceRecord(
                left=cur_left, right=min(new.right,left), node=node, children=new.children, time=time, population=new.population) )
            cur_left=min(new.right,left)
            k+=1
            continue
        combined_children=tuple(sorted(children+new.children))
        combined_rec=msprime.CoalescenceRecord(
            left=cur_left, right=min(new.right,right), node=new.node, children=combined_children, time=new.time, population=new.population)
        if cur_left == left:
            # print("equal left")
            if new.right < right:
                # print("overlap right")
                mod_rec=msprime.CoalescenceRecord(
                        left=new.right, right=right, node=node, children=children, time=time, population=population )
                existing[k]=mod_rec
                k+=1
                existing.insert(k,combined_rec)
                k+=1
            else:
                # print("dangling right")
                existing[k]=combined_rec
                k+=1
        else:
            # here we know that left < cur_left < right
            # print("overlap left")
            mod_rec=msprime.CoalescenceRecord(
                    left=left, right=cur_left, node=node, children=children, time=time, population=population )
            existing[k]=mod_rec
            k+=1
            existing.insert(k,combined_rec)
            k+=1
            if new.right < right:
                # print("overlap right")
                existing.insert(k,msprime.CoalescenceRecord(
                    left=new.right, right=right, node=node, children=children, time=time, population=new.population))
                k+=1
        cur_left=min(new.right,right)
    # add whatever's left at the end
    if cur_left < new.right:
        existing.insert(k,msprime.CoalescenceRecord(
            left=cur_left, right=new.right, node=new.node, children=new.children, time=new.time, population=new.population) )
    # print("getting")
    # for x in existing:
    #     print("   ", x)
    return None

