'''
combine partial information into full records

Each diploid individual contains two chromosomes, `i` and `j`, inherited from
its two parents:
    `i: [0--i_p1--x_i\--i_p2--1]`
    `j: [0--j_p2--x_j\--j_p2--1]`
Each meiosis takes the two chromosomes of the diploid individual, labels one `r`
And one `l`, and draws a breakpoint x, producing an offspring `k: [0--l--x_k\--r--1]`
that combines the two parent chromosomes. This event yields partial inforamtion on
potential coalescence events (using `msprime` record tuple l, r, n, (c), p, t
(and ignoring p and t). WLOG for several gametes produced from one individual,
we have

- pick x_k1, x_k2, ...x_kn. together with the breakpoints of the parents, we
have, in order along the line:
- [0, sorted(x_k1, x_k2, ..., x_kn, x_l, x_r, 1.0]

'''
from sortedcontainers import SortedSet
from msprime.trees import CoalescenceRecord


# class CoalescenceRecord(CoalescenceRecord):
#     # over ride coalescence records so it has a more readable str representation
#     __slots__ = ()
#     legend = "t    : [L,  R) - node <-(children)\n"
#
#     def __str__(self):
#         s = "{} : [{}, {}) - {} <-{}".format(
#             round(self.time, 3),
#             round(self.left, 3), round(self.right, 2),
#             self.node, tuple(self.children),
#         )
#
#         return s
#

ChildrenSet = lambda : SortedSet(key=lambda x: x)

def merge_records(l: list, debug: bool=False):
    '''
    merge records to provide a list of complete, and of incomplete records

    args:
        l: a list of CoalesenceRecords of length n
        debug (False)

    return:
        complete: an iterator over complete records
        incomplete: an iterator over incomplete records
    '''
    records = []
    pop = l[0].population
    node = l[0].node
    children =   ChildrenSet() # use set for easy removal and addition
    time = []
    left = None
    endpoints, childs, times, add, inds = _prepare_records_to_merge(l)
    for k, ind in enumerate(inds):  # inds aren't used
        if add[k]:
            if left is not None and children != ChildrenSet():
                assert time != []
                if left != endpoints[k]:
                    records.append(CoalescenceRecord(time=max(time), left=left,
                                                    right=endpoints[k], node=node,
                                                    children=tuple(children), population=pop))
            elif left is not None and children == ChildrenSet():
                assert time == []
                break

            # add children, time, start a new record
            for c in childs[k]:
                children.add(c)
            time.append(times[k])
            left = endpoints[k]
        elif not add[k]:
            # close off the old record
            assert children != ChildrenSet(), [str(ll) for ll in l]
            assert time != []
            if left != endpoints[k]:
                # CRAPPY FIX!
                # if we don't re-add kids here then if moving from (8,9,11)
                # to 8,9 we'll remove the 8 with current example
                for c in childs[k]:
                    children.add(c)
                records.append(CoalescenceRecord(time=max(time), left=left,
                                                 right=endpoints[k], node=node,
                                                 children=tuple(children), population=pop))
            # save children but remove the current one
            for c in childs[k]:
                children.remove(c)
            time.remove(times[k])  # is there a better way to do this?
            # start a new record
            left = endpoints[k]
        else:
            raise ValueError("elements of `add` should be True or False")

    comp = _complete_records(records)
    incomp = _incomplete_records(records)

    if debug:  # debugging
        print("\nWith inputs...")
        for ll in l:
            print(ll)
        complete = [c for c in comp]
        incomplete = [c for c in incomp]
        assert len(complete) + len(incomplete) == len(records), \
            "comp: {}\n incomp:{}\n all:{}\n".format(complete,
                                                     incomplete,
                                                     records)
        print("Merged")
        for c in complete:
            print(c)
        print("... leaving")
        for i in incomplete:
            print(i)
        print("... (all records created")
        for r in records:
            print(r)
        print("...)")
        return complete, incomplete

    return comp, incomp


def fix_incomplete(record, added_child=None):
    children = sorted(record.children + [added_child, ])
    return CoalescenceRecord(time=record.time,
                             left=record.left,
                             right=record.right,
                             node=record.node,
                             children=children,
                             population=record.population,
                             )


def _complete_records(records):
    '''complete (children > 1)'''
    for record in records:
        if len(record.children) > 1:
            yield record


def _incomplete_records(records):
    '''incomplete (children == 1)'''
    for record in records:
        if len(record.children) == 1:
            yield record

def _prepare_records_to_merge(l: list):
    '''
    args:
        l: a list of objects of length n
           they must have properties: left, time

    returns a zip object containing three length 2n lists:
        endpoints: a sorted list consisting of all left and rigth endpoints in l
        childs: list of children from the records
        times: list of times
        add: list of True/False indicating whether to add children (i.e.,
             whether the endpoint is left
        ind: list of indices in 1...n, which each appearing twice so that so
             that if add[k] is True then l[ind[k]].left = endpoints[k] and if
             add[k] is False then l[ind[k]].right = endpoints[k]. For convenience
             both the children and time elements of these nodes are provided also above
    '''
    endpoint_index = 0 # where, in the tuples below, does endpoint_appear?
    endpoint_tuples = ( ((rec.left, rec.children, rec.time, True, ind),
                         (rec.right, rec.children, rec.time, False, ind)) for
                       ind, rec in enumerate(l) )
    master_list = ( item for sublist in endpoint_tuples for item in sublist )
    # zip up (ie construct a list of first element, a list of secoond, etc,
    # the tuples after sorting them by the endpoints
    return zip(*sorted(master_list, key=lambda x: x[endpoint_index]))
