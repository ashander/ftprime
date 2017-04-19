import msprime
from collections import OrderedDict


class ARGrecorder(OrderedDict):
    '''
    Keys are individual IDs, and values are tuples, whose first entry is a Node,
    which is a tuple of (birth time, population, is_sample),
    and the second entry is a ordered list of nonoverlapping Edgesets.
    This inherits from OrderedDict so that if individuals are added by birth
    order then records can easily be output ordered by time.

    Note that the 'time' fields must all be *strictly* greater than 0,
    and is in units of generations *ago*, i.e., 'reverse time'.
    '''

    def __init__(self):
        super(ARGrecorder, self).__init__()
        self.num_nodes = 0

    def __str__(self):
        ret = self.node_table().__str__()
        ret += "\n---------\n"
        ret += self.edgeset_table().__str__()
        return ret

    def add_individual(self, name, time, population=msprime.NULL_POPULATION,
                       is_sample=False):
        '''Add a new individual.
        We need to add individuals when they are *born*,
        rather than the first time they reproduce, to ensure
        that records are output in order by birth time of the parent.
        '''
        if name not in self:
            self[name] = (msprime.Node(time=time, population=population,
                                       name=name, is_sample=is_sample), [])
            self.num_nodes = max(self.num_nodes, 1+int(name))

    def add_record(self, left, right, parent, children):
        '''
        Add records corresponding to a reproduction event in which children (a
        tuple of IDs) inherit from parent (a single ID) on the interval
        [left,right).
        '''
        # unneeded but helpful for debugging
        if parent not in self.keys():
            raise ValueError("Parent " + str(parent) +
                             "'s birth time has not been recorded with " +
                             ".add_individual().")
        # time = self[parent][0]
        new_rec = msprime.Edgeset(
                parent=parent,
                children=children,
                left=left,
                right=right)
        merge_records(new_rec, self[parent][1])

    def edgesets(self):
        '''
        Returns an iterator of Edgesets, sorted in reverse order by
        the order the parents were input; so if parents are input in order
        by birth time, the output will be in the time order required by
        msprime.
        '''
        for y in reversed(self.values()):
            for z in y[1]:
                yield z

    def edgeset_table(self):
        '''
        Returns an EdgesetTable instance.
        '''
        table = msprime.EdgesetTable()
        for es in self.edgesets():
            table.add_row(
                left=es.left, right=es.right, parent=es.parent,
                children=es.children)
        return table

    def node_table(self):
        '''
        Returns a NodeTable instance, whose k-th row corresponds to node k.
        '''
        table = msprime.NodeTable()
        empty_node = msprime.Node(time=0.0)
        for k in range(self.num_nodes):
            try:
                node = self[k][0]
            except KeyError:
                node = empty_node
            table.add_row(flags=node.flags, time=node.time,
                          population=node.population)
        return table

    def dump_sample_table(self, out):
        '''
        Write out the table of info about the samples.
        '''
        out.write("id\tflags\tpopulation\ttime\n")
        for idx in self:
            node = self[idx][0]
            if node.is_sample():
                out.write("{}\t{}\t{}\t{}\n".format(idx, node.flags,
                                                    node.population, node.time))

    def tree_sequence(self, sites=None, mutations=None):
        '''
        Produce a tree sequence from the ARG.
        If sites and mutations is present (and are appropriate tables)
        then these will be added to the tree.
        '''
        if mutations is not None and sites is not None:
            ts = msprime.load_tables(
                nodes=self.node_table(),
                edgesets=self.edgeset_table(),
                sites=sites, mutations=mutations)
        else:
            ts = msprime.load_tables(
                    nodes=self.node_table(),
                    edgesets=self.edgeset_table())
        return ts

    def add_samples(self, samples, length, populations=None):
        '''
        Add phony records that stand in for sampling the IDs in `samples`,
        whose populations are given in `populations` (default: NULL), on a
        chromosome of total length `length`.
        '''
        if populations is None:
            populations = [msprime.NULL_POPULATION for x in samples]
        for k, parent in enumerate(samples):
            self.add_individual(k, time=0.0, population=populations[k],
                                is_sample=True)
            self.add_record(
                    left=0.0,
                    right=length,
                    parent=parent,
                    children=(k,))


def merge_records(new, existing):
    '''
    Incorporate a new record (l,r,x,c,t[x])
      into a list of existing ones (a,b,x,C,t[x]) sorted on left endpoint.
    Keeping them in sorted order simplifies the procedure
      (makes it so we don't have to split the new record).
    '''
    k = 0
    cur_left = new.left
    # print("MR: -----")
    # print("adding", new)
    # print("    to", existing)
    while (k < len(existing)) and (cur_left < new.right):
        left = existing[k].left
        right = existing[k].right
        parent = existing[k].parent
        children = existing[k].children
        # print("k:",k)
        # print("existing:",existing[k])
        # print("cur_left:",cur_left)
        if new.parent != parent:
            raise ValueError("Trying to merge records with different parents.")
        if right <= cur_left:
            # no overlap
            # print("no overlap")
            k += 1
            continue
        if cur_left < left:
            # print("dangling left")
            existing.insert(k, msprime.Edgeset(
                left=cur_left,
                right=min(new.right, left),
                parent=parent,
                children=new.children))
            cur_left = min(new.right, left)
            k += 1
            continue
        combined_children = tuple(sorted(children+new.children))
        combined_rec = msprime.Edgeset(
            left=cur_left,
            right=min(new.right, right),
            parent=new.parent,
            children=combined_children)
        if cur_left == left:
            # print("equal left")
            if new.right < right:
                # print("overlap right")
                mod_rec = msprime.Edgeset(
                        left=new.right,
                        right=right,
                        parent=parent,
                        children=children)
                existing[k] = combined_rec
                k += 1
                existing.insert(k, mod_rec)
                k += 1
            else:
                # print("dangling right")
                existing[k] = combined_rec
                k += 1
        else:
            # here we know that left < cur_left < right
            # print("overlap left")
            mod_rec = msprime.Edgeset(
                    left=left,
                    right=cur_left,
                    parent=parent,
                    children=children)
            existing[k] = mod_rec
            k += 1
            existing.insert(k, combined_rec)
            k += 1
            if new.right < right:
                # print("overlap right")
                existing.insert(k, msprime.Edgeset(
                    left=new.right,
                    right=right,
                    parent=parent,
                    children=children))
                k += 1
        cur_left = min(new.right, right)
    # add whatever's left at the end
    if cur_left < new.right:
        existing.insert(k, msprime.Edgeset(
            left=cur_left,
            right=new.right,
            parent=new.parent,
            children=new.children))
    # print("getting")
    # for x in existing:
    #     print("   ", x)
    return None
