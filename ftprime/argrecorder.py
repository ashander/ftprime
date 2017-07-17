import msprime
import numpy as np
from _msprime import NODE_IS_SAMPLE

NULL_ID = -1

class ARGrecorder(object):
    '''
    To record the ARG, this keeps track of
        - a NodeTable and an EdgesetTable
        - information about what has happened since the last simplification
        - a mapping from (input) individual IDs to (output) node IDs
    The reason for maintaining this mapping, rather than having (input)
    individual IDs always equal to the (output) node IDs is to allow periodic
    simplification, which decouples the two.

    Between simplification steps, the EdgesetTable is appended to, so it is
    always "up to date".  However, the NodeTable is *not* kept up to date,
    because its `time` fields are recorded in *time ago*; we also keep track of
        - a list of birth times of individual IDs
    which are translated to time-ago at each simplification step, and appended to the Node Table.

    The internal state is stored using
        - self.node_ids[k] : gives the output Node ID corresponding to the
          input individual ID k+self.input_min
    and self.birth_times and self.populations stores the corresponding birth
    times (in forwards-time) and populations.  These need to contain information
    about every input individual who might be seen in a new edgeset (ie, who
    might be a parent).  

    Assumes that input individual IDs are *nondecreasing*.
    '''

    def __init__(self, max_indivs=1000):
        super(ARGrecorder, self).__init__()
        # maximum number of individuals that can be kept track of in the internal state
        self.n = max_indivs
        # this is output node ID of the next new node
        self.num_nodes = 0
        self.max_time = 0.0
        # last time simplify() was run
        self.last_update_time = 0.0
        self.ts = msprime.TreeSequence
        # edges since the last simplify():
        #  need to keep these separate only because they get appended in reverse order
        self.new_edgesets = msprime.EdgesetTable()
        # minimum currently active input ID
        self.input_min = 0
        # list of corresponding output node IDs
        self.node_ids = np.empty(self.n, dtype='int32')
        self.node_ids.fill(NULL_ID)
        # and populations
        self.populations = np.empty(self.n, dtype='int32')
        self.populations.fill(NULL_ID)
        # list of corresponding (forwards-time) birth times
        self.birth_times = np.empty(self.n, dtype=float)
        self.birth_times.fill(np.nan)

    def __str__(self):
        ret = "Input ID offset:\n"
        ret += self.input_min.__str__()
        ret += "Node IDs:\n"
        ret += self.node_ids.__str__()
        ret += "Populations:\n"
        ret += self.populations.__str__()
        ret += "Birth times:\n"
        ret += self.birth_times.__str__()
        ret += "NodeTable:\n"
        ret += self.node_table().__str__()
        ret += "\n---------\n"
        ret += self.edgeset_table().__str__()
        return ret

    def ix(self, x):
        """
        Translate input IDs into indices in the internal arrays.
        """
        return x + self.input_min

    def ox(self, x):
        """
        Translate indices in the internal arrays to input IDs.
        """
        return x - self.input_min

    def add_individual(self, name, time, population=msprime.NULL_POPULATION):
        '''
        Add a new individual, recording its birth time and assigning it a new output Node ID.
        '''
        i = self.ix(name)
        if self.node_ids[i] == NULL_ID:
            self.node_ids[i] = self.num_nodes
            self.num_nodes += 1
            self.birth_times[i] = time
            self.populations[i] = population
            self.max_time = max(self.max_time, time)
        else:
            raise ValueError("Attempted to add " + str(name) + 
                             ", who already exits, as a new individual.")

    def add_record(self, left, right, parent, children):
        '''
        Add records corresponding to a reproduction event in which children (a
        tuple of IDs) inherit from parent (a single ID) on the interval
        [left,right).
        '''
        # unneeded but helpful for debugging
        i = self.ix(parent)
        if self.node_ids[i] == NULL_ID:
            raise ValueError("Parent " + str(parent) +
                             "'s birth time has not been recorded with " +
                             ".add_individual().")
        # time = self.birth_times[parent]
        out_parent = self.node_ids[self.ix(parent)]
        out_children = tuple([self.node_ids[self.ix(u)] for u in children])
        self.new_edgesets.add_row(parent=out_parent,
                                  children=out_children,
                                  left=left,
                                  right=right)

    def update_nodes(self, nodes):
        """
        Updates the NodeTable to include information from `birth_times` and `populations`.
        """
        # number of previous nodes
        m = nodes.num_rows
        new_times = np.empty(self.num_nodes, dtype=float)
        new_populations = np.empty(self.num_nodes, dtype='int32')
        # new flags will all be 0
        new_flags = np.zeros(self.num_nodes, dtype='uint32')
        new_times[:m] = nodes.time + self.max_time - self.max_update_time
        new_times[m:] = np.nan
        new_populations[:m] = nodes.populations
        new_populations[m:] = msprime.NULL_POPULATION
        new_flags[:m] = nodes.flags
        # which ones need to be added
        add_these = (self.node_ids != NULL_ID) & (self.node_ids >= m)
        # indices in new array of these
        ii = self.node_ids[add_these]
        new_times[ii] = self.birth_times[add_these]
        new_populations[ii] = self.populations[add_these]
        nodes.set_columns(time=new_times,
                          flags=new_flags,
                          population=new_populations)

    def update_edgesets(self, edgesets):
        """
        Updates the EdgesetTable to include edges stored in new_edgesets: these
        just need to be reversed and appended.
        """
        m = self.new_edgesets.num_rows
        n = edgesets.num_rows
        mc = edgesets.children_length)
        nc = sum(self.new_edgesets.children_length)
        new_left = np.empty(m+n, dtype=float)
        new_right = np.empty(m+n, dtype=float)
        new_parent = np.empty(m+n, dtype='int32')
        new_children = np.empty(mc+nc, dtype='int32')
        new_children_length = np.empty(m+n, dtype='uint32')
        new_left[:m] = edgesets.left
        new_left[m+n:m:-1] = self.new_edgesets.left
        new_right[:m] = edgesets.right
        new_right[m+n:m:-1] = self.new_edgesets.right
        new_parent[:m] = edgesets.parent
        new_parent[m+n:m:-1] = self.new_edgesets.parent
        new_children_length[:m] = edgesets.children_length
        new_children_length[m+n:m:-1] = self.new_edgesets.children_length
        new_children[:mc] = edgesets.children
        # can't do this as above because of requirement that children are sorted
        k = 0
        for j in self.new_edgesets.children_length:
            new_children[mc+nc-k:mc+nc-(k+j+1):-1] = self.new_edgesets.children[k:k+j]
            k += j
        edgesets.set_columns(left=new_left,
                          right=new_right,
                          parent=new_parent,
                          children=new_children,
                          children_length=new_children_length)

    def simplify(self, samples):
        """
        `samples` should be a list of all "currently living" input individual
        IDs: i.e., anyone who might be a parent or a sample in the future.
        Updates the current internal state by
            1. adding pending nodes and edges to the tree sequence
            2. running .simplify() on the tables to remove information
                irrelevant to anyone not in `samples` 
            3. updating the current state
        """
        # update the TreeSequence
        nodes = msprime.NodeTable()
        edgesets = msprime.EdgesetTable()
        self.ts.dump_tables(nodes=nodes, edgesets=edgesets)
        self.update_nodes(nodes)
        self.update_edgesets(edgesets)
        self.ts = msprime.load_tables(nodes=nodes, edgesets=edgesets)
        self.ts.simplify(samples=self.node_ids[self.ix(samples)])
        # update the internal state
        old_populations = self.populations[self.ix(samples)]
        old_birth_times = self.birth_times[self.ix(samples)]
        self.num_nodes = self.ts.num_nodes
        self.last_update_time = self.max_time
        self.new_edgesets.reset()
        self.input_min = min(samples)
        # new indices
        ii = self.ix(samples) # note this differs from use above as input_min has changed
        self.node_ids.fill(NULL_ID)
        self.node_ids[ii] = self.ts.samples()
        self.populations.fill(NULL_ID)
        self.populations[ii] = old_populations
        self.birth_times.fill(np.nan)
        self.birth_times[ii] = old_birth_times
        ## run some checks
        self.ts.dump_tables(nodes=nodes)
        assert all(self.populations[ii] == nodes.population[self.ts.samples()])
        assert all(self.max_time - self.birth_times[ii] == nodes.time[self.ts.samples()])


    def node_table(self):
        '''
        Returns a NodeTable instance, whose k-th row corresponds to node k.
        Translates time to 'time ago' by subtracting time from max_time, which
        by default is the largest time seen so far.
        '''
        nodes = msprime.NodeTable()
        self.ts.dump_tables(nodes=nodes)
        return(nodes)

    def edgeset_table(self):
        '''
        Returns a NodeTable instance, whose k-th row corresponds to node k.
        Translates time to 'time ago' by subtracting time from max_time, which
        by default is the largest time seen so far.
        '''
        edgesets = msprime.EdgesetTable()
        self.ts.dump_tables(edgesets=edgesets)
        return(edgesets)

    def dump_sample_table(self, out):
        '''
        Write out the table of info about the samples.
        '''
        out.write("id\tflags\tpopulation\ttime\n")
        for idx in self.ts.samples():
            node = ts.node(idx)
            out.write("{}\t{}\t{}\t{}\n".format(idx, 
                                                node.flags,
                                                node.population, 
                                                node.time))

    def mark_samples(self, samples):
        """
        Mark these individuals as samples internally (but do not simplify).
        """
        tables = self.ts.dump_tables()
        nodes = tables.edges
        new_flags = nodes.flags
        new_flags &= ~NODE_IS_SAMPLE
        new_flags[samples] |= NODE_IS_SAMPLE
        nodes.fill_columns(time=nodes.time,
                           population=nodes.population,
                           flags=new_flags)
        self.ts.load_tables(**tables)

    def tree_sequence(self):
        return self.ts
