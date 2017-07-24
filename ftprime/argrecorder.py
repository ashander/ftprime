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

    Must be initialized with a tree sequence ``ts``, this will serve as
    the history of this first generation of individuals.  The ``samples`` of
    this tree sequence correspond to the individuals in the first generation.

    If different than the sample IDs of ``ts``, ``first_generation`` has the
    node IDs of the first generation.
    '''

    def __init__(self, ts, first_generation=None, time=0.0, max_indivs=1000):
        """
        :param ts TreeSequence: The tree sequence defining relationships between
            the first genreation of individuals.
        :param list first_generation: The input IDs of the first generation that
            should be mapped onto the samples in the tree sequence ``ts``.
        :param float time: The (forwards) time of the "present" at the start of
            the simulation.
        :param int max_indivs: The largest number of individuals (past or present)
            that we will track at any given time.
        """
        if first_generation is None:
            first_generation = ts.samples()
        # maximum number of individuals that can be kept track of in the internal state
        self.n_max = max_indivs
        # this is output node ID of the next new node
        self.num_nodes = ts.num_nodes # n
        # this is the largest (forwards) time seen so far
        self.max_time = time  # T
        # last (forwards) time we updated node times
        self.last_update_time = time  # T_0
        # number of nodes that have the time right
        self.last_update_node = ts.num_nodes
        assert len(ts.samples()) == len(first_generation)
        # minimum currently active input ID
        self.input_min = 0  # m_0
        # list of output node IDs indexed by input labels
        self.node_ids = np.empty(self.n_max, dtype='int32') # L
        self.reset_node_ids(first_generation, ts.samples())
        # the actual tables that get updated
        #  DON'T actually store ts, just the tables:
        self.nodes = msprime.NodeTable()
        self.edgesets = msprime.EdgesetTable()
        self.sites = msprime.SiteTable()
        self.mutations = msprime.MutationTable()
        ts.dump_tables(nodes=self.nodes, edgesets=self.edgesets,
                       sites=self.sites, mutations=self.mutations)
        # list of site positions, maintained as site tables don't have
        #   efficient checking for membership
        self.site_positions = self.sites.positions

    def __str__(self):
        ret = "Input ID offset m_0:\n"
        ret += str(self.input_min)
        ret += "Max time so far:\n"
        ret += str(self.max_time)
        ret += "Last update time:"
        ret += str(self.last_update_time)
        ret += "Node IDs:\n"
        ret += self.node_ids.__str__()
        ret += "Tables:\n"
        for x in (self.nodes, self.edgesets, self.sites, self.mutations):
            ret += x.__str__()
        ret += "\n---------\n"
        return ret

    def reset_node_ids(self, generation, samples):
        """
        Reset the internal mapping of input to node IDs except for each
        individual in generation is mapped to the corresponding one in samples.
        """
        self.node_ids.fill(NULL_ID)
        for j, x in zip(generation, samples):
            self.node_ids[self.input_min + j] = x

    def add_individual(self, input_id, time, population=msprime.NULL_POPULATION):
        '''
        Add a new individual, recording its birth time and assigning it a new
        output Node ID.
        '''
        i = input_id - self.input_min
        if self.node_ids[i] == NULL_ID:
            self.nodes.add_row(flags=0, population=population, time=time)
            self.node_ids[i] = self.num_nodes
            self.num_nodes += 1
            self.max_time = max(self.max_time, time)
        else:
            raise ValueError("Attempted to add " + str(input_id) +
                             ", who already exits, as a new individual.")

    def add_record(self, left, right, parent, children):
        '''
        Add records corresponding to a reproduction event in which children (a
        tuple of IDs) inherit from parent (a single ID) on the interval
        [left,right).
        '''
        # unneeded but helpful for debugging
        out_parent = self.node_ids[parent - self.input_min]
        if out_parent == NULL_ID:
            raise ValueError("Parent " + str(parent) +
                             "'s birth time has not been recorded with " +
                             ".add_individual().")
        out_children = tuple([self.node_ids[u - self.input_min] for u in children])
        self.edgesets.add_row(parent=out_parent,
                              children=out_children,
                              left=left,
                              right=right)

    def add_mutation(self, position, node, derived_state):
        """
        Adds a new mutation to mutation table, and a new site if necessary as well.
        """
        if position not in self.site_positions:
            site = self.sites.num_rows
            self.sites.add_row(position=position, ancestral_state=ancestral_state)
            self.site_positions.append(position)
        else:
            site = self.site_positions.index(position)
        self.mutations.add_row(site=site, node=node, derived_state=derived_state)

    def update_times(self):
        """
        Update times in the NodeTable.
        """
        dt = self.max_time - self.last_update_time
        self.nodes.time[:self.last_update_node] =
            self.nodes.time[:self.last_update_node] + dt
        self.nodes.time[self.last_update_node:] =
            self.max_time - self.nodes.time[self.last_update_node:]
        self.last_update_time = self.max_time
        self.last_update_node = self.nodes.num_rows

    def simplify(self, samples):
        """
        `samples` should be a list of all "currently living" input individual
        IDs: i.e., anyone who might be a parent or a sample in the future.
        """
        # update times
        self.update_times()
        ts = msprime.load_tables(nodes=self.nodes, edgesets=self.edgesets,
                                 sites=self.sites, mutations=self.mutations)
        sample_nodes = [self.node_ids[x - self.input_min] for x in samples]
        ts.simplify(samples=sample_nodes)
        # update the internal state
        self.num_nodes = ts.num_nodes
        self.last_update_node = ts.num_nodes
        self.input_min = min(samples)
        # update index map
        self.reset_node_ids(samples, ts.samples()):
        ts.dump_tables(nodes=self.nodes, edgesets=self.edgesets,
                       sites=self.sites, mutations=self.mutations)

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
        new_flags = self.nodes.flags
        new_flags &= ~NODE_IS_SAMPLE
        new_flags[samples] |= NODE_IS_SAMPLE
        self.nodes.fill_columns(time=self.nodes.time,
                                population=self.nodes.population,
                                flags=new_flags)
        self.ts.load_tables(**tables)

    def tree_sequence(self, samples):
        """
        Return the tree sequence for a given set of input samples.
        """
        self.update_times()
        # problem: this may not work if simplify hasn't just been run
        ts = msprime.load_tables(nodes=self.nodes, edgesets=self.edgesets,
                                 sites=self.sites, mutations=self.mutations)
        sample_nodes = [self.node_ids[x - self.input_min] for x in samples]
        ts.simplify(samples=sample_nodes)
        return ts
