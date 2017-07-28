import msprime
import numpy as np

NULL_ID = -1

def null_tree_sequence():
    return msprime.load_tables(nodes=msprime.NodeTable(),
                               edgesets=msprime.EdgesetTable())

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
    which are translated to time-ago at each simplification step, and appended
    to the Node Table.

    The internal state is stored using
        - ``self.node_ids[k]`` : the output Node ID corresponding to the input
          individual ID ``k``.

    Assumes that input individual IDs are *nondecreasing*.

    Must be initialized with a tree sequence ``ts``, this will serve as
    the history of this first generation of individuals.  The ``samples`` of
    this tree sequence correspond to the individuals in the first generation.
    '''

    def __init__(self, ts, node_ids, time=0.0):
        """
        :param ts TreeSequence: The tree sequence defining relationships between
            the first genreation of individuals.
        :param list node_ids: A dict indexed by input IDs so that
            ``node_ids[k]`` is the node ID of the node corresponding to sample
            ``k`` in the initial ``ts``.  Must specify this for every individual
            that may be a parent moving forward.
        :param float time: The (forwards) time of the "present" at the start of
            the simulation.
        """
        # this is output node ID of the next new node
        self.num_nodes = ts.num_nodes # n
        self.sequence_length = ts.sequence_length
        # this is the largest (forwards) time seen so far
        self.max_time = time  # T
        # last (forwards) time we updated node times
        self.last_update_time = time  # T_0
        # number of nodes that have the time right
        self.last_update_node = ts.num_nodes
        # dict of output node IDs indexed by input labels
        self.node_ids = dict(node_ids)
        # the actual tables that get updated
        #  DON'T actually store ts, just the tables:
        self.nodes = msprime.NodeTable()
        self.edgesets = msprime.EdgesetTable()
        self.sites = msprime.SiteTable()
        self.mutations = msprime.MutationTable()
        self.migrations = msprime.MigrationTable()
        ts.dump_tables(nodes=self.nodes, edgesets=self.edgesets,
                       sites=self.sites, mutations=self.mutations)
        # list of site positions, maintained as site tables don't have
        #   efficient checking for membership
        self.site_positions = {p:k for k, p in enumerate(self.sites.position)}
        # for bookkeeping
        self.num_simplifies = 0

    def __str__(self):
        ret = "\n---------\n"
        ret += "Max time so far:\n"
        ret += str(self.max_time) + "\n"
        ret += "Last update time:\n"
        ret += str(self.last_update_time) + "\n"
        ret += "Node IDs:\n"
        ret += str(self.node_ids) + "\n"
        ret += "Nodes:\n"
        ret += str(self.nodes) + "\n"
        ret += "Edgesets:\n"
        ret += str(self.edgesets) + "\n"
        ret += "Sites:\n"
        ret += str(self.sites) + "\n"
        ret += "Mutations:\n"
        ret += str(self.mutations) + "\n"
        ret += "Migrations:\n"
        ret += str(self.migrations) + "\n"
        ret += "\n---------\n"
        return ret

    def check_ids(self, input_ids):
        """
        Check that all ``input_ids`` are valid.
        """
        for u in input_ids:
            if u not in self.node_ids:
                raise ValueError("Input ID " + str(u) + " not recorded.")

    def get_nodes(self, input_ids):
        """
        Return the output node IDs corresponding to a list of input IDs.
        """
        return [self.node_ids[j] for j in input_ids]

    def add_individual(self, input_id, time,
                       flags=msprime.NODE_IS_SAMPLE,
                       population=msprime.NULL_POPULATION):
        '''
        Add a new individual, recording its (forwards) birth time and assigning
        it a new output Node ID.  New individuals are marked as samples by
        default because we (expect to) have their full history.
        '''
        if input_id not in self.node_ids:
            self.nodes.add_row(flags=flags, population=population,
                               time=time)
            self.node_ids[input_id] = self.num_nodes
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
        if parent not in self.node_ids:
            raise ValueError("Parent " + str(parent) +
                             "'s birth time has not been recorded with " +
                             ".add_individual().")
        out_parent = self.node_ids[parent]
        out_children = tuple([self.node_ids[u] for u in children])
        self.edgesets.add_row(parent=out_parent,
                              children=out_children,
                              left=left,
                              right=right)

    def add_mutation(self, position, node, derived_state, ancestral_state):
        """
        Adds a new mutation to mutation table, and a new site if necessary as well.
        """
        if position not in self.site_positions:
            site = self.sites.num_rows
            self.sites.add_row(position=position, ancestral_state=ancestral_state)
            self.site_positions[position] = site
        else:
            site = self.site_positions[position]
        self.mutations.add_row(site=site, node=node, derived_state=derived_state)

    def update_times(self):
        """
        Update times in the NodeTable.
        """
        dt = self.max_time - self.last_update_time
        times = self.nodes.time
        times[:self.last_update_node] = times[:self.last_update_node] + dt
        times[self.last_update_node:] = self.max_time - times[self.last_update_node:]
        self.nodes.set_columns(flags=self.nodes.flags,
                               population=self.nodes.population,
                               time=times)
        self.last_update_time = self.max_time
        self.last_update_node = self.nodes.num_rows

    def simplify(self, samples):
        """
        Simplifies the underlying tables and returns the resulting tree
        sequence.  To get the tree sequence for a set of samples use
        :meth:``ARGrecorder.tree_sequence`` instead.  `samples` should be a
        list of all "currently living" input individual IDs: i.e., anyone who
        might be a parent or a sample in the future.
        """
        self.check_ids(samples)
        self.update_times()
        msprime.sort_tables(nodes=self.nodes, edgesets=self.edgesets,
                            sites=self.sites, mutations=self.mutations,
                            migrations=self.migrations)
        ts = msprime.load_tables(nodes=self.nodes, edgesets=self.edgesets,
                                 sites=self.sites, mutations=self.mutations,
                                 migrations=self.migrations)
        sample_nodes = self.get_nodes(samples)
        ts = ts.simplify(samples=sample_nodes)
        # update the internal state
        self.num_nodes = ts.num_nodes
        self.last_update_node = ts.num_nodes
        # update index map
        self.node_ids = {k : v for k, v in zip(samples, ts.samples())}
        ts.dump_tables(nodes=self.nodes, edgesets=self.edgesets,
                       sites=self.sites, mutations=self.mutations,
                       migrations=self.migrations)
        self.num_simplifies += 1
        return ts

    def tree_sequence(self, samples=None):
        """
        Return the simplified tree sequence for a given set of input samples,
        *without* simplifying the tables stored internally. (This *does* sort
        them, however.) To simplify the tables stored internally as well, use
        :meth:``ARGrecorder.simplify``.
        """
        if samples is None:
            samples = self.sample_ids()
        else:
            self.check_ids(samples)
        self.update_times()
        msprime.sort_tables(nodes=self.nodes, edgesets=self.edgesets,
                            sites=self.sites, mutations=self.mutations,
                            migrations=self.migrations)
        ts = msprime.load_tables(nodes=self.nodes, edgesets=self.edgesets,
                                 sites=self.sites, mutations=self.mutations,
                                 migrations=self.migrations)
        sample_nodes = self.get_nodes(samples)
        return ts.simplify(samples=sample_nodes)

    def sample_ids(self):
        """
        Return a list of the input IDs corresponding to the samples in the
        internal tables.
        """
        flags = self.nodes.flags
        out = []
        for input_id in self.node_ids:
            j = self.node_ids[input_id]
            if flags[j] & msprime.NODE_IS_SAMPLE:
                out.append(input_id)
        return out

    def dump_sample_table(self, out):
        '''
        Write out the table of info about the samples.
        '''
        flags = self.nodes.flags
        population = self.nodes.population
        time = self.nodes.time
        out.write("id\tflags\tpopulation\ttime\n")
        for j in range(self.nodes.num_rows):
            if flags[j] & msprime.NODE_IS_SAMPLE:
                out.write("{}\t{}\t{}\t{}\n".format(j,
                                                    flags[j],
                                                    population[j],
                                                    time[j]))

    def phony_samples(self, samples, dt=1):
        '''
        Add phony records to the end of the NodeTable that stand in for
        sampling the IDs in `samples`. Will be recorded at living dt units of
        time after the sampled individual.  This is necessary if any of the
        samples are parents to other ones, to avoid ``msprime``'s restriction
        on not sampling internal nodes.  These will be given arbitrary input
        IDs, that will probably break if the simulation is run further.
        '''
        self.check_ids(samples)
        times = self.nodes.time
        populations = self.nodes.population
        sample_nodes = self.get_nodes(samples)
        sample_times = [times[i] for i in sample_nodes]
        sample_populations = [populations[i] for i in sample_nodes]
        new_samples = [None for _ in samples]
        self.sequence_length = max(self.edgesets.right)
        j = 0
        for k in range(len(samples)):
            # find an unused input id
            j = 1 + max(self.node_ids.keys())
            new_samples[k] = j
            self.add_individual(input_id=new_samples[k],
                                time=sample_times[k] + dt,
                                flags=msprime.NODE_IS_SAMPLE,
                                population=sample_populations[k])
            self.add_record(left=0.0,
                            right=self.sequence_length,
                            parent=samples[k],
                            children=(new_samples[k],))
        self.mark_samples(new_samples)
        return new_samples

    def mark_samples(self, samples):
        """
        Mark these individuals as samples internally (but do not simplify).
        """
        self.check_ids(samples)
        sample_nodes = self.get_nodes(samples)
        sample_flag = np.array(msprime.NODE_IS_SAMPLE, dtype='uint32')
        new_flags = self.nodes.flags & ~sample_flag
        new_flags[sample_nodes] |= sample_flag
        self.nodes.set_columns(time=self.nodes.time,
                               population=self.nodes.population,
                               flags=new_flags)
