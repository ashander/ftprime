import msprime
import time as timer  # otherwise name clash
import numpy as np

NULL_ID = -1

def null_tree_sequence():
    return msprime.load_tables(nodes=msprime.NodeTable(),
                               edges=msprime.EdgeTable())

class ARGrecorder(object):
    '''
    To record the ARG, this keeps track of
        - a NodeTable and an EdgeTable
        - information about what has happened since the last simplification
        - a mapping from (input) individual IDs to (output) node IDs
    The reason for maintaining this mapping, rather than having (input)
    individual IDs always equal to the (output) node IDs is to allow periodic
    simplification, which decouples the two.

    Between simplification steps, the EdgeTable is appended to, so it is
    always "up to date".  However, the NodeTable is *not* kept up to date,
    because its `time` fields are recorded in *time ago*; we also keep track of
        - a list of birth times of individual IDs
    which are translated to time-ago at each simplification step, and appended
    to the Node Table.

    The internal state is stored using
        - ``self.node_ids[k]`` : the output Node ID corresponding to the input
          individual ID ``k``.

    Must be initialized with a set of tables which will serve as the history of
    this first generation of individuals.

    To use an ARGrecorder instance to record the ARG coming from an
    forwards-time simulator:

    0. Initialize with a tree sequence describing history before the simulation
        begins, and the mapping from input IDs (in the forwards-time simulator) to
        output IDs (node IDs in the tree sequence).

    1. Each time a new individual (chromosome) is born, record:
        a. their birth time with ``add_individual``,

        b. how their chromosome was inherited from her parents with ``add_record``.

    2. Periodically, run ``simplify(samples)`` to remove unnecessary
        information from the recorded tables.  ``samples`` should be a list of
        input IDs of all individuals whose history may be needed in the future:
        the current generation, and any ancestral samples.

    Note: at any time, individuals for whom we have complete information are
    marked as samples in the Node Table; however, this is not consulted when
    calling ``simplify``.

    '''

    def __init__(self, node_ids=None, tables=None, ts=None, time=0.0,
                 sequence_length=None, timings=None):
        """
        The tables passed in define history before the simulation begins.  If
        these are missing, then the input IDs specified in ``node_ids`` must be
        ``0...n-1``.

        :param dict node_ids: A dict indexed by input IDs so that
            ``node_ids[k]`` is the node ID of the node corresponding to sample
            ``k`` in the initial ``ts``.  Must specify this for every individual
            that may be a parent moving forward.
        :param TableCollection tables: A collection of tables describing
            prehistory of the simulation.
        :param TreeSequence ts: An alternative method to specifying past history.
        :param float time: The (forwards) time of the "present" at the start of
            the simulation.
        :param float sequence_length: The total length of the sequence (derived
            from input if not provided).
        :param ftprime.benchmarker.Timings timings:  An object to record timing
        information.
        """
        if timings is not None:
            self.timings = timings
            start = timer.process_time()
        else:
            self.timings = None

        # this is the largest (forwards) time seen so far
        self.max_time = time  # T
        # dict of output node IDs indexed by input labels
        if node_ids is None:
            self.node_ids = {}
        else:
            self.node_ids = dict(node_ids)
        # the actual tables that get updated
        #  DON'T actually store ts, just the tables:
        if ts is not None:
            tables = ts.dump_tables()
        elif tables is None:
            tables = msprime.TableCollection(sequence_length=sequence_length)
            if ts is None:
                for j, k in enumerate(sorted(self.node_ids.keys())):
                    assert j == self.node_ids[k]
                    tables.nodes.add_row(population=msprime.NULL_POPULATION, time=time)
        self.tables = tables
        if sequence_length is not None:
            if ts is not None:
                if sequence_length != ts.sequence_length:
                    raise ValueError("Provided sequence_length does not match",
                                     "that of tree sequence ts.")
            self.sequence_length = sequence_length
        elif ts is not None:
            self.sequence_length = ts.sequence_length
        else:
            if edges.num_rows > 0:
                self.sequence_length = max(edges.right)
            else:
                raise ValueError("If prior history is not specified, sequence",
                                 "length must be provided.")
        # last (forwards) time we updated node times
        self.last_update_time = time  # T_0
        # number of nodes that have the time right
        self.last_update_node = self.tables.nodes.num_rows
        # list of site positions, maintained as site tables don't have
        #   efficient checking for membership
        self.site_positions = {p:k for k, p in enumerate(self.tables.sites.position)}
        # for bookkeeping
        self.num_simplifies = 0
        if self.timings is not None:
            self.timings.time_prepping += timer.process_time() - start

    def __str__(self):
        ret = "\n---------\n"
        ret += "Max time so far:\n"
        ret += str(self.max_time) + "\n"
        ret += "Last update time:\n"
        ret += str(self.last_update_time) + "\n"
        ret += "Node IDs:\n"
        ret += str(self.node_ids) + "\n"
        ret += "Tables:\n"
        ret += str(self.tables) + "\n"
        return ret

    def __call__(self, parent, time, population, child, left, right):
        """
        Does both ``add_individual()`` and ``add_record steps()``.

        :param int parent: Input ID of the parent.
        :param float time: The time of birth of the child.
        :param int population: The population ID where the child is born.
        :param int child: Input ID of the child.
        :param float left: Left end of the segment that child inherits from parent.
        :param float right: Right end of the segment that child inherits from parent.
        """
        if child not in self.node_ids:
            self.add_individual(input_id=child, time=time, population=population)
        self.add_record(left=left, right=right, parent=parent, children=(child,))

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

        :param int input_id: The input ID of the new individual.
        :param float time: The time of birth of the individual.
        :param flags int: Any msprime flags to record (probably not needed).
        :param population int: The population ID of birth of the indivdiual
            (may be omitted).  
        '''
        if (population != msprime.NULL_POPULATION 
                and population > self.tables.populations.num_rows):
            if population < 0:
                raise ValueError("Illegal population: " + str(population))
            while self.tables.populations.num_rows <= population:
                self.tables.populations.add_row()
        if input_id not in self.node_ids:
            self.node_ids[input_id] = self.tables.nodes.num_rows
            self.tables.nodes.add_row(flags=flags, population=population,
                                      time=time)
            self.max_time = max(self.max_time, time)
        else:
            # nothing bad happens if we try to add an individual more than once,
            # but this is helpful for debugging
            raise ValueError("Attempted to add " + str(input_id) +
                             ", who already exits, as a new individual.")

    def add_record(self, left, right, parent, children):
        '''
        Add records corresponding to a reproduction event in which children (a
        tuple of IDs) inherit from parent (a single ID) on the interval
        [left,right).

        :param float left: The left endpoint of the chromosomal segment inherited.
        :param float right: The right endpoint of the chromosomal segment inherited.
        :param int parent: The input ID of the parent.
        :param list children: An iterable of input IDs of the children.
        '''
        if parent not in self.node_ids:
            raise ValueError("Parent " + str(parent) +
                             "'s birth time has not been recorded with " +
                             ".add_individual().")
        out_parent = self.node_ids[parent]
        out_children = tuple([self.node_ids[u] for u in children])
        for child in out_children:
            self.tables.edges.add_row(parent=out_parent,
                                      child=child,
                                      left=left,
                                      right=right)
    def add_mutation(self, position, node, derived_state, ancestral_state):
        """
        Adds a new mutation to mutation table, and a new site if necessary as well.
        :param float position: The chromosomal position of the mutation.
        :param int node: The input ID of the individual on whose chromosome the
            mutation occurred.  
        :param string derived_state: The allele resulting from the mutation.
        :param string ancestral_state: The original allele that the mutation
            replaces (only used if this is the first mutation at this position).
        """
        if position not in self.site_positions:
            site = self.tables.sites.num_rows
            self.tables.sites.add_row(position=position, ancestral_state=ancestral_state)
            self.site_positions[position] = site
        else:
            site = self.site_positions[position]
        self.tables.mutations.add_row(site=site, node=self.node_ids[node], 
                                      derived_state=derived_state)

    def update_times(self):
        """
        Update the times in the NodeTable.  This is necessary because input
        times, as recorded in the NodeTable, are in forwards time, but the
        NodeTable must be in reverse time (time since the end of the
        simulation).  Therefore, this needs to (a) add an increment to any
        already-updated times in the NodeTable, and (b) reverse any times added
        since the last update.
        """
        dt = self.max_time - self.last_update_time
        times = self.tables.nodes.time
        times[:self.last_update_node] = times[:self.last_update_node] + dt
        times[self.last_update_node:] = self.max_time - times[self.last_update_node:]
        self.tables.nodes.set_columns(flags=self.tables.nodes.flags,
                                      population=self.tables.nodes.population,
                                      time=times)
        self.last_update_time = self.max_time
        self.last_update_node = self.tables.nodes.num_rows

    def simplify(self, samples):
        """
        Simplifies the underlying tables.  `samples` should be a list of all
        "currently living" input individual IDs: i.e., anyone who might be a
        parent or a sample in the future.

        Note: to get the tree sequence for a set of samples use
        :meth:``ARGrecorder.tree_sequence``.

        This resets the underlying map from input IDs to output IDs. The
        individals in ``samples`` will be assigned output IDs
        ``0,...,len(samples)-1``.

        :param list samples: A list of the input IDs whose entire history
            should be kept; information not relevant to the history of these
            samples will be discarded.
        """
        self.check_ids(samples)
        self.update_times()
        sample_nodes = self.get_nodes(samples)
        if self.timings is not None:
            start = timer.process_time()
        self.tables.sort()
        if self.timings is not None:
            start2 = timer.process_time()
            self.timings.time_sorting += start2 - start
        self.tables.simplify(sample_nodes)
        if self.timings is not None:
            self.timings.time_simplifying += timer.process_time() - start2
        # update the internal state
        self.last_update_node = self.tables.nodes.num_rows
        # update index map: sample[k] now maps to k
        self.node_ids = {k : v for v, k in enumerate(samples)}
        self.num_simplifies += 1

    def tree_sequence(self, samples=None):
        """
        Return the simplified tree sequence for a given set of input samples,
        *without* simplifying the tables stored internally. (This *does* sort
        them and label samples, however.) To simplify the tables stored
        internally as well, use :meth:``ARGrecorder.simplify``.

        :param list samples: A list of the input IDs whose history is recorded
            in the resulting tree sequence.  If this is missing, all available
            individuals will be used.
        :return TreeSequence: The simplified tree sequence recording the
            history of ``samples``; in this tree sequence, ``sample[k]``
            corresponds to Node ID ``k``.
        """
        if samples is None:
            samples = self.sample_ids()
        else:
            self.check_ids(samples)
        self.update_times()
        if self.timings is not None:
            start = timer.process_time()
        self.tables.sort()
        self.mark_samples(samples)
        if self.timings is not None:
            self.timings.time_sorting += start - timer.process_time()
        ts = self.tables.tree_sequence()
        sample_nodes = self.get_nodes(samples)
        return ts.simplify(samples=sample_nodes)

    def sample_ids(self):
        """
        Return a list of the input IDs corresponding to the samples in the
        internal tables.

        :return list: A list of input IDs.
        """
        flags = self.tables.nodes.flags
        out = []
        for input_id in self.node_ids:
            j = self.node_ids[input_id]
            if flags[j] & msprime.NODE_IS_SAMPLE:
                out.append(input_id)
        return out

    def mark_samples(self, samples):
        """
        Mark these individuals as samples internally (but do not simplify).

        :param list samples: A list of input IDs that should be marked as with
            the msprime flag for samples in the underlying NodeTable.  This
            does not affect what happens at the next ``simplify``, and is
            provided mainly for convenience.  
        """
        self.check_ids(samples)
        sample_nodes = self.get_nodes(samples)
        sample_flag = np.array(msprime.NODE_IS_SAMPLE, dtype='uint32')
        new_flags = self.tables.nodes.flags & ~sample_flag
        new_flags[sample_nodes] |= sample_flag
        self.tables.nodes.set_columns(time=self.tables.nodes.time,
                                      population=self.tables.nodes.population,
                                      flags=new_flags)
