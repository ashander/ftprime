import msprime
import itertools
import time as timer  # otherwise name clash
import numpy as np

NULL_ID = -1

node_dt = np.dtype([('flags', np.uint32),
                    ('time', np.float),
                    ('population', np.int32)])

edge_dt = np.dtype([('left', np.float),
                    ('right', np.float),
                    ('parent', np.int32),
                    ('child', np.int32)])


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

        b. how their chromosome was inherited from her parents with ``add_record``,

        c. and any new mutations with ``add_mutation``. 

    2. Periodically, run ``simplify(samples)`` to remove unnecessary
        information from the recorded tables.  ``samples`` should be a list of
        input IDs of all individuals whose history may be needed in the future:
        the current generation, and any ancestral samples.

    Note: at any time, individuals for whom we have complete information are
    marked as samples in the Node Table; however, this is not consulted when
    calling ``simplify``.

    '''
    __nodes = None
    __edges = None

    def __init__(self, node_ids=None, ts=None, time=0.0,
                 sequence_length=None, timings=None):
        """
        The tree sequence ``ts`` passed in defines history before the
        simulation begins.  If this is missing, then the input IDs specified in
        ``node_ids`` must be ``0...n-1``.

        :param dict node_ids: A dict indexed by input IDs so that
            ``node_ids[k]`` is the node ID of the node corresponding to sample
            ``k`` in the initial ``ts``.  Must specify this for every individual
            that may be a parent moving forward.
        :param TreeSequence ts: The tree sequence containing prior history.
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
        self.__nodes = msprime.NodeTable()
        self.__edges = msprime.EdgeTable()
        self.__sites = msprime.SiteTable()
        self.__mutations = msprime.MutationTable()
        self.__migrations = msprime.MigrationTable()
        if ts is None:
            for j, k in enumerate(sorted(self.node_ids.keys())):
                assert j == self.node_ids[k]
            default = (msprime.NODE_IS_SAMPLE, time, msprime.NULL_POPULATION)
            nt = np.fromiter(zip(itertools.cycle((default, ))),
                             count=len(self.node_ids), dtype=node_dt)
            self.__nodes.set_columns(flags=nt['flags'],
                                     time=nt['time'],
                                     population=nt['population'])
            self.last_sorted_edge = 0
        if ts is not None:
            ts.dump_tables(nodes=self.__nodes, edges=self.__edges,
                           sites=self.__sites, mutations=self.__mutations,
                           migrations=self.__migrations)
            self.last_sorted_edge = self.edges.num_rows
        if sequence_length is not None:
            if ts is not None:
                if sequence_length != ts.sequence_length:
                    raise ValueError("Provided sequence_length does not match",
                                     "that of tree sequence ts.")
            self.sequence_length = sequence_length
        elif ts is not None:
            self.sequence_length = ts.sequence_length
        else:
            raise ValueError("If prior history is not specified, sequence",
                             "length must be provided.")
        # last (forwards) time we updated node times
        self.last_update_time = time  # T_0
        # number of nodes that have the time right
        self.last_update_node = self.nodes.num_rows
        # list of site positions, maintained as site tables don't have
        #   efficient checking for membership
        self.site_positions = {p:k for k, p in enumerate(self.sites.position)}
        # for bookkeeping
        self.num_simplifies = 0
        if self.timings is not None:
            self.timings.time_prepping += timer.process_time() - start

    def __repr__(self):
        ret = "\n---------\n"
        ret += ' '.join(["Nodes:", repr(self.nodes), "\n"])
        ret += ' '.join(["Edges:", repr(self.edges), "\n"])
        ret += ' '.join(["Sites:", repr(self.sites), "\n"])
        ret += ' '.join(["Mutations:", repr(self.mutations), "\n"])
        # ret += "\n---------\n"
        return ret

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
        ret += "Edges:\n"
        ret += str(self.edges) + "\n"
        ret += "Sites:\n"
        ret += str(self.sites) + "\n"
        ret += "Mutations:\n"
        ret += str(self.mutations) + "\n"
        ret += "Migrations:\n"
        # ret += str(self.migrations) + "\n"
        # ret += "\n---------\n"
        return ret

    def __call__(self, parents, times, populations, childs, lefts, rights):
        """
        Does both ``add_individuals()`` and ``add_records()`` steps.

        :param iterable of int parent: The input IDs of the parent.
        :param iterable of int populations: The population ID of birth of the
            children.
        :param iterable of float times: The time of birth of the children.
        :param iterable of int childs: An iterable of input IDs of the
            children.
        :param iterable of float left: The left endpoint of the chromosomal
            segment that child inherits from parent.
        :param iterable of float right: The right endpoint of the chromosomal
            segment that child inherits from parent.
        """

        def _unique(childs, times, pops):
            seen = set()
            for child, time, pop in zip(childs, times, pops):
                if child not in seen:
                    yield (child, time, pop)
                seen.add(child)

        new_nodes, new_times, new_pops = zip(*_unique(childs,
                                                      times,
                                                      populations))
        self.add_individuals(input_ids=new_nodes, times=new_times,
                             populations=new_pops)
        self.add_records(lefts=lefts, rights=rights, parents=parents,
                         childs=childs)

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

    def add_individuals(self, input_ids, times, flagss=None, populations=None):
        '''
        Add new individuals, resulting in records as in `.add_individual`, but
        more efficiently.

        :param iterable of int input_ids: The input ID of the new individual.
        :param iterable of float times: The time of birth of the individual.
        :param iterable of int flagss: Any msprime flags to record.
        :param iterable of int populations: The population ID of birth of the
            individual.

        :raises ValueError for input_ids that have already been entered
        '''
        if flagss is None:
            flagss = itertools.cycle((msprime.NODE_IS_SAMPLE, ))
        if populations is None:
            populations = itertools.cycle((msprime.NULL_POPULATION, ))

        # inflate to a list so we have a length; then set up the map from
        # node_ids to tree sequence samples
        num_new_nodes = 0
        for u, n in zip(input_ids, itertools.count(start=self.nodes.num_rows)):
            if u in self.node_ids:
                raise ValueError("Input ID " + str(u) +
                                 " already exits, and cannot be added as new.")
            self.node_ids[u] = n
            num_new_nodes += 1

        nt = np.fromiter(zip(flagss, times, populations), dtype=node_dt,
                         count=num_new_nodes)
        self.__nodes.append_columns(flags=nt['flags'],
                                    time=nt['time'],
                                    population=nt['population'])
        self.max_time = max(self.max_time, max(times))

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
        self.add_individuals((input_id, ), (time, ), (flags, ), (population, ))

    def add_records(self, lefts, rights, parents, childs):
        '''
        Add multiple records corresponding to a reproduction event in which
        childs (each a single ID), inherit from parent (also each  single ID) on
        the interval [left, right).

        :param iterable of float left: The left endpoint of the chromosomal
            segment inherited.
        :param iterable of float right: The right endpoint of the chromosomal
            segment inherited.
        :param iterable of int parent: The input ID of the parent.
        :param iterable of int: An iterable of input IDs of the
            children.
        '''
        def _valid_parents(parents):
            for p in parents:
                if p not in self.node_ids:
                    raise ValueError("Parent " + str(p) +
                                     "'s birth time has not been recorded with " +
                                     ".add_individual().")
                yield p

        out_parents = (self.node_ids[parent]
                       for parent in _valid_parents(parents))
        out_childrens = (self.node_ids[u] for u in childs)
        edges = zip(lefts, rights, out_parents, out_childrens)
        et = np.fromiter(edges, dtype=edge_dt)
        self.__edges.append_columns(parent=et['parent'],
                                    child=et['child'],
                                    left=et['left'],
                                    right=et['right'])

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
        kids = tuple(c for c in children)
        self.add_records(itertools.repeat(left, len(kids)),
                         itertools.repeat(right, len(kids)),
                         itertools.repeat(parent, len(kids)),
                         kids)

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
            site = self.sites.num_rows
            self.sites.add_row(position=position, ancestral_state=ancestral_state)
            self.site_positions[position] = site
        else:
            site = self.site_positions[position]
        self.mutations.add_row(site=site, node=self.node_ids[node], 
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
        times = self.nodes.time
        times[:self.last_update_node] = times[:self.last_update_node] + dt
        times[self.last_update_node:] = self.max_time - times[self.last_update_node:]
        self.__nodes.set_columns(flags=self.nodes.flags,
                                 population=self.nodes.population,
                                 time=times)
        self.last_update_time = self.max_time
        self.last_update_node = self.nodes.num_rows

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
            before = timer.process_time()

        if self.last_sorted_edge < self.edges.num_rows:
            # begin modified block from @jeromekelleher and @molpopgen
            # Copy the already sorted edges to local arrays
            left = self.edges.left[:self.last_sorted_edge]
            right = self.edges.right[:self.last_sorted_edge]
            parent = self.edges.parent[:self.last_sorted_edge]
            child = self.edges.child[:self.last_sorted_edge]
            # Get the new edges and reverse them. After this, we know that all edges
            # are correctly sorted with respect to time. We then sort each time slice
            # individually, reducing the overall cost of the sort.
            new_left = self.edges.left[:self.last_sorted_edge:-1]
            new_right = self.edges.right[:self.last_sorted_edge:-1]
            new_parent = self.edges.parent[:self.last_sorted_edge:-1]
            new_child = self.edges.child[:self.last_sorted_edge:-1]
            parent_time = self.__nodes.time[new_parent]
            breakpoints = np.where(parent_time[1:] != parent_time[:-1])[0] + 1
            self.__edges.reset()
            if self.timings is not None:
                self.timings.time_appending += timer.process_time() - before
                before = timer.process_time()
            start = 0
            for end in itertools.chain(breakpoints, [-1]):
                assert np.all(parent_time[start: end] == parent_time[start])
                self.__edges.append_columns(left=new_left[start: end],
                                            right=new_right[start: end],
                                            parent=new_parent[start: end],
                                            child=new_child[start: end])
                msprime.sort_tables(nodes=self.__nodes,
                                    edges=self.__edges,
                                    sites=self.__sites,
                                    mutations=self.__mutations,
                                    edge_start=start)
                start = end
            if self.timings is not None:
                self.time_sorting += timer.process_time() - before
                before = timer.process_time()

            # NOTE: matches fwdpy11_argexample code from
            # https://github.com/molpopgen/fwdpy11_arg_example/pull/8
            # Append the old sorted edges to the table.
            self.__edges.append_columns(left=left, right=right, parent=parent,
                                        child=child)
        if self.timings is not None:
            before = timer.process_time()
        msprime.simplify_tables(samples=sample_nodes, nodes=self.__nodes,
                                edges=self.__edges, sites=self.__sites,
                                mutations=self.__mutations,
                                sequence_length=self.sequence_length)
        #                       migrations=self.migrations)
        self.last_sorted_edge = self.edges.num_rows
        # end modified block from @jeromekelleher and @molpopgen
        if self.timings is not None:
            self.timings.time_simplifying += timer.process_time() - before
        # update the internal state
        self.last_update_node = self.nodes.num_rows
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
        msprime.sort_tables(nodes=self.__nodes, edges=self.__edges,
                            sites=self.__sites, mutations=self.__mutations)
        #                   migrations=self.migrations)
        self.last_sorted_edge = self.edges.num_rows
        self.mark_samples(samples)
        if self.timings is not None:
            self.timings.time_sorting += start - timer.process_time()
        ts = msprime.load_tables(nodes=self.nodes, edges=self.edges,
                                 sites=self.sites, mutations=self.mutations,
                                 sequence_length=self.sequence_length)
        #                        migrations=self.migrations)
        sample_nodes = self.get_nodes(samples)
        return ts.simplify(samples=sample_nodes)

    def sample_ids(self):
        """
        Return a list of the input IDs corresponding to the samples in the
        internal tables.

        :return list: A list of input IDs.
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
        new_flags = self.nodes.flags & ~sample_flag
        new_flags[sample_nodes] |= sample_flag
        self.__nodes.set_columns(time=self.nodes.time,
                                 population=self.nodes.population,
                                 flags=new_flags)

    @property
    def nodes(self):
        return self.__nodes

    @property
    def edges(self):
        return self.__edges

    @property
    def sites(self):
        return self.__sites

    @property
    def mutations(self):
        return self.__mutations

    @property
    def migrations(self):
        return self.__migrations
