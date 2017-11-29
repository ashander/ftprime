from .argrecorder import (
    ARGrecorder,
    edge_dt,
    node_dt,
)
import msprime
import random
import numpy as np
import time as timer
from .benchmarker import Timings


class RecombCollector:
    '''
    Collect and parse recombination events as output by simuPOP's Recombinator,
    which outputs like:

    offspringID parentID startingPloidy rec1 rec2 ....

    ... coming in *pairs*.

    This keeps track of *time* - so when used, the time must be updated -
    in simuPOP, by adding rc.increment_time() to the PreOps.  It should be in PreOps
    because if so
        - the initial generation is recorded at time 0.0
        - calling increment_time() increases the time by 1
        - the first generation is recorded at time 1.0
        - ...

    '''
    def __init__(self, ts, node_ids, locus_position, benchmark=False,
                 mode='text'):
        """
        :param TreeSequence ts: A tree sequence describing the history of each
            chromosome in the population before the simulation starts.
        :param dict node_ids: A dict indexed by (individual ID, ploidy)
            ``node_ids[(k,0)]`` is the node ID of the node corresponding to the
            maternally inherited chromosome of sample ``k`` in the initial ``ts``,
            and ``node_ids[(k,1)]`` is the node ID for the paternal chromosome.
            Must specify this for every individual that may be a parent moving forward.
        :param list locus_position: A list of coordinates on the genome of the loci
            that simuPOP is keeping track of.  There must be a locus at the beginning
            and also at the end of the chromosome.
            :param bool benchmark: Whether to store benchmark information in the
            ARGrecorder.
        :param str mode: can be 'text or 'binary' then bstrs must be passed to
            `.collect_recombs`.

        """
        if mode == 'text':
            self.split = '\n'
        elif mode == 'binary':
            self.split = b'\n'
        else:
            raise ValueError("mode must be 'str' or 'binary'")
        self.sequence_length = ts.sequence_length
        self.locus_position = locus_position
        self.last_child = -1
        self.time = 0.0

        if locus_position[0] != 0.0 or locus_position[-1] != self.sequence_length:
            raise ValueError("locus_position (and lociPos) must include a locus\
                              at each end of the chromosome.")

        haploid_node_ids = {self.i2c(x[0], x[1]):node_ids[(x[0], x[1])] 
                            for x in node_ids}
        if not benchmark:
            self.args = ARGrecorder(node_ids=haploid_node_ids, ts=ts)
        else:
            self.args = ARGrecorder(node_ids=haploid_node_ids, ts=ts,
                                    timings=Timings())

        # will record IDs of diploid samples here when they are chosen
        # but note we don't keep anything else about them here (time, location)
        # as this is recorded by the ARGrecorder
        self.diploid_samples = None

    @property
    def mode(self):
        if self.split == '\n':
            return 'text'
        elif self.split == b'\n':
            return 'binary'

    def i2c(self, k, p):
        """
        Get the "chromosome ID", used in the underlying ARGrecorder.

        :param int k: The individual ID.
        :param int p: The chromosome index (0 or 1).

        :return numpy.int32: The "input ID" used by the underlying ARGrecorder
        corresponding to this chromsome.
        """
        # individual ID to chromsome ID
        if not (p == np.int32(0) or p == np.int32(1)):
            raise ValueError("Chromosome ID must be 0 (paternal) or 1 (maternal).")
        return np.int32(2) * k + p

    def i2n(self, k, p):
        """
        Get the output node ID used in the underlying tree sequence tables.

        :param int k: The individual ID.
        :param int p: The chromosome index (0 or 1).

        :return int: The node ID for this chromosome in the output tables.
        """
        return self.args.node_ids[self.i2c(k,p)]

    def increment_time(self):
        self.time += 1.0

    def collect_recombs(self, lines):
        """
        Collects recombinations arriving in text form like:
            offspringID parentID startingPloidy rec1 rec2 ....
        in *pairs* for (paternal, maternal) chromosomes, as output by
        ``simuPOP.Recombinator()``. A parental chromosome inherited without a
        crossover would be recorded with no recombinations.
        """
        if self.args.timings is not None:
            before = timer.process_time()
        nodes, edges_data = zip(*self._nodes_edges_iter(lines))
        node_ids, nodes_data = zip(*nodes)
        # Below we use flatten to avoid [[(el1), (el2), ... (el3)]]
        nd = np.column_stack(nodes_data).flatten()
        self.args.add_individuals(node_ids, nd)
        self.args.add_records([edge
                              for edges in edges_data
                              for edge in edges])
        if self.args.timings is not None:
            self.args.timings.time_appending += timer.process_time() - before

    def _nodes_edges_iter(self, lines):
        for line in lines.strip().split(self.split):
            # print("A: "+line)
            # child, parent, ploid,*rec = [int(x) for x in line.split()]
            linex = np.fromstring(line, sep=' ', dtype=np.int32)
            child = linex[0]
            parent = linex[1]
            ploid = linex[2]
            rec = linex[3:]
            # lines come in pairs: maternal/paternal.
            if child == self.last_child:
                child_p = np.int32(1)
            else:
                child_p = np.int32(0)
                self.last_child = child
            # if self.ind_to_time(child) > self.ind_to_time(parent):
            #     raise ValueError(str(child)+" at "+
            #            str(self.ind_to_time(child))+'
            #            " does not come after " + str(parent)+
            #            " at "+str(self.ind_to_time(parent)))

            start = 0.0
            child_chrom = self.i2c(child, child_p)
            assert type(child_chrom) is np.int32, str(type(child_chrom))
            # print("  existing parent:",parent, ploid,"-i2c->",
            # self.i2c(parent, ploid))

            # print("  adding child:",child, child_p,"-i2c->",child_chrom,
            # self.time)
            ndata = (child_chrom,
                     np.array(
                         (msprime.NODE_IS_SAMPLE, self.time, msprime.NULL_POPULATION),
                         dtype=node_dt))
            ep_iter = ((self.locus_position[r], self.locus_position[r + 1])
                       for r in rec
                       # do this check to avoid a simuPOP bug
                       if r < len(self.locus_position) - 1)
            edata = []
            for i, eps in enumerate(ep_iter):
                breakpoint = random.uniform(*eps)
                # print("--- ",start, self.locus_position[r],
                #       "< = ",breakpoint,"<=",self.locus_position[r+1])
                edata.append((start,
                              breakpoint,
                              self.i2c(parent, ploid),
                              child_chrom))
                start = breakpoint
                ploid = ((ploid + 1) % 2)
            # print("--- ",start, self.sequence_length," |")
            edata.append((start,
                          self.sequence_length,
                          self.i2c(parent, ploid),
                          child_chrom))
            yield ndata, edata

    def tree_sequence(self, samples):
            """
            Returns a tree sequence, that retains only information relevant
            to the diploid individuals listed in `samples`.

            :param list samples: A list of diploid input individual IDs.
            """
            haploid_ids = [self.i2c(i,p) for i in samples for p in (0,1)]
            return self.args.tree_sequence(haploid_ids)

    def simplify(self, samples):
        """
        Simplify the underlying tree sequence, retaining only information relevant
        to the diploid individuals listed in `samples`.

        :param list samples: A list of diploid input individual IDs.
        """
        haploid_ids = [self.i2c(i,p) for i in samples for p in (0,1)]
        self.args.simplify(haploid_ids)

    def add_locations(self, input_ids, locations):
        """
        Assign the `population` field of each individual in `input_ids` to the corresponding
        entry in `locations`.

        :param list input_ids: A list of input diploid individual IDs.
        :param list locations: A list of population IDs.
        """
        populations = self.args.nodes.population
        for i, loc in zip([int(j) for j in input_ids], locations):
            for p in (0,1):
                h = self.i2c(i, p)
                node_id = self.args.node_ids[h]
                populations[node_id] = int(loc)
        self.args.nodes.set_columns(flags=self.args.nodes.flags,
                                    time=self.args.nodes.time,
                                    population=populations)

    def dump_sample_table(self, out):
        """
        TODO record here somehow which chromosomes here are in the same individual
        """
        self.args.dump_sample_table(out)
