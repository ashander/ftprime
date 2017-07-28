from .argrecorder import ARGrecorder
import random

# first is 'maternal', second is 'paternal'
mapa_labels = ( 2, 1 )  # NOQA


def ind_to_chrom(ind, mapa):
    '''
    Returns the unique *chromosome ID* corresponding to
    the chromosome of individual 'ind' inherited from
    parent 'mapa' (either 2=maternal or 1=paternal).  (Chromosome IDs are ints.)
    This takes
       0, paternal -> 0
       0, maternal -> 1
       1, paternal -> 2
       1, maternal -> 3
       2, paternal -> 4
       2, maternal -> 5
    ...
    '''
    return int(2 * ind + mapa - 1)


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

    def __init__(self, ts, node_ids, locus_position):
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
        """
        self.sequence_length = ts.sequence_length
        self.locus_position = locus_position
        self.last_child = -1
        self.time = 0.0

        if locus_position[0] != 0.0 or locus_position[-1] != self.sequence_length:
            raise ValueError("locus_position (and lociPos) must include a locus\
                              at each end of the chromosome.")

        haploid_node_ids = {self.i2c(x[0], x[1]):node_ids[(x[0], x[1])] 
                            for x in node_ids}
        self.args = ARGrecorder(ts, node_ids=haploid_node_ids)
        # will record IDs of diploid samples here when they are chosen
        # but note we don't keep anything else about them here (time, location)
        # as this is recorded by the ARGrecorder
        self.diploid_samples = None

    def i2c(self, k, p):
        # individual ID to chromsome ID
        out = ind_to_chrom(int(k), mapa_labels[p])
        return out

    def increment_time(self):
        self.time += 1.0

    def collect_recombs(self, lines):
        for line in lines.strip().split('\n'):
            # print("A: "+line)
            # child, parent, ploid,*rec = [int(x) for x in line.split()]
            linex = [int(x) for x in line.split()]
            child = linex[0]
            parent = linex[1]
            ploid = linex[2]
            rec = linex[3:]
            # lines come in pairs: maternal/paternal.
            if child == self.last_child:
                child_p = 1
            else:
                child_p = 0
                self.last_child = child
            # if self.ind_to_time(child) > self.ind_to_time(parent):
            #     raise ValueError(str(child)+" at "+
            #            str(self.ind_to_time(child))+'
            #            " does not come after " + str(parent)+
            #            " at "+str(self.ind_to_time(parent)))

            start = 0.0
            child_chrom = self.i2c(child, child_p)
            # print("  existing parent:",parent, ploid,"-i2c->",
            # self.i2c(parent, ploid))

            # print("  adding child:",child, child_p,"-i2c->",child_chrom,
            # self.time)
            self.args.add_individual(child_chrom, self.time)
            for r in rec:
                # do this check to avoid a simuPOP bug
                if r < len(self.locus_position) - 1:
                    breakpoint = random.uniform(self.locus_position[r],
                                                self.locus_position[r + 1])
                    # print("--- ",start, self.locus_position[r],
                    #       "< = ",breakpoint,"<=",self.locus_position[r+1])
                    self.args.add_record(
                            left=start,
                            right=breakpoint,
                            parent=self.i2c(parent, ploid),
                            children=(child_chrom,))
                    start = breakpoint
                    ploid = ((ploid + 1) % 2)
            # print("--- ",start, self.sequence_length," |")
            self.args.add_record(
                    left=start,
                    right=self.sequence_length,
                    parent=self.i2c(parent, ploid),
                    children=(child_chrom,))

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

    def add_diploid_samples(self, nsamples, sample_ids):
        """
        Choose randomly a set of individuals, label these as samples,
        and simplify the tree sequence.  This is irreversible.
        """
        sample_indices = random.sample(range(len(sample_ids)), nsamples)
        self.diploid_samples = [int(sample_ids[k]) for k in sample_indices]
        self.simplify_to_samples()

    def simplify_to_samples(self):
        """
        Label samples and simplify the tree sequence.  This is irreversible.
        """
        chrom_samples = [ind_to_chrom(x, a)
                         for x in self.diploid_samples
                         for a in mapa_labels]
        self.args.simplify(samples=chrom_samples)
