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
    ... coming in *pairs*

    This needs:
    namples - number of *diploid* samples
    first_gen - list of individuals IDs in the initial population
    ancestor_age - number of generations before beginning of simulation that common ancestor lived
    length - length of chromosome
    locus_position - list of positions of the loci simuPOP refers to along the chromosome

    This keeps track of *time* - so when used, the time must be updated -
    in simuPOP, by adding rc.increment_time() to the PreOps.  It should be in PreOps
    because if so
        - the initial generation is recorded at time 0.0
        - calling increment_time() increases the time by 1
        - the first generation is recorded at time 1.0
        - ...
    '''

    def __init__(self, first_gen, ancestor_age, length,
                 locus_position):
        self.ancestor_age = ancestor_age
        self.length = length
        self.locus_position = locus_position
        self.last_child = -1
        self.time = 0.0

        if locus_position[0] != 0.0 or locus_position[-1] != length:
            raise ValueError("locus_position (and lociPos) must include a locus at each end of the chromosome.")

        self.args = ARGrecorder()

        self.universal_ancestor = 0
        # will record IDs of diploid samples here when they are chosen
        # but note we don't keep anything else about them here (time, location)
        # as this is recorded by the ARGrecorder
        self.diploid_samples = None
        self.args.add_individual(name=self.universal_ancestor,
                                 time= (-1) * self.ancestor_age)
        # add initial generation
        first_haps = [self.i2c(k, p) for k in first_gen for p in [0, 1]]
        first_haps.sort()
        self.args.add_record(
                left=0.0,
                right=self.length,
                parent=self.universal_ancestor,
                children=tuple(first_haps))
        for k in first_haps:
            # print("Adding:",k, p, self.time)
            self.args.add_individual(k, self.time)

    def i2c(self, k, p):
        # individual ID to chromsome ID
        # 1+ for the universal ancestor at slot 0
        out = 1 + ind_to_chrom(k, mapa_labels[p])
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
            # print("--- ",start, self.length," |")
            self.args.add_record(
                    left=start,
                    right=self.length,
                    parent=self.i2c(parent, ploid),
                    children=(child_chrom,))

    def add_diploid_samples(self, nsamples, sample_ids, populations):
        # sample_ids is the list of diploid IDs to draw the samples from
        # NOTE: does NOT remove previous samples
        assert(len(sample_ids) == len(populations))
        sample_indices = random.sample(range(len(sample_ids)), nsamples)
        self.diploid_samples = [sample_ids[k] for k in sample_indices]
        # print("Samples ("+str(nsamples)+" of them): "+
        #       str(self.diploid_samples)+"\n")
        # need chromosome ids
        chrom_samples = [ind_to_chrom(x, a)
                         for x in self.diploid_samples
                         for a in mapa_labels]
        # locations - repeated twice as it's for haploids
        sample_populations = [populations[k]
                              for k in sample_indices
                              for _ in range(2)]
        self.args.add_samples(samples=chrom_samples, length=self.length,
                              populations=sample_populations)
