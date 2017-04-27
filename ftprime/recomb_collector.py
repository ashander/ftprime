from .argrecorder import ARGrecorder
import random

# first is 'maternal', second is 'paternal'
mapa_labels = ( 2, 1 )  # NOQA


def ind_to_chrom(ind, mapa):
    '''
    Returns the unique *chromosome ID* corresponding to
    the chromosome of individual 'ind' inherited from
    parent 'mapa' (either 1 or 2).  (Chromosome IDs are ints.)
    '''
    return int(2*ind+mapa-1)


class RecombCollector:
    '''
    Collect and parse recombination events as output by simuPOP's Recombinator,
    which outputs like:
    offspringID parentID startingPloidy rec1 rec2 ....
    ... coming in *pairs*

    This needs:
    namples - number of *diploid* samples
    generations - number of generations simulation will be run for
    N - number of individuals in the population per generation
    ancestor_age - number of generations before beginning of simulation that common ancestor lived
    length - length of chromosome
    locus_position - list of positions of the loci simuPOP refers to along the chromosome

    This keeps track of *time* - so when used, the time must be updated -
    in simuPOP, by adding rc.increment_time() to the PreOps.  It should be in PreOps
    because if so
        - the initial generation is recorded at time generations + 1
        - calling increment_time() decrements the time by 1
        - the first generation is recorded at time generations
        - ...
        - the final generation is at time 1
        - the samples are at time 0
    '''

    def __init__(self, generations, N, ancestor_age, length,
                 locus_position):
        self.generations = generations
        self.N = N
        self.ancestor_age = ancestor_age
        self.length = length
        self.locus_position = locus_position
        self.last_child = -1
        self.time = float(generations) + 1
        self.nsamples = 0

        self.args = ARGrecorder()

        self.universal_ancestor = 0
        # will record IDs of diploid samples here when they are chosen
        # but note we don't keep anything else about them here (time, location)
        # as this is recorded by the ARGrecorder
        self.diploid_samples = None
        max_time = float(1 + self.generations + self.ancestor_age)
        self.args.add_individual(name=self.universal_ancestor,
                                 time=max_time)
        # add initial generation
        # uses the fact that simupop will number inds in first gen like 1..n
        # TODO - make sure this matches up with ind ids in first gen
        first_gen = [self.i2c(k, p) for k in range(1, self.N+1) for p in [0, 1]]
        first_gen.sort()
        self.args.add_record(
                left=0.0,
                right=self.length,
                parent=self.universal_ancestor,
                children=tuple(first_gen))
        for k in range(1, self.N + 1):
            for p in [0, 1]:
                # print("Adding:",k, p,self.i2c(k, p), self.time)
                self.args.add_individual(self.i2c(k, p), self.time)

    def i2c(self, k, p):
        # individual ID to chromsome ID
        # "1+" is for the universal common ancestor added in initialization
        out = 1 + ind_to_chrom(k, mapa_labels[p])
        return out

    def increment_time(self):
        self.time -= 1.0

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
        self.nsamples = nsamples
        # TODO: remove previous samples
        assert(len(sample_ids) == len(populations))
        sample_indices = random.sample(range(len(sample_ids)), nsamples)
        self.diploid_samples = [sample_ids[k] for k in sample_indices]
        # print("Samples ("+str(self.nsamples)+" of them): "+
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
