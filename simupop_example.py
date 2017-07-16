import simuPOP as sim
import math
import random
from ftprime import RecombCollector

popsize = 100
nsamples = 2
length = 1000
nloci = 11
recomb_rate = 0.001

# 1. Must have a locus at each end of the chromosome
locus_position = [length * x/(nloci-1) for x in range(nloci)]
                                                                                 
init_geno = [sim.InitGenotype(freq=[0.9, 0.1], loci=sim.ALL_AVAIL)]

# 2. Must include 'ind_id' infoField for sim.Recombinator
pop = sim.Population(
        size=[popsize],
        loci=[nloci],
        lociPos=locus_position,
        infoFields=['ind_id'])

# 3. Must tag the first generation so we can use it when initializing RecombCollector
id_tagger = sim.IdTagger()
# note: if doing this more than once in the same python session must do this to reset:
#   creating a new IdTagger doesn't reset it
id_tagger.reset(startID=1)

# tag the first generation so we can pass it to rc
id_tagger.apply(pop)

# 4. Initialize, including how long until a common ancestor before the present
rc = RecombCollector(first_gen=pop.indInfo("ind_id"), ancestor_age=10, 
                     length=length, locus_position=locus_position)

recombinator = sim.Recombinator(intensity=recomb_rate,
                                output=rc.collect_recombs,
                                infoFields="ind_id")

pop.evolve(
    initOps=[
        sim.InitSex()
    ]+init_geno,
    preOps=[
        # 5. Must keep time up to date in the RecombCollector
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        # Must return true or false. True keeps whole population (?)
    ],
    matingScheme=sim.RandomMating(
        ops=[id_tagger, recombinator]),
    postOps=[
        sim.PyEval(r"'Gen: %2d\n' % (gen, )", step=1)
    ],
    gen=20
)

# 6. After, must specify which individuals are sampled.
#    This needs to know the locations of the individuals to be able to record this.
locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_diploid_samples(nsamples, pop.indInfo("ind_id"), locations)

# 7. Output the tree sequence from the recorded information.
full_ts = rc.args.tree_sequence()

# 8. Optionally, simplify it to remove redundant information.
ts = full_ts.simplify()

# 9. Write out the tree sequence in a file.
ts.dump("simulated.ts")
