import simuPOP as sim
import math
import random
from ftprime import RecombCollector
import msprime

popsize = 10
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


# 4. Simulate population history before the start of the simulation,
#    and which individuals in this history correspond to indvidiuals
#    in the simuPOP simulation.
first_gen = pop.indInfo("ind_id")
init_ts = msprime.simulate(2*len(first_gen), # this is haploid simulation
                           length=max(locus_position))
haploid_labels = [(k,p) for k in first_gen 
                        for p in (0,1)]
node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}

# 5. Initialize
rc = RecombCollector(ts=init_ts, node_ids=node_ids,
                     locus_position=locus_position)

recombinator = sim.Recombinator(intensity=recomb_rate,
                                output=rc.collect_recombs,
                                infoFields="ind_id")


# 6. Run the simulation, simplifying the underlying tree sequence
#    every 5 generations (helps reduce memory usage for bigger sims)
simplify_interval = 5

pop.evolve(
    initOps=[
        sim.InitSex()
    ]+init_geno,
    preOps=[
        # 5. Must keep time up to date in the RecombCollector
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        # Must return true or false. True keeps whole population (?)
        sim.PyOperator(lambda pop: rc.simplify(pop.indInfo("ind_id")) or True, step=simplify_interval),
    ],
    matingScheme=sim.RandomMating(
        ops=[id_tagger, recombinator]),
    postOps=[
        sim.PyEval(r"'Gen: %2d\n' % (gen, )", step=1),
    ],
    gen=20
)


# 7. Do one last simplify().
rc.simplify(pop.indInfo("ind_id"))

# 8. Copy location information over from simuPOP for the final generation.
locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_locations(pop.indInfo("ind_id"), locations)

# 9. Choose which individuals are sampled.
rc.add_diploid_samples(nsamples, pop.indInfo("ind_id"))

# 10. Output the tree sequence from the recorded information.
ts = rc.args.tree_sequence()

# 11. Write out the tree sequence in a file.
ts.dump("simulated.ts")
