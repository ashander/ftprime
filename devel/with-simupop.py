import simuPOP as sim
import msprime
import random
from meiosis import MeiosisTagger
# from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch3_sec4.html
import simuOpt
simuOpt.setOptions(optimized=False, debug='DBG_WARNING')
sim.turnOnDebug('DBG_ALL')
# sim.turnOnDebug('DBG_POPULATION,DBG_INDIVIDUAL')

popSize=20
generations=10
replications=1
nsamples=5

meioser=MeiosisTagger(nsamples)


reproduction = sim.RandomMating(
    ops=[
        # sim.IdTagger(),                   # new ID for offspring
        sim.PyTagger(func=meioser.new_offspring),
        # PyTagger operates on infofields defined by names of the arguments of func
    ]
)

# note that 'pop' is a temporary object.
pop = sim.Population(size=popSize, infoFields=['ind_id'])
simu = sim.Simulator(pop)

simu.evolve(
    # everyone initially will have the same allele frequency
    initOps = [
        sim.InitSex(),
        # sim.IdTagger(),
        # sim.PyOperator(func=init_labels),
        meioser
    ],
    matingScheme = reproduction,
    gen = generations
)

# temporary reference to 0-th replicate population
pop=simu.population(0)
pop_ids = [ ind.info('ind_id') for ind in pop.individuals() ]

print("final individuals:")
for ind in pop.individuals():
    print("id: ", ind.info('ind_id'))

for x in meioser.records.dump_records():
    print(x)

samples=random.sample(pop_ids,nsamples)
times=[meioser.time[x] for x in samples]
meioser.records.add_samples(samples=samples,times=times,populations=[0 for x in samples])

ts=meioser.records.tree_sequence()

for t in ts.trees():
    print(t)
