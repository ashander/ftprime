import simuPOP as sim
from meiosis import MeiosisTagger

def random_breakpoint() :
    return min(1.0,max(0.0, 2*np_random.random()-0.5)) 

popSize=20
generations=10
replications=1
nsamples=5

def make_labeler(nsamples):
    next_id=float(nsamples)
    while True:
        yield next_id
        next_id += 1.0


meioser=MeiosisTagger(nsamples)


reproduction = sim.RandomMating(
    ops=[
        # sim.IdTagger(),                   # new ID for offspring
        sim.PyTagger(func=meioser.new_offspring),
        # PyTagger operates on infofields defined by names of the arguments of func
    ]
)

pop = sim.Population(size=popSize, infoFields=['ind_id'])
simu=sim.Simulator(pop)

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

for ind in pop.individuals():
    print("id: ", ind.info('ind_id'))
