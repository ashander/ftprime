import simuPOP as sim
import random
from ftprime import MeiosisTagger, ind_to_chrom, mapa_labels
# from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch3_sec4.html
import simuOpt
simuOpt.setOptions(optimized=False, debug='DBG_WARNING')
# sim.turnOnDebug('DBG_ALL')
# sim.turnOnDebug('DBG_POPULATION,DBG_INDIVIDUAL')
sim.setOptions(seed=111)

def test_simupop_runs():
    popSize=5
    generations=3
    replications=1
    nsamples=5

    meioser=MeiosisTagger(nsamples,ngens=1+generations)
    def step_gen(pop,param):
        # function to increment internal clock in MeiosisTagger
        # as in http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec14.html#python-operator-pyoperator
        dt,=param
        print("Time was", meioser.gen," and now is ",meioser.gen+dt)
        meioser.gen+=dt
        return True

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
        preOps = sim.PyOperator(func=step_gen,param=(1,)),
        matingScheme = reproduction,
        gen = generations
    )

    # temporary reference to 0-th replicate population
    pop=simu.population(0)
    pop_ids = pop.indInfo("ind_id")

    print("final individuals:")
    for ind in pop.individuals():
        print("id: ", ind.info('ind_id'))

    locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
    meioser.add_diploid_samples(pop_ids,locations)

    for x in meioser.records.edgesets():
        print(x)

    for x in meioser.records.nodes():
        print(x)

    # print("Python trees:")
    # for t in trees(list(meioser.records.edgesets())):
    #     print(t)

    print("msprime trees:")
    # this is location,time for the samples
    ts=meioser.records.tree_sequence()

    for t in ts.trees():
        print(t)