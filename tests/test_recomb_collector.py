import simuPOP as sim
import random
from ftprime import RecombCollector, ind_to_chrom, mapa_labels
import math
# from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch3_sec4.html
import simuOpt
simuOpt.setOptions(optimized=False, debug='DBG_WARNING')
# sim.turnOnDebug('DBG_ALL')
# sim.turnOnDebug('DBG_POPULATION,DBG_INDIVIDUAL')
sim.setOptions(seed=111)

def test_simupop_runs():
    popsize=5
    generations=3
    replications=1
    nsamples=5
    length=100
    nloci=10
    locus_position=list(range(0,length,length/nloci))
    recomb_rate=1.0

    rc = RecombCollector(
            nsamples=nsamples, generations=generations, N=popsize, 
            ancestor_age=10, length=length, locus_position=locus_position)

    init_geno = sim.InitGenotype(freq=0.1,loci=list(range(nloci)))

    pop = sim.Population(
            size=[popsize]*npops, 
            loci=[nloci], 
            lociPos=locus_position,
            infoFields=['ind_id'])

    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.IdTagger(),
        ]+init_geno,
        preOps=[],
        matingScheme=sim.RandomMating(
            ops=[
                sim.IdTagger(),
                sim.Recombinator(intensity=recomb_rate,
                    output=rc.collect_recombs,
                    infoFields="ind_id"),
            ] ),
        postOps=[
            sim.PyEval(r"'Gen: %2d\n' % (gen, )", step=1)
        ],
        gen = generations
    )

