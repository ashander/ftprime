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

def check_record_order(args):
    for ind in args:
        if len(args[ind][1])>0:
            x=0.0
            for cr in args[ind][1]:
                if x>cr.left or cr.left >= cr.right:
                    print("bad record order:",ind,args[ind])
                    raise ValueError
                x=cr.right

def test_simupop_runs():
    for popsize in [5,10,20]:
        print("Popsize: ",popsize)
        generations=3
        replications=1
        nsamples=2
        length=10
        nloci=5
        locus_position=list(range(0,length,int(length/nloci)))
        recomb_rate=0.05

        rc = RecombCollector(
                nsamples=nsamples, generations=generations, N=popsize, 
                ancestor_age=10, length=length, locus_position=locus_position)

        init_geno = [sim.InitGenotype(freq=[0.9,0.1],loci=sim.ALL_AVAIL)]

        id_tagger=sim.IdTagger(begin=0)
        id_tagger.reset(startID=1)  # must reset - creating a new one doesn't

        pop = sim.Population(
                size=[popsize],
                loci=[nloci], 
                lociPos=locus_position,
                infoFields=['ind_id'])

        pop.evolve(
            initOps=[
                sim.InitSex(),
                id_tagger
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

        rc.add_samples()

        check_record_order(rc.args)

        ts = rc.args.tree_sequence(samples=[(0,0) for x in range(rc.nsamples)])

        print("coalescence records:")
        for x in ts.records():
            print(x)

        ts.simplify(samples=list(range(nsamples)))

        print("trees:")
        for x in ts.trees():
            print(x)
