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
                    print(ind,args[ind])
                    raise ValueError
                x=cr.right

def test_simupop_runs():
    popsize=5
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

    pop = sim.Population(
            size=[popsize],
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


# Troublesome:
# 
# 14.0     2___    .     2___    . 
#         / \  \   .    / \  \   . 
# 4.0    8   7  15 .   8   7  15 . 
#        |   |   | .   |   |   | . 
# 3.0   22   18 17 .  22   18 17 . 
#        |    |  | .   |    |  | . 
# 2.0   27   28 30 .  27   28 30 . 
#       /\       | .  /\       | . 
# 1.0  39 41     | . 39 41     | . 
#       |        | .  |        | . 
# 0.0   0        1 .  0        1 . 
#                  .             . 
#     0.0 --------0.4 ---------- 
#  
# 
# 
# records=[
#     cr(left=0.0, right=10.0, node=39, children=(0,), time=1.0, population=0),
#     cr(left=0.0, right=10.0, node=30, children=(1,), time=2.0, population=0),
#     cr(left=0.0, right=0.4, node=27, children=(39, 41), time=2.0, population=0),
#     cr(left=0.4, right=0.5, node=27, children=(41,), time=2.0, population=0),
#     cr(left=0.4, right=0.7, node=27, children=(39,), time=2.0, population=0),
#     cr(left=0.5, right=0.7, node=28, children=(41,), time=2.0, population=0),
#     cr(left=0.7, right=4.5, node=28, children=(39, 41), time=2.0, population=0),
#     cr(left=4.5, right=10.0, node=28, children=(39, 41), time=2.0, population=0),
#     cr(left=0.0, right=10.0, node=22, children=(27,), time=3.0, population=0),
#     cr(left=0.0, right=10.0, node=17, children=(30,), time=3.0, population=0),
#     cr(left=0.0, right=10.0, node=18, children=(28,), time=3.0, population=0),
#     cr(left=0.0, right=10.0, node=15, children=(17,), time=4.0, population=0),
#     cr(left=0.0, right=4.5, node=7, children=(18,), time=4.0, population=0),
#     cr(left=4.5, right=10.0, node=7, children=(18, 22,), time=4.0, population=0),
#     cr(left=0.0, right=4.5, node=8, children=(22,), time=4.0, population=0),
#     cr(left=0.0, right=10.0, node=2, children=(7, 8, 15), time=14.0, population=0)
#     ]
# 
# 
# _ts = _msprime.TreeSequence()
# _ts.load_records(records=records)
# ts = msprime.TreeSequence(_ts)
# 
# ts.simplify(samples=[0,1])
# 
