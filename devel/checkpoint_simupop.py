#!/usr/bin/env python3

import gzip
import sys, os
import math
import time
import random
from checkpointed_recomb_collector import CheckpointedRecombCollector
import msprime
import argparse

description='For testing.'

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--level", "-l", type=int, help="level of computation to stop at.")
parser.add_argument("--generations", "-T", type=int, help="number of generations to run for.")
parser.add_argument("--popsize", "-N", type=int, help="size of each subpopulation")
parser.add_argument("--length", "-L", type=float, help="number of bp in the chromosome")
parser.add_argument("--recomb_rate", "-r", type=float, help="recombination rate")
parser.add_argument("--nsamples", "-k", type=int, help="number of *diploid* samples, total")
parser.add_argument("--ancestor_age", "-A", type=float, help="time to ancestor above beginning of sim")
parser.add_argument("--mut_rate", "-U", type=float, help="mutation rate")
parser.add_argument("--seed", "-d", type=int, help="random seed", default=random.randrange(1,1000))
parser.add_argument("--outdir", "-o", help="name of output directory", default=None)

args = parser.parse_args()

import simuPOP as sim

sim.setRNG(seed=args.seed)
random.seed(args.seed)

if args.outdir is None:
    args.outdir = "sim_" + str(args.seed)
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
args.treefile = args.outdir + "/sim.trees"
args.vcffile = args.outdir + "/sim.vcf"
args.logfile = args.outdir + "/sim.log"
args.samples_file = args.outdir + "/samples.tsv"

vcffile = open(args.vcffile, "w")
logfile = open(args.logfile, "w")
samples_file = open(args.samples_file,"w")

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

pop = sim.Population(
        size=[args.popsize],
        loci=[2], 
        lociPos=[0, args.length],
        infoFields=['ind_id','fitness'])

id_tagger = sim.IdTagger()
id_tagger.apply(pop)

# record recombinations
rc = CheckpointedRecombCollector(
        first_gen=pop.indInfo("ind_id"), 
        ancestor_age=args.ancestor_age, 
        length=args.length, 
        locus_position=[0,args.length],
        level=args.level)

pop.evolve(
    initOps=[
        sim.InitSex(),
    ],
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            id_tagger,
            sim.Recombinator(intensity=args.recomb_rate,
                output=rc.collect_recombs,
                infoFields="ind_id"),
        ] ),
    postOps=[
        sim.PyEval(r"'Gen: %2d\n' % gen ", step=50)
    ],
    gen = args.generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

locations = [pop.subPopIndPair(x)[0] for x in range(pop.popSize())]
rc.add_diploid_samples(nsamples=args.nsamples, 
                       sample_ids=pop.indInfo("ind_id"),
                       populations=locations)

del pop

logfile.write("Samples:\n")
logfile.write(str(rc.diploid_samples)+"\n")
logfile.write("----------\n")
logfile.flush()

logfile.write("All done!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.close()

sys.exit()


## DON'T DO THIS STUFF

rc.args.dump_sample_table(out=samples_file)

ts = rc.args.tree_sequence()
del rc

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

# nodefile = open("nodes.txt","w")
# edgefile = open("edges.txt","w")
# ts.dump_text(nodes=nodefile, edges=edgefile)
# nodefile.flush()
# edgefile.flush()

minimal_ts = ts.simplify()
del ts

logfile.write("Simplified; now writing to treefile.\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

minimal_ts.dump(args.treefile)

mut_seed=args.seed
logfile.write("Done, now generating mutations with seed "+str(mut_seed)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

rng = msprime.RandomGenerator(mut_seed)
nodes = msprime.NodeTable()
edges = msprime.EdgeTable()
sites = msprime.SiteTable()
mutations = msprime.MutationTable()
minimal_ts.dump_tables(nodes=nodes, edges=edges)
mutgen = msprime.MutationGenerator(rng, args.mut_rate)
mutgen.generate(nodes, edges, sites, mutations)

# print(nodes, file=logfile)
# print(edges, file=logfile)
# print(sites, file=logfile)
# print(mutations, file=logfile)

mutated_ts = msprime.load_tables(
    nodes=nodes, edges=edges, sites=sites, mutations=mutations)

del minimal_ts

logfile.write("Generated mutations!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("Mean pairwise diversity: {}\n".format(mutated_ts.get_pairwise_diversity()/mutated_ts.get_sequence_length()))
logfile.write("Sequence length: {}\n".format(mutated_ts.get_sequence_length()))
logfile.write("Number of trees: {}\n".format(mutated_ts.get_num_trees()))
logfile.write("Number of mutations: {}\n".format(mutated_ts.get_num_mutations()))

mutated_ts.write_vcf(vcffile,ploidy=1)

logfile.write("All done!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.close()
