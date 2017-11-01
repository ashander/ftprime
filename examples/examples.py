#!/usr/bin/env python3.5
description = '''
Simulate a chromosome with a single selected loci using simuPOP.
'''

import gzip
import sys
from argparse import ArgumentParser
import math
import time
import random
from ftprime import RecombCollector
import msprime

REPORTING_STEP = 50

parser = ArgumentParser(description=description)
parser.add_argument("-T","--generations", dest="generations", type=int,
        help="number of generations to run for")
parser.add_argument("-N","--popsize", dest="popsize", type=int,
        help="size of the population", default=100)
parser.add_argument("-r","--recomb_rate", dest="recomb_rate", type=float,
        help="recombination rate", default=1e-7)
parser.add_argument("-L","--length", dest="chrom_length", type=int,
        help="number of bp in the chromosome", default=100)
parser.add_argument("-U","--neut_mut_rate", dest="neut_mut_rate", type=float,
        help="neutral mutation rate", default=1e-7)
parser.add_argument("-l","--nselloci", dest="nselloci", type=int,
        help="number of selected loci", default=1)
parser.add_argument("-u","--sel_mut_rate", dest="sel_mut_rate", type=float,
        help="mutation rate of selected alleles", default=1e-7)
parser.add_argument("-a","--gamma_alpha", dest="gamma_alpha", type=float,
        help="alpha parameter in gamma distributed selection coefficient", default=.23)
parser.add_argument("-b","--gamma_beta", dest="gamma_beta", type=float, 
        help="beta parameter in gamma distributed selection coefficient", default=5.34)
parser.add_argument("-k","--nsamples", dest="nsamples", type=int,
        help="number of *diploid* samples, total", default=100)
parser.add_argument("-o","--outfile", dest="outfile", type=str,
        help="name of output PED file (default: not output)", default=None)
parser.add_argument("--gc", "-G", dest="simplify_interval", type=int,
        help="Interval between simplify steps.", default=500)
parser.add_argument("-g","--logfile", dest="logfile", type=str,
        help="name of log file (or '-' for stdout)", default="-")
parser.add_argument("-s","--selloci_file", dest="selloci_file", type=str,
        help="name of file to output selected locus information", default="sel_loci.txt")
parser.add_argument("--treefile","-t", type=str, dest="treefile",
        help="name of output file for trees (default: not output)",default=None)

args = parser.parse_args()

if args.generations is None:
    parser.print_help()
    sys.exit()

# some simupop options involving mutation type
import simuOpt
simuOpt.setOptions(alleleType='mutant')
import simuPOP as sim

def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print("Something not right here.")
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj

logfile = fileopt(args.logfile, "w")
selloci_file = args.selloci_file

logfile.write("Options:\n")
logfile.write(str(args)+"\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

# locations of the loci along the chromosome?
# hard code defaults for simupop:
# >The default positions are 1, 2, 3, 4, ... on each
# >chromosome.
locus_position = list(range(0, args.chrom_length))

# which loci are under selection?
selected_loci = math.ceil(args.chrom_length / 2)
try:
    sl  = set(selected_loci)
except TypeError:
    sl = set((selected_loci, ))
neutral_loci = list(set(range(1,args.chrom_length)) - sl)

###
# random selection coefficients:
# modified from http://simupop.sourceforge.net/manual_svn/build/userGuide_ch5_sec9.html

class GammaDistributedFitness:
    def __init__(self, alpha, beta):
        # mean is alpha/beta
        self.coefMap = {}
        self.alpha = alpha
        self.beta = beta
    def __call__(self, loc, alleles):
        # because s is assigned for each locus, we need to make sure the
        # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s = self.coefMap[loc]
        else:
            s = random.gammavariate(self.alpha, self.beta)
            self.coefMap[loc] = s
        # print(str(loc)+":"+str(alleles)+"\n")
        # needn't return fitness for alleles=(0,0) as simupop knows that's 1
        if 0 in alleles:
            return 1. - s
        else:
            return 1. - 2.*s

pop = sim.Population(
        size=args.popsize,
        loci=[args.chrom_length],
        lociPos=locus_position,
        infoFields=['ind_id','fitness'])


# set up recomb collector
# NB: we have to simulate an initial tree sequence
id_tagger = sim.IdTagger()
id_tagger.apply(pop)
first_gen = pop.indInfo("ind_id")
length = max(locus_position)

# Since we want to have a finite site model, we force the recombination map
# to have exactly `length` loci with a fixed recombination rate between them.
rcmb_map = msprime.RecombinationMap.uniform_map(length, args.recomb_rate, length)
init_ts = msprime.simulate(2*len(first_gen), recombination_map=rcmb_map)
print('\n Init ts')
print(init_ts.dump_tables())
haploid_labels = [(k,p) for k in first_gen
                        for p in (0,1)]
node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}
rc = RecombCollector(ts=init_ts, node_ids=node_ids,
                     locus_position=locus_position)
print('\n In recomb collector')
print(rc.args.tree_sequence().dump_tables())
sys.exit(0)
# initially, population is monogenic
init_geno=[sim.InitGenotype(freq=1.0)]

pop.evolve(
    initOps=[
        sim.InitSex(),
    ]+init_geno,
    preOps=[
        sim.PyOperator(lambda pop: rc.increment_time() or True),
        sim.SNPMutator(u=args.neut_mut_rate,v=0,loci=neutral_loci),
        sim.SNPMutator(u=args.sel_mut_rate,v=0,loci=selected_loci),
        sim.PyMlSelector(GammaDistributedFitness(args.gamma_alpha, args.gamma_beta),
            loci=selected_loci, output=">>"+selloci_file),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            id_tagger,
            sim.Recombinator(rates=args.recomb_rate, output=rc.collect_recombs,
                             infoFields="ind_id"),
        ] ),
    postOps=[
        sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=REPORTING_STEP,
                 vars=['numOfSegSites', 'numOfFixedSites']),
        sim.PyEval(r"'Gen: %2d #seg/#fixed sites: %d / %d\n' % (gen, numOfSegSites, numOfFixedSites)", step=REPORTING_STEP),
        sim.PyOperator(lambda pop: rc.simplify(pop.indInfo("ind_id")) or True,
                       step=args.simplify_interval),
    ],
    gen = args.generations
)

logfile.write("Done simulating!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

logfile.write("Collecting samples:\n")
logfile.write("  " + str(args.nsamples) + " of them")
# logfile.write("  " + "ids:" + str(pop.indInfo("ind_id")))

diploid_samples = random.sample(pop.indInfo("ind_id"), args.nsamples)
rc.simplify(diploid_samples)

del pop

logfile.write("Samples:\n")
#logfile.write(str(rc.diploid_samples)+"\n")
logfile.write("----------\n")
logfile.flush()

ts = rc.args.tree_sequence()
del rc

logfile.write("Loaded into tree sequence!\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

if args.treefile is not None:
    ts.dump(args.treefile)

logfile.write("Writing out samples.\n")
logfile.write(time.strftime('%X %x %Z')+"\n")
logfile.write("----------\n")
logfile.flush()

logfile.write("All done!\n")
logfile.close()
