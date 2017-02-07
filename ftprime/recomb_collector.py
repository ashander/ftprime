import simuPOP as sim
from .argrecorder import ARGrecorder
from .meiosistagger import ind_to_chrom, mapa_labels
import math
import random

class RecombCollector:
    '''
    Collect and parse recombination events as output my simuPOP's Recombinator,
    which outputs like:
    offspringID parentID startingPloidy rec1 rec2 ....
    ... coming in *pairs*

    This needs:
    namples - number of *diploid* samples
    generations - number of generations simulation will be run for
    N - number of individuals in the population per generation
    ancestor_age - number of generations before beginning of simulation that common ancestor lived
    length - length of chromosome
    locus_position - list of positions of the loci simuPOP refers to along the chromosome
    '''

    def __init__(self, nsamples, generations, N, ancestor_age, length, locus_position):
        self.nsamples=nsamples
        self.generations=generations
        self.N=N
        self.ancestor_age=ancestor_age
        self.length=length
        self.locus_position=locus_position
        self.last_child=-1

        self.args=ARGrecorder()
        self.universal_ancestor=2*nsamples
        self.args.add_individual(name=self.universal_ancestor,time=float(1+self.generations+self.ancestor_age))
        # add initial generation
        first_gen = [self.i2c(k,p) for k in range(1,self.N+1) for p in [0,1]]
        first_gen.sort()
        self.args.add_record(0.0,self.length,self.universal_ancestor,tuple(first_gen))
        for k in range(1,self.N+1): 
            for p in [0,1]:
                self.args.add_individual(self.i2c(k,p),self.ind_to_time(k))

    def ind_to_time(self,k):
        # simuPOP has nonoverlapping gens so we can map indiv ID to time
        return 1+self.generations-math.floor((k-1)/self.N)

    def i2c(self,k,p):
        # individual ID to chromsome ID
        # "1+" is for the universal common ancestor added in initialization
        out=1+2*self.nsamples+ind_to_chrom(k,mapa_labels[p])
        return out
    
    def collect_recombs(self,lines):
        for line in lines.strip().split('\n'):
            # print("A: "+line)
            # child,parent,ploid,*rec = [int(x) for x in line.split()]
            linex = [int(x) for x in line.split()]
            child=linex[0]
            parent=linex[1]
            ploid=linex[2]
            rec=linex[3:]
            # lines come in pairs: maternal/paternal.
            if child==self.last_child:
                child_p=1
            else:
                child_p=0
                self.last_child=child
            # print(child,child_p,parent,ploid)
            if self.ind_to_time(child) > self.ind_to_time(parent):
                raise ValueError(str(child)+" at "+str(self.ind_to_time(child))+" does not come after " + str(parent)+" at "+str(self.ind_to_time(parent)))
            start=0.0
            child_chrom=self.i2c(child, child_p)
            # print(".. Adding",child_chrom)
            self.args.add_individual(child_chrom, self.ind_to_time(child))
            for r in rec:
                breakpoint=random.uniform(self.locus_position[r],self.locus_position[r+1])
                # print("--- ",start,breakpoint)
                self.args.add_record(
                        left=start,
                        right=breakpoint,
                        parent=self.i2c(parent,ploid),
                        children=(child_chrom,))
                start=breakpoint
                ploid=((ploid+1)%2)
            # print("--- ",start,self.length," |")
            self.args.add_record(
                    left=start,
                    right=self.length,
                    parent=self.i2c(parent,ploid),
                    children=(child_chrom,))

    def add_samples(self):
        # some messing around to fill up the required samples
        pop_ids = range(1+self.generations*self.N,1+(1+self.generations)*self.N)
        diploid_samples=random.sample(pop_ids,self.nsamples)
        # print("Samples ("+str(self.nsamples)+" of them): "+str(diploid_samples)+"\n")
        # need chromosome ids
        chrom_samples = [ ind_to_chrom(x,a) for x in diploid_samples for a in mapa_labels ]
        self.args.add_samples(samples=chrom_samples,length=self.length)

