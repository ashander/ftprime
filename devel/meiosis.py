import simuPOP as sim
import numpy.random as np_random
from itertools import count
from pedrecorder import PedigreeRecorder
from math import floor

# first is 'maternal', second is 'paternal'
mapa_labels = ( 2, 1 )

def random_breakpoint() :
    return min(1.0,max(0.0, 2*np_random.random()-0.5)) 

def make_labeler(nsamples):
    # simuPOP's infoFields are coerced to floats
    next_id=float(nsamples)
    while True:
        yield next_id
        next_id += 1.0

def ind_to_chrom(ind,mapa):
    # mapa is either 1 or 2 (see above)
    # and chrom IDs are ints...
    return int(2*ind+mapa-1)

def chrom_to_ind(chrom):
    # but ind IDs are floats (simuPOP)
    return floor(float(chrom)/2)

class MeiosisTagger(sim.PyOperator):
    def __init__(self, nsamples, ngens, *args, **kwargs):
        self.nsamples = nsamples  # number of INDIVIDUALS
        self.gen=0   # the current time (in forwards-time)
        self.ngens = ngens # maximum number of generations
        self.labeler=make_labeler(nsamples)  # returns INDIVIDUAL labels
        self.records=PedigreeRecorder()  # indexed by CHROMOSOMES
        self.time={}  # indexed by INDIVIDUALS
        sim.PyOperator.__init__(self, func=self.init_labels, *args, **kwargs)

    def init_labels(self,pop):
        '''
        This initializes the labels.
        '''
        for ind in pop.individuals():
            next_id = next(self.labeler)
            for mapa in mapa_labels: # mapa_labels.values():
                self.records.add_individual(ind_to_chrom(next_id,mapa))
            self.time[next_id]=self.ngens
            ind.setInfo(next_id,'ind_id')
        return True

    def new_offspring(self,ind_id):
        '''
        This gives new indivdiuals a unique ID
        but also along the way records parental IDs and recombination locations
        in records.
        Note that self.gen must be kept up to date externally!
        '''
        # print("new_off:",self.records)
        child=next(self.labeler)
        ## Figure out how to get time in here from simuPOP
        self.time[child]=self.ngens-self.gen
        # which parent and which chrom in offspring
        for parent,mapa in zip(ind_id,mapa_labels): #mapa_labels.values()):
            bp=random_breakpoint()
            # order of inheritance from parent
            lparent,rparent=np_random.permutation(mapa_labels) #mapa_labels.values())
            # print('meiosis:',child,parent,mapa,lparent,rparent,bp)
            chrom_id=ind_to_chrom(child,mapa)
            self.records.add_individual(chrom_id)
            if bp > 0.0 :
                self.records.add_record(left=0.0, right=bp, parent=ind_to_chrom(parent,lparent), children=(chrom_id,), time=self.time[parent], population=0)
            if bp < 1.0 :
                self.records.add_record( left=bp, right=1.0, parent=ind_to_chrom(parent,rparent), children=(chrom_id,), time=self.time[parent], population=0)
        # return label for new individual
        return child
