import simuPOP as sim
import numpy.random as np_random
from itertools import count
from .argrecorder import ARGrecorder
from math import floor

# first is 'maternal', second is 'paternal'
mapa_labels = ( 2, 1 )

def random_breakpoint() :
    return min(1.0,max(0.0, 2*np_random.random()-0.5))

def make_labeler(nsamples):
    '''
    Returns sequential integers starting from nsamples,
    suitable for assigning unique IDs to individuals.
    (These are floats, as that is what simuPOP uses for IDs.)
    '''
    next_id=float(nsamples)
    while True:
        yield next_id
        next_id += 1.0

def ind_to_chrom(ind,mapa):
    '''
    Returns the unique *chromosome ID* corresponding to
    the chromosome of individual 'ind' inherited from 
    parent 'mapa' (either 1 or 2).  (Chromosome IDs are ints.)
    '''
    return int(2*ind+mapa-1)

def chrom_to_ind(chrom):
    '''
    Provides the 2-to-1 inverse map from chromosome ID back to
    the corresponding individual ID (a float, as above).
    '''
    return floor(float(chrom)/2)

class MeiosisTagger(sim.PyOperator):
    '''
    This class provides an interface between simuPOP and the ARGrecorder,
    which requires we keep track of:

        - `gen` : keeps track of the current time (must be updated)

        - `ngens` : total number of generations that we will simulate for,
            so that at the end we can convert to reverse-time; this must be
            strictly larger than the birth generation of any sampled individual.

        - `labeler` : used to output unique individual IDs

        - `time` : a dictionary of birth times of individuals,
            which only needs to be currently alive individuals
            but we don't bother to remove dead ones from it
        
    All the work is done in `new_offspring`; from simuPOP's point of view
    all that does is return a new individual ID; but internally it also
    generates recombination events and stores them in an ARGrecorder.
    '''
    
    def __init__(self, nsamples, ngens, *args, **kwargs):
        self.nsamples = nsamples  # number of INDIVIDUALS
        self.gen=0   # the current time (in forwards-time)
        self.ngens = ngens # maximum number of generations
        self.labeler=make_labeler(nsamples)  # returns INDIVIDUAL labels
        self.records=ARGrecorder()  # indexed by CHROMOSOMES
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
