import simuPOP as sim
import numpy.random as np_random
from itertools import count
from pedrecorder import PedigreeRecorder
from math import floor

# first is 'maternal', second is 'paternal'
mapa_labels = [ 2, 1 ]

def random_breakpoint() :
    return min(1.0,max(0.0, 2*np_random.random()-0.5)) 

def make_labeler(nsamples):
    # simuPOP's infoFields are coerced to floats
    next_id=float(nsamples)
    while True:
        yield next_id
        next_id += 1.0

def ind_to_chrom(ind,mapa):
    # mapa is either 1 or 2
    return 2*ind+mapa-1

def chrom_to_ind(chrom):
    return floor(chrom/2)

class MeiosisTagger(sim.PyOperator):
    def __init__(self, nsamples, *args, **kwargs):
        self.nsamples = nsamples  # number of INDIVIDUALS
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
            self.time[next_id]=99
            ind.setInfo(next_id,'ind_id')
        return True

    def new_offspring(self,ind_id):
        '''
        This gives new indivdiuals a unique ID
        but also along the way records parental IDs and recombination locations
        in records.
        '''
        # print("new_off:",self.records)
        child=next(self.labeler)
        ## Figure out how to get time in here from simuPOP
        self.time[child]=99
        # which parent and which chrom in offspring
        for parent,mapa in zip(ind_id,mapa_labels): #mapa_labels.values()):
            bp=random_breakpoint()
            # order of inheritance from parent
            lparent,rparent=np_random.permutation(mapa_labels) #mapa_labels.values())
            print('meiosis:',child,parent,mapa,lparent,rparent,bp)
            chrom_id=ind_to_chrom(child,mapa)
            self.records.add_individual(chrom_id)
            if bp > 0.0 :
                self.records.add_record(left=0.0, right=bp, parent=ind_to_chrom(parent,lparent), children=(chrom_id,), time=self.time[parent], population=0)
            if bp < 1.0 :
                self.records.add_record( left=bp, right=1.0, parent=ind_to_chrom(parent,rparent), children=(chrom_id,), time=self.time[parent], population=0)
        # return label for new individual
        return child

def meiosis(ind_id):
    """choose a starting chromo and breakpoint from each parent

    father_chromo -- tuple of parents' father_chromos
    father_breakpoint -- tuple of parents' father_breakpoints
    mother_chromo -- tuple of parents' mother_chromos
    mother_breakpoint -- tuple of parents' mother_breakpoints

    The arguments are unused. Meiosis does two things: (1) Generate a random
    starting chromosome from both father and mother (father_chromo,
    mother_chromo) with 1 = paternal and 2 = maternal in each case,
    which agrees with simuPOPs convention for .sex()
    (see above definition of `labels`)

    (2) generate a breakpoint for each parent, which is a binomial mixture of
    50% 0 (for no recombination) and 50 % (0, 1] (for recombination at some
    point along the chromosome)."""

    print("ind_id:", ind_id)

    # these two specify which parent the left end of the chromosome comes from:
    # possible values are 1,2 which must be synchronized with the above `labels`
    f_chromo, m_chromo = np_random.random_integers(1, 2, 2)

    # these two random integers in [0, 1] are used as indices for the
    # mixture to construct the breakpoint
    f_idx, m_idx = np_random.random_integers(0, 1, 2)

    m_bp_if_recombine, f_bp_if_recombine = np_random.ranf(2)
    # (0, 1] (np.random.ranf is [0, 1))
    f_breakpoint = [1.0, f_bp_if_recombine][f_idx]
    m_breakpoint = [1.0, m_bp_if_recombine][m_idx]

    # np.random.choice chooses from the array passed with uniform probability
    # return f_chromo, f_breakpoint, m_chromo, m_breakpoint
    return ind_id[0]

