import simuPOP as sim
import numpy.random as np_random
from itertools import count
from pedrecorder import PedigreeRecorder

# *_chromo fields are in the dict of values;
#  ... this MUST be synchronized with the set that `f_chromo` and `m_chromo` are chosen from below:
labels = { 'paternal' : 1, 'maternal' : 2 }

def random_breakpoint() :
    return min(1.0,max(0.0, 2*np_random.random()-0.5)) 

def make_labeler(nsamples):
    next_id=float(nsamples)
    while True:
        yield next_id
        next_id += 1.0


class MeiosisTagger(sim.PyOperator):
    def __init__(self, nsamples, *args, **kwargs):
        self.nsamples = nsamples
        self.labeler=make_labeler(nsamples)
        self.records=PedigreeRecorder()
        self.time={}
        sim.PyOperator.__init__(self, func=self.init_labels, *args, **kwargs)

    def init_labels(self,pop):
        '''
        This initializes the labels.
        '''
        for ind in pop.individuals():
            next_id = next(self.labeler)
            self.records.add_individual(next_id)
            self.time[next_id]=99
            ind.setInfo(next_id,'ind_id')
        return True

    def new_offspring(self,ind_id):
        '''
        This gives new indivdiuals a unique ID
        but also along the way records parental IDs and recombination locations
        in records.
        '''
        print("new_off:",self.records)
        child=next(self.labeler)
        bp=random_breakpoint()
        if np_random.random() < 0.5:
            lparent,rparent=ind_id
        else:
            rparent,lparent=ind_id
        print('meiosis:',child,lparent,rparent,bp)
        ## Figure out how to get time in here from simuPOP
        self.time[child]=99
        self.records.add_individual(child)
        if bp > 0.0 :
            self.records.add_record(left=0.0, right=bp, parent=lparent, children=(child,), time=self.time[lparent], population=0)
        if bp < 1.0 :
            self.records.add_record( left=bp, right=1.0, parent=rparent, children=(child,), time=self.time[rparent], population=0)
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

