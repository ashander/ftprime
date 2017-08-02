from .wf import *
from .wf_with_call import *

import random

def random_breakpoint() :
    return min(1.0, max(0.0, 2*random.random()-0.5))

def random_mutations(rate) :
    nmuts = np.random.poisson(lam=rate)
    return [random.random() for _ in range(nmuts)]

