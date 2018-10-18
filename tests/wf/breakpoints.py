import numpy as np


def random_breakpoint():
    return min(1.0, max(0.0, 2*np.random.uniform()-0.5))


def random_mutations(rate):
    nmuts = np.random.poisson(lam=rate)
    return np.random.uniform(size=nmuts)
