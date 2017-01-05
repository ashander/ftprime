
import unittest
import itertools
import random

import msprime
import _msprime
from branch_stats import *

##
# This tests implementation of the algorithm described in branch-lengths-methods.md
##

def build_tree_sequence(records, mutations=[]):
    ts = _msprime.TreeSequence()
    ts.load_records(records)
    ts.set_mutations(mutations)
    return msprime.TreeSequence(ts)


def path_length(tr,x,y):
    L = 0
    mrca = tr.mrca(x,y)
    for u in x,y:
        while u != mrca:
            L += tr.branch_length(u)
            u = tr.parent(u)
    return L

def branch_length_diversity(ts,x,y):
    S = 0
    for tr in ts.trees():
        S += path_length(tr,x,y)*tr.length
    return S/ts.sequence_length

def branch_length_Y(ts,x,y,z):
    S = 0
    for tr in ts.trees():
        xy_mrca = tr.mrca(x,y)
        xz_mrca = tr.mrca(x,z)
        yz_mrca = tr.mrca(y,z)
        if xy_mrca == xz_mrca:
            #   /\
            #  / /\        
            # x y  z       
            S += path_length(tr,x,yz_mrca)*tr.length
        elif xy_mrca == yz_mrca:
            #   /\
            #  / /\        
            # y x  z       
            S += path_length(tr,x,xz_mrca)*tr.length
        elif xz_mrca == yz_mrca:
            #   /\
            #  / /\        
            # z x  y       
            S += path_length(tr,x,xy_mrca)*tr.length
    return S/ts.sequence_length

class BranchStatsTestCase(unittest.TestCase):
    """
    Tests of branch statistic computation.
    """
    random_seed = 123456

    def check_pairwise_diversity(self,ts):
        samples = random.sample(ts.samples(),2)
        A = [ [samples[0]], [samples[1]] ]
        def f(x):
            return (x[0]>0) != (x[1]>0)
        self.assertEqual(
                branch_stats_node_iter(ts,A,f,method='mutations'),
                ts.pairwise_diversity(samples=samples) )
        self.assertAlmostEqual(
                branch_stats_node_iter(ts,A,f,method='length'),
                branch_length_diversity(ts,samples[0],samples[1]) )
        self.assertEqual(
                branch_stats(ts,A,f),
                branch_length_diversity(ts,samples[0],samples[1]) )

    def check_Y_stat(self,ts):
        samples = random.sample(ts.samples(),3)
        A = [ [samples[0]], samples[1:3] ]
        def f(x):
            return ( (x[0]==1) and (x[1]==0) ) or ( (x[0]==0) and (x[1]==2) )
        self.assertEqual(
                branch_stats(ts,A,f),
                branch_length_Y(ts,A[0][0],A[1][0],A[1][1]) )
        self.assertEqual(
                branch_stats_node_iter(ts,A,f,method='length'),
                branch_length_Y(ts,A[0][0],A[1][0],A[1][1]) )

    def test_pairwise_diversity(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        self.check_pairwise_diversity(ts)

    def test_Y_stat(self):
        ts = msprime.simulate(10, random_seed=self.random_seed)
        self.check_Y_stat(ts)

    def test_case_1(self):
        ## With mutations:
        #
        # 1.0          6
        # 0.7         / \                                    5
        #            /   X                                  / \
        # 0.5       X     4                4               /   4
        #          /     / \              / \             /   X X
        # 0.4     X     X   \            X   3           X   /   \
        #        /     /     X          /   / X         /   /     \
        # 0.0   0     1       2        1   0   2       0   1       2
        #          (0.0, 0.2),        (0.2, 0.8),       (0.8, 1.0)
        #
        records = [
            msprime.CoalescenceRecord(
                left=0.2, right=0.8, node=3, children=(0,2), 
                time=0.4, population=0),
            msprime.CoalescenceRecord(
                left=0.0, right=0.2, node=4, children=(1,2), 
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.2, right=0.8, node=4, children=(1,3), 
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.8, right=1.0, node=4, children=(1,2), 
                time=0.5, population=0),
            msprime.CoalescenceRecord(
                left=0.8, right=1.0, node=5, children=(0,4), 
                time=0.7, population=0),
            msprime.CoalescenceRecord(
                left=0.0, right=0.2, node=6, children=(0,4),  
                time=1.0, population=0)]

        ll_ts = _msprime.TreeSequence()
        ll_ts.load_records(records)
        ts = msprime.TreeSequence(ll_ts)

        ts.set_mutations([ 
            (0.1, 0),
            (0.15, 0),
            (0.15,1),
            (0.05,4),
            (0.1,2),
            (0.3,1),
            (0.6,2),
            (0.9,0),
            (0.95,1),
            (0.95,2) ])

        self.check_pairwise_diversity(ts)
        self.check_Y_stat(ts)

        # divergence between 0 and 1 
        A = [ [0], [1] ]
        def f(x):
            return (x[0]>0) != (x[1]>0)

        # branch lengths:
        right_ans = 2*( 1 * (0.2-0) + 0.5 * (0.8-0.2) + 0.7 * (1.0-0.8) )
        self.assertEqual(branch_stats(ts,A,f),right_ans)
        self.assertEqual(branch_stats_node_iter(ts,A,f),right_ans)

        # Y-statistic for (0/12)
        A = [ [0], [1,2] ]
        def f(x):
            return ( (x[0]==1) and (x[1]==0) ) or ( (x[0]==0) and (x[1]==2) )
        # branch lengths:
        right_ans = 0.2*(1 + 0.5) + 0.6*(0.4) + 0.2*(0.7+0.2)
        self.assertEqual(branch_length_Y(ts,0,1,2),right_ans)
        self.assertEqual(branch_stats(ts,A,f),right_ans)
        self.assertEqual(branch_stats_node_iter(ts,A,f),right_ans)
