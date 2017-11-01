# coding: utf-8
import msprime
ts = msprime.load('neutral_ts')
tsc = msprime.TreeStatCalculator(ts)
print('mean pairwise tmrca matrix among all 6 haploid samples for 10 windows on genome')
for t in tsc.mean_pairwise_tmrca(windows=[float(i) for i in range(10)],
                                 sample_sets=[[0,1,2,3,4,5], [0,1,2,3,4,5]]):
    print(t)
