import msprime
import ftprime

nloci = int(1e4 / (4 * 100 * 1e-7))

locus_position = list(range(0, nloci))
init_ts = msprime.simulate(2*100,
                           Ne=100,
                           recombination_rate=1e-7/ 2.0,
                           length=max(locus_position))
haploid_labels = [(k,p) for k in range(100)  
                        for p in (0,1)]
node_ids = {x:j for x, j in zip(haploid_labels, init_ts.samples())}
rc = ftprime.RecombCollector(ts=init_ts, node_ids=node_ids,
                     locus_position=locus_position,
                     benchmark=True, mode='text')
lines = '101 31 1 79980 1220688 14497527 18316131\n101 10 1 3639834 18480896\n'
