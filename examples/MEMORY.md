I ran several simulations from the command line for reps `$NSIM = 1 10 100 1000 20000`, using the command:

```sh
/usr/bin/time -v python examples.py -T 21 -N 100 -r 0.01 -l 10 -k 10 -n $NSIM
```

I recorded the maximum resident size in `memory.tsv`. The data:

```sh
reps	max_memory_kbytes
1	54380
10	62216
100	88468
1000	354252
2000	657656
```

Although the script calls a function that both exits and includes `del pop` and
`del rc`, the data (i.e., comparing the memory use from 1000 to 2000 reps)
indicate that total memory use increases linearly with the number of reps.

## More
the provenance output by the script was

```python
Namespace(chrom_length=100, gamma_alpha=0.23, gamma_beta=5.34, generations=20, logfile='-', neut_mut_rate=1e-07, nsamples=10, nselloci=10, nsims=2000, outfile=None, popsize=100, recomb_rate=0.01, sel_mut_rate=1e-07, selloci_file='sel_loci.txt', simplify_interval=500, treefile=None)
```
