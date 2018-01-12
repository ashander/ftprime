# Minimal `simuPOP` example

-  [simupop_example.py](simupop_example.py): A simple walk-through of using with simuPOP.

# Single locus example

After this, `tree.hdf5` will contain a tree-sequence with a single selected
loci in the middle of block of length `L`.

```sh
python3 examples.py -N 5000 -L 10000 -g log.txt -T 50 -t tree.hdf5
```

The selection coefficient is drawn from a gamma distribution (parameters can be passed to the script) but one could also edit the script, replacing the callable `GammaDistributedFitness` with a simple function using logic like [this section](https://github.com/petrelharp/ftprime_ms/blob/ed3423c298dc245ca5b94e1f653292ea719cb7ae/sims/run-singlelocus.py#L121-L124) and hard-coding `s`.


# Setup

Follow the steps in [the main README](../README.md#development) to set up a conda environment for `ftprime`.
