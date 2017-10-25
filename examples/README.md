# Single locus example

After this, `tree.hdf5` will contain a tree-sequence with a single selected
loci in the middle of block of length `L`.

```sh
python3 examples.py -N 5000 -L 10000 -g log.txt -T 50 -t tree.hdf5
```

The selection coefficient is drawn from a gamma distribution (parameters can be passed to the script) but one could also edit the script, replacing the callable `GammaDistributedFitness` with a simple function using logic like [this section](https://github.com/petrelharp/ftprime_ms/blob/ed3423c298dc245ca5b94e1f653292ea719cb7ae/sims/run-singlelocus.py#L121-L124) and hard-coding `s`.

# Setup

```sh
conda config --add channels conda-forge
conda update -q conda
conda info -a
if [[ ! -d "/opt/conda/envs/ftprime" ]] ; then
	conda create -q -n ftprime python=3.5
fi
conda env update -q -n ftprime -f ../environment.yml
source activate ftprime
pip install --editable=git+https://github.com/ashander/ftprime.git@master#egg=ftprime
```
