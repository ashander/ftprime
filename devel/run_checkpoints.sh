#!/bin/bash


SEED=$RANDOM

for LEVEL in 1 2 3
do
    OUTDIR=sim_${SEED}_level_${LEVEL} 
    mkdir -p $OUTDIR
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
        ./checkpoint_simupop.py -l ${LEVEL} -T 1000 -N 400 -L 1e6 -r 1e-6 -k 20 -A 10 -U 1e-6 -d $SEED -o $OUTDIR \
            &> $OUTDIR/time.log
done

tail -n 1 sim_${SEED}_level_*/time.log
