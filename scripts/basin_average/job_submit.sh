#!/bin/bash
#PBS -q regular 
#PBS -A P48500028
#PBS -N remap
#PBS -l walltime=2:00:00
#PBS -j oe 
#PBS -l select=1:ncpus=6:mpiprocs=6

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

./run_remap_gcm.py
