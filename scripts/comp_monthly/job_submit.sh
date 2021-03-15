#!/bin/bash
#PBS -q regular 
#PBS -A P48500028
#PBS -N monthly 
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -j oe 

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

#mpirun -n 36 ./comp_monthly.py
./comp_monthly_gcm_met.py
