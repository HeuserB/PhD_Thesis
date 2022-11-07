#!/bin/sh
squeue -u bh326 > log
MAXJOBS=2
RUNNINGJOBS=$ expr $(sed -n '$=' log) - 1
rm log
LOGDIR=/data/bh326/LAMMPS/free_expansion_dunn
python file_comparison.py $LOGDIR