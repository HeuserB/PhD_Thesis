#!/bin/sh
mkdir inputs_free_expansion_restart

cp create_inputs_restart.py in.restart_free_expansion_dunn_template lammps.SLURM_free_expansion_dunn_restart check_status/check_failed.sh check_status/file_comparison.py inputs_free_expansion_restart

mkdir inputs_free_expansion_restart/restart_inputs
cp run_jobs.sh run_more.sh inputs_free_expansion_restart/restart_inputs
mkdir inputs_free_expansion_restart/restart_inputs/done/

tar -czvf free_expansion_restart.tar.gz inputs_free_expansion_restart
rm -r inputs_free_expansion_restart
scp free_expansion_restart.tar.gz bh326@titan1.hpc.uni-rostock.de:/home/bh326/free_expansion_dunn/
rm free_expansion_restart.tar.gz

ssh bh326@titan1.hpc.uni-rostock.de 'cd ~/free_expansion_dunn/  && tar -xf free_expansion_restart.tar.gz'