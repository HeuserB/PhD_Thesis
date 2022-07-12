#!/bin/sh
mkdir inputs_free_expansion_restart
cp ../free_expansion_airebo_m/run_jobs.sh inputs_free_expansion_restart/
cp ../free_expansion_airebo_m/run_more.sh inputs_free_expansion_restart/
mkdir inputs_free_expansion_restart/done/
python3 create_inputs.py
tar -czvf free_expansion_restart.tar.gz inputs_free_expansion_restart
rm -r inputs_free_expansion_restart
scp free_expansion_restart.tar.gz bh326@titan1.hpc.uni-rostock.de:/home/bh326/free_expansion_airebo_m/
rm free_expansion_restart.tar.gz

ssh bh326@titan1.hpc.uni-rostock.de 'cd ~/free_expansion_airebo_m/  && tar -xf free_expansion_restart.tar.gz'