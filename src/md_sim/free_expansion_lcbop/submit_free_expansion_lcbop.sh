#!/bin/sh
mkdir inputs_free_expansion
cp run_jobs.sh inputs_free_expansion/
cp run_more.sh inputs_free_expansion/
mkdir inputs_free_expansion/done/
python3 create_inputs.py
tar -czvf free_expansion.tar.gz inputs_free_expansion
rm -r inputs_free_expansion
scp free_expansion.tar.gz bh326@titan1.hpc.uni-rostock.de:/home/bh326/free_expansion_lcbop/
rm free_expansion.tar.gz

ssh bh326@titan1.hpc.uni-rostock.de 'cd ~/free_expansion_lcbop/  && tar -xf free_expansion.tar.gz'