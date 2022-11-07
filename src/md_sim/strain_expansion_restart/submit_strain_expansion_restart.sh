#!/bin/sh
mkdir inputs_strain_expansion
mkdir inputs_strain_expansion/done/
python3 create_inputs_strain_expansion_restart.py
tar -czvf strain_expansion.tar.gz inputs_strain_expansion
rm -r inputs_strain_expansion
scp strain_expansion.tar.gz bh326@titan1.hpc.uni-rostock.de:/home/bh326/strain_expansion/
rm strain_expansion.tar.gz

ssh bh326@titan1.hpc.uni-rostock.de 'cd ~/strain_expansion/  && tar -xf strain_expansion.tar.gz'