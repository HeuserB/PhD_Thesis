#!/bin/sh
mkdir inputs_strain_expansion
cp run_jobs.sh inputs_strain_expansion/
cp run_more.sh inputs_strain_expansion/
mkdir inputs_strain_expansion/done/
python3 create_inputs_strain_expansion.py
mv inputs_strain_expansion /home/benjamin/Documents/LAMMPS/strain_expansion/