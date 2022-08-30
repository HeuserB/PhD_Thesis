import sys
import matplotlib.pyplot as plt
#import numpy as np
import os

sys.path.append("../../src/md_sim")
from structure_analysis import combine_hdf5

#store_as_hdf("/mnt/TOSHIBA_EXT/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/free_expansion_airebo_m/", "airebo_m", 100)


written_files   = []
#dir             = '/mnt/TOSHIBA_EXT/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/free_expansion_airebo_m/'
dir             = '/mnt/TOSHIBA_EXT/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/free_expansion_lcbop/'
for file in os.listdir(dir):
    # check only text files
    if file.endswith('.h5'):
        written_files.append(dir + file)

combine_hdf5('/mnt/TOSHIBA_EXT/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/free_expansion_lcbop/lcbop.h5',h5files= written_files, run_group='lcbop')