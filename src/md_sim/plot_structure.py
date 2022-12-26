import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import logging
import sys
import pandas as pd

logger                  = logging.getLogger()

def load_hdf5(hdf_file):
    f   = h5.File(hdf_file, "r")
    logger.info(f"File loading complete.")
    equiliration, expansion = pd.DataFrame(), pd.DataFrame()
    logger.info(f"Datasets created.")
    dfs    = [equiliration, expansion]
    dict    = {0:'ThermoData/equilibration', 1:'ThermoData/expansion'}
    steps   = [0,0]   
    n_atoms = [0,0] 
    for purpose in [0,1]:
        keys    = []
        steps[purpose]   =   f[dict[purpose]].attrs['n_steps']
        n_atoms[purpose]   =   f[dict[purpose]].attrs['n_atoms']
        for key in f[dict[purpose]].keys():
            keys.append(key)
            dfs[purpose][key]    = f[dict[purpose]][key][()]
        logger.info(f'Set dataframe columns  of {dict[purpose]} as: {keys}')
    return dfs, steps, n_atoms

def diamond_percentage(hdf_file, n_steps_ex):
    D_perc  =   np.empty(n_steps_ex)
    f   = h5.File(hdf_file, "r")
    n_atoms = f['ThermoData/equilibration'].attrs['n_atoms']
    logger.info(f'Loaded {n_atoms} atoms.')
    maxDiamond   = n_atoms - np.sum(f['StructureAnalysis/StructureType'][120] == 0)
    print(np.max(f['StructureAnalysis/StructureType'][10000]))
    f.close()
    ### structure type 0: unordered, 1: diamond, 2: diamond first neighbour, 3: diamond second neighbour ###
    print(maxDiamond)
    print(n_atoms - maxDiamond)

if __name__=="__main__":    
    logger.setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) == 0:
        logger.error('Please provide a .h5 file!')
    dfs, steps, n_atoms =   load_hdf5(sys.argv[1])
    diamond_percentage(sys.argv[1], steps[1])