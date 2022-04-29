
from dataclasses import dataclass
from distutils.command.config import dump_file
from ase import Atoms
from matplotlib import gridspec
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

import sys
import numpy as np
import re
NA = np.newaxis
from scipy.optimize import curve_fit
import pandas as pd
import os.path

import matplotlib.pyplot as plt

start_pattern       = "Step\s+[^\n]*\n"
end_pattern         = "Loop\s+time\s+of"
timestep_pattern    = re.compile('\ntimestep\s+([0-9]*.[0-9]*)\s')
atom_pattern        = re.compile('\nCreated ([0-9]*) atoms')
units_pattern       = re.compile('units (.*)\n')
time_unit           = {'metal' : 1.0e-12, 'real' : 1.0e-15,'si' : 1., 'cgs' :   1., 'electron' : 1e-15 , 'micro': 1.0e-6}


def load_thermo(file : str):
    """
    :param file:            name of the logfile
    :return:                data, timesteps, number of atoms
    """

    f                   = open(file)
    txt                 = f.read()
    f.close()

    num_of_sets         = len(re.findall(start_pattern, txt))
    #print(f"Found {num_of_sets} different run commands")
    dt                  = [np.float(i) for i in timestep_pattern.findall(txt)]
    n_atoms             = np.sum([np.int(i) for i in atom_pattern.findall(txt)])
    #print(f"A total of {n_atoms} were created.")
    unit                = units_pattern.findall(txt)[0]
    #print(f'Unit system is {unit}')

    if len(dt) == 1 :
        dt              = np.repeat(dt,num_of_sets)
        #print(f"All runs have the same timestep timesteps: {dt}")
    elif len(dt) == num_of_sets:
        print(f"Different timesteps are {dt}")
    else:
        print(f"Timesteps do not match run numbers!")
        return
   
    dfs                 = []
    
    for i in range(num_of_sets):
        txt_w       = re.split(start_pattern,txt)[i+1]
        keys        = ["Step"]
        [keys.append(i) for i in re.split('\s+',(re.split('\n',re.split("Step\s",txt)[i+1])[0]).strip())]

        txt_w           = re.split(end_pattern, txt_w)[0]
        txt_w           = txt_w.strip()
        data            = np.fromstring(txt_w,sep="\t")
        data            = data.reshape(-1,len(keys))
        df              = pd.DataFrame(data=data, columns=keys)

        dfs.append(df)

    return dfs, dt * time_unit[unit], n_atoms

def energy_drift(df : pd.DataFrame, timestep : float, n_atoms : int):
    """
    :param df:            DataFrame of the run
    :param timestep:      timestep for the run in seconds
    :param n_atoms:       number of atoms in the simualtion  
    :return:              maximal energy drift over time unit for the whole run
    """
    energies        = df['TotEng'].to_numpy()
    ### calculate forward difference ###
    deltaE          = (np.abs((energies[1:] - energies[:-1])))/timestep

    time                = df['Step'] * timestep
    E_tot               = df['TotEng']

    linear_model        = np.polyfit(time[200:],E_tot[200:],1)
    drift               = linear_model[0]/n_atoms
    #print(f'Drift is {drift* 1e-12} eV/ps/atom')
    drift_kJ            = drift * 1.6021774232052328e-22
    drift_kJ_s_mol      = drift_kJ * 6.02214076e23
    drift_kJ_ns_mol     = drift_kJ_s_mol * 1e-9
    #print(f'Drift is {drift_kJ_ns_mol} kJ/ns/mol')

    return drift

def load_all_runs(dir : str):
    """
    :param dir:             directory where runs should be searched 
    :return:                maximal energy drift over time unit for the whole run
    """
    run_pattern             = re.compile('grid_[0-9]_[0-9]+')
    subrun_pattern          = re.compile('([0-9]+)K_([0-9]*x){2}([0-9]+).log$')
    grid_pattern            = re.compile('([0-9])_([0-9]+)')
    dirs                    = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    run_dirs, grid_const    = [], {}
    data                    = pd.DataFrame(data=None, index=None, columns=['lattice_constant','temperature','cells','timestep_eq',\
                                'timestep_ex','drift','n_atoms','dE','dumpfile'])

    for d in dirs:
        if run_pattern.match(d):
            run_dirs.append(d)
            tmp             = (grid_pattern.findall(d))[0]
            a               = float(tmp[0]) + float(tmp[1])*1e-2
            grid_const.update({d:a})
    
    for d in run_dirs:
        subdir              = os.path.join(dir, d)
        files               = os.listdir(subdir)
        for file in files:
            if subrun_pattern.match(file):
                dump_file                   = re.sub('.log','',file)
                T, cells                    = int(subrun_pattern.findall(file)[0][0]), int(re.split('x',subrun_pattern.findall(file)[0][1])[0])
                dfs, dt, n_atoms            = load_thermo(os.path.join(subdir,file))
                drift                       = energy_drift(dfs[0], dt[0], n_atoms)
                data.loc[len(data.index)]   = [grid_const[d],T, cells, dt[0], dt[1], drift, n_atoms, \
                                            (dfs[0]['TotEng'].values)[200] - (dfs[0]['TotEng'].values)[-1], os.path.join(subdir,dump_file)]
                #print(f"Directory: {d} file : {file}, grid constant: {grid_const[d]} T: {T}, cells: {cells}")
    return data