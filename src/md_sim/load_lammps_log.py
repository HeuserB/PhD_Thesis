#from dataclasses import dataclass
#from distutils.command.config import dump_file
#from ase import Atoms
#from matplotlib import gridspec
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

#import sys
import numpy as np
import re
NA = np.newaxis
from scipy.optimize import curve_fit
import pandas as pd
import os.path

import matplotlib.pyplot as plt

start_pattern       = "Step\s+[^\n]*\n"
end_pattern         = "Loop\s+time\s+of"
timestep_pattern    = re.compile('\ntimestep\s+([0-9]+.[0-9]+)\s')
atom_pattern        = re.compile('\nCreated ([0-9]*) atoms')
units_pattern       = re.compile('units (.*)\n')
run_pattern         = re.compile('\nrun\s+([0-9]+)\s')
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
    #[print(f'i is: {i}') for i in timestep_pattern.findall(txt)] 
    dt                  = [np.float(i) for i in timestep_pattern.findall(txt)]
    n_atoms             = np.sum([np.int(i) for i in atom_pattern.findall(txt)])
    #print(f"A total of {n_atoms} were created.")
    unit                = units_pattern.findall(txt)[0]
    if unit == 'box':
        unit = 'metal'


    if len(dt) == 1 :
        dt              = np.repeat(dt,num_of_sets)
        #print(f"All runs have the same timestep timesteps: {dt}")
    elif len(dt) == num_of_sets:
        print(f"Different timesteps are {dt}")
    elif len(dt)==0:
        print(f"Timestep retreived from read dump. Manually setting to 0.0001 ps")
        dt          =    0.0001
    else:
        print(f"Timesteps do not match run numbers!")
        return
    run_steps           = [int(i) for i in run_pattern.findall(txt)]
    dfs                 = []
    
    for i in range(num_of_sets):
        txt_w       = re.split(start_pattern,txt)[i+1]
        keys        = ["Step"]
        [keys.append(j) for j in re.split('\s+',(re.split('\n',re.split("Step\s",txt)[i+1])[0]).strip())]
        txt_w           = re.split(end_pattern, txt_w)[0]
        txt_w           = txt_w.strip()
        data            = np.fromstring(txt_w,sep="\t")
        data            = data.reshape(-1,len(keys))
        df              = pd.DataFrame(data=data, columns=keys)

        dfs.append(df)

    return dfs, dt * time_unit[unit], n_atoms, run_steps

def energy_drift(df : pd.DataFrame, timestep : float, n_atoms : int):
    """
    :param df:            DataFrame of the run
    :param timestep:      timestep for the run in seconds
    :param n_atoms:       number of atoms in the simualtion  
    :return:              maximal energy drift over time unit for the whole run
    """
    energies            = df['TotEng'].to_numpy()
    ### calculate forward difference ###
    deltaE              = (np.abs((energies[1:] - energies[:-1])))/timestep

    time                = df["Step"] * timestep
    E_tot               = df['TotEng']

    try:
        linear_model        = np.polyfit(time[200:],E_tot[200:],1)
    except:
        print("Not 200 thermo entries in dataset. Switching to last 20")
        linear_model        = np.polyfit(time[20:],E_tot[20:],1)
    drift               = linear_model[0]/n_atoms
    #print(f'Drift is {drift* 1e-12} eV/ps/atom')
    drift_kJ            = drift * 1.6021774232052328e-22
    drift_kJ_s_mol      = drift_kJ * 6.02214076e23
    drift_kJ_ns_mol     = drift_kJ_s_mol * 1e-9
    #print(f'Drift is {drift_kJ_ns_mol} kJ/ns/mol')

    return drift

def load_run(file : str):
    grid_pattern                = re.compile('([0-9])_([0-9]+)')
    subrun_pattern              = re.compile('([0-9]+)K_([0-9]*x){2}([0-9]+)[_restart]*.log$')
    if re.match('.*restart.*',file):
        print('File is a restart file!')
        return
    dir                         = os.path.dirname(file)
    print(dir)
    grid_const                  = {}
    thermo_dfs                  = []
    tmp                         = (grid_pattern.findall(dir))[0]
    a                           = float(tmp[0]) + float(tmp[1])*1e-2
    grid_const.update({dir:a})
    data                        = pd.DataFrame(data=None, index=None, columns=['lattice_constant','temperature','cells','timestep_eq',\
                                    'timestep_ex','drift','n_atoms','dE','dumpfile','run_steps_eq', 'run_steps_ex'])
    dump_file                   = re.sub('.log','',file)
    T, cells                    = int(subrun_pattern.findall(file)[0][0]), int(re.split('x',subrun_pattern.findall(file)[0][1])[0])
    dfs, dt, n_atoms, run_steps = load_thermo(file)
    thermo_dfs.append(dfs)
    drift                       = energy_drift(dfs[0], dt[0], n_atoms)
    data.loc[len(data.index)]   = [grid_const[dir],T, cells, dt[0], dt[-1], drift, n_atoms, \
                                    (dfs[0]['TotEng'].values)[200] - (dfs[0]['TotEng'].values)[-1], dump_file,\
                                    run_steps[0], run_steps[-1]]
    return data, thermo_dfs

def load_all_runs(dir : str):
    """
    :param dir:             directory where runs should be searched 
    :return:                pandas dataframe with all run information
    """
    run_pattern             = re.compile('grid_[0-9]_[0-9]+')
    subrun_pattern          = re.compile('([0-9]+)K_([0-9]*x){2}([0-9]+)(_[0-9]+ns)*.log')
    grid_pattern            = re.compile('([0-9])_([0-9]+)')
    dirs                    = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    run_dirs, grid_const    = [], {}
    thermo_dfs              = []
    data                    = pd.DataFrame(data=None, index=None, columns=['lattice_constant','temperature','cells','timestep_eq',\
                                'timestep_ex','drift','n_atoms','dE','dumpfile','run_steps_eq', 'run_steps_ex'])

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
            strain_time_pattern         =   re.compile('.*100ns.*')
            if strain_time_pattern.match(file):
                pass
            else:
                continue
            if subrun_pattern.match(file):
                dump_file                   = re.sub('.log','.dump',file)
                print(dump_file)
                T, cells                    = int(subrun_pattern.findall(file)[0][0]), int(re.split('x',subrun_pattern.findall(file)[0][1])[0])
                dfs, dt, n_atoms, run_steps = load_thermo(os.path.join(subdir,file))
                thermo_dfs.append(dfs)
                #print(dfs[0]["Step"] * dt[0])
                drift                       = energy_drift(dfs[0], dt[0], n_atoms)
                try:
                    data.loc[len(data.index)]   = [grid_const[d],T, cells, dt[0], dt[1], drift, n_atoms, \
                                            (dfs[0]['TotEng'].values)[200] - (dfs[0]['TotEng'].values)[-1], os.path.join(subdir,dump_file),\
                                            run_steps[0], run_steps[1]]
                except:
                    data.loc[len(data.index)]   = [grid_const[d],T, cells, dt[0], dt[1], drift, n_atoms, \
                                            (dfs[0]['TotEng'].values)[20] - (dfs[0]['TotEng'].values)[-1], os.path.join(subdir,dump_file),\
                                            run_steps[0], run_steps[1]]
                print(f"Directory: {d} file : {file}, grid constant: {grid_const[d]} T: {T}, cells: {cells}")
    return data, thermo_dfs