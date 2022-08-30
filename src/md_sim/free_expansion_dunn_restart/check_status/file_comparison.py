#!/usr/bin/python

import sys
import numpy as np
import re
import os

def list_infiles(log_dir):
    logfiles        = []
    for root, dir_names, file_names in os.walk(log_dir):
        for id, filename in enumerate(file_names):
            if re.search('.log',filename):
                logfiles.append(os.path.join(root,filename))
    return logfiles

def check_step(logfile):
    f           = open(logfile)
    text        = f.read()
    f.close()
    text        = text.strip()
    last_line   = re.split('\n',text)[-1]
    try:
        step        = int(re.split('\s+', last_line.strip())[1])
    except:
        print("Initial simulation failed!")
        step        = 0
    return int(step)

def list_steps(log_dir):
    logfiles        = list_infiles(log_dir)
    steps           = [None] * len(logfiles) 
    for id, logfile in enumerate(logfiles):
        print(logfile)
        step        = check_step(logfile)
        steps[id]   = step

    return logfiles, steps

def write_input_restart(failed_logfile, step):
    from create_inputs_restart import input_temp

    f = open("in.restart_free_expansion_dunn_template",'r')
    in_txt = f.read()
    f.close()

    f = open("lammps.SLURM_free_expansion_dunn_restart",'r')
    slurm_txt = f.read()
    f.close()

    print(failed_logfile)
    T               = re.findall('[0-9]*K',failed_logfile)[0]
    T               = int(''.join([a for a in T if a.isdigit()]))

    grid            = re.findall('grid_[0-9]*_[0-9]*',failed_logfile)[0]
    grid            = re.split('_',grid)[-2:]
    grid_constant   = float(grid[0] + '.' + grid[1]) 
    
    n_cells         = re.findall('[0-9]x',failed_logfile)[0]
    n_cells         = int(''.join([a for a in n_cells if a.isdigit()]))

    print(T, grid_constant, n_cells)
    input_temp(in_txt, slurm_txt, T, grid_constant, n_cells, step,t_init=17000.0,seed=123456789)

if __name__ == "__main__":
    failed_logfiles, steps = list_steps(sys.argv[1])
    write_input_restart(failed_logfiles[0],steps[0])