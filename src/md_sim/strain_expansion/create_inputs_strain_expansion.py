import numpy as np
import re

f = open("in.strain_expansion_template",'r')
in_txt = f.read()
f.close()

f = open("lammps.SLURM_strain_expansion",'r')
slurm_txt = f.read()
f.close()

#liquidline = np.loadtxt('../../../../W_PhD_Data/Phase_data/Liquid_line.csv', delimiter=',')

def input_temp(in_txt, slurm_txt, temp, grid_constant, potential, expansion_time, n_cells=8, seed=123456789,local=False):
    if (n_cells%2 != 0):
        print('Use an even number of grid cells!')
        return
    pair_style          = {'edip': 'edip', 'airebo':'airebo/morse 3.0 1 1','lcbop':'lcbop'}
    parameter_files     = {'edip': '/home/bh326/POTENTIALS/C.edip', 'airebo_m':'/home/bh326/POTENTIALS/CH.airebo-m','lcbop':'/home/bh326/POTENTIALS/C.lcbop'}
    if local == True:
        parameter_files  = {'edip': '/home/benjamin/Documents/LAMMPS/C.edip'}
    parameter_file      = parameter_files[potential]
    grid = re.sub('\.','_',str('%.2f' %grid_constant))
    temp = round(temp,1)
    in_txt = re.sub(r'\$temperature',str(temp),in_txt)
    in_txt = re.sub(r'\$grid_constant',str(grid_constant),in_txt)
    in_txt = re.sub(r'\$n_cells', str(n_cells), in_txt)
    in_txt = re.sub(r'\$seed', str(seed), in_txt)
    in_txt = re.sub(r'\$potential', potential, in_txt)
    in_txt = re.sub(r'\$parameter_file', parameter_file, in_txt)
    in_txt = re.sub(r'\$expansion_time', str(expansion_time), in_txt)

    #velfile = 'dump.melt_velocity_' + str(int(temp)) + "K_" + str(int(press/1e5)) + "GPa_PP"
    dumpfile = '/data/bh326/LAMMPS/strain_expansion/' + potential +  '/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + '_' + str(expansion_time) + 'ns.dump'
    logfile = '/data/bh326/LAMMPS/strain_expansion/' + potential +  '/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + '_' + str(expansion_time) + 'ns.log'
    if local == True:
        dumpfile = '/home/benjamin/Documents/LAMMPS/strain_expansion/' + potential +  '/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + '_' + str(expansion_time) + 'ns.dump'
        logfile = '/home/benjamin/Documents/LAMMPS/strain_expansion/' + potential +  '/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + '_' + str(expansion_time) + 'ns.log'


    in_txt = re.sub(r'\$dumpfile',dumpfile,in_txt)
    in_txt = re.sub(r'\$logfile',logfile,in_txt)

    outfile = 'inputs_strain_expansion/in.strain_expansion_' + potential + "_" + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid)+ '_' + str(expansion_time) + 'ns'
    
    f = open(outfile, 'w')
    f.write(in_txt)
    f.close()

    slurm_txt = re.sub(r'\$temp',str(int(temp)),slurm_txt)
    slurm_txt = re.sub(r'\$n_cell',str(n_cells),slurm_txt)
    slurm_txt = re.sub(r'\$potential',potential,slurm_txt)
    slurm_txt = re.sub(r'\$grid_constant',re.sub('\.','_',str(grid)),slurm_txt)
    slurm_txt = re.sub(r'\$expansion_time', str(expansion_time), slurm_txt)
    infile = re.sub('inputs_strain_expansion/', '', outfile)
    slurm_txt = re.sub(r'\$infile',infile,slurm_txt)
    
    outfile_slurm  = 'inputs_strain_expansion/lammps.strain_expansion_' + potential + "_" + str(temp) +'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid) + '_' + str(expansion_time) + 'ns'

    f = open(outfile_slurm, 'w')
    f.write(slurm_txt)
    f.close()

grids = [3.2,3.25,3.3,3.35,3.40,3.45]
temperatures = [2000,3000,4000,5000]
potentials = ['edip']
expansion_times = [100]
for potential in potentials:
    for grid in grids:
        for temperature in temperatures:
            for expansion_time in expansion_times:
                input_temp(in_txt, slurm_txt, temperature, grid, potential, expansion_time,local=False)