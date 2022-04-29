import numpy as np
import re

f = open("in.free_expansion_airebo_m_template",'r')
in_txt = f.read()
f.close()

f = open("lammps.SLURM_free_expansion_airebo_m",'r')
slurm_txt = f.read()
f.close()

#liquidline = np.loadtxt('../../../../W_PhD_Data/Phase_data/Liquid_line.csv', delimiter=',')

def input_temp(in_txt, slurm_txt, temp, grid_constant, n_cells, t_init=17000.0,seed=123456789):
    if (n_cells%2 != 0):
        print('Use an even number of grid cells!')
        return
    x_min, x_max, y_min, y_max, z_min, z_max = np.repeat(n_cells * grid_constant /2.,6)*np.array([-1.,1.,-1.,1.,-1.,1.])
    grid = re.sub('\.','_',str('%.2f' %grid_constant))
    temp = round(temp,1)
    in_txt = re.sub(r'\$t_init',str(t_init),in_txt)
    in_txt = re.sub(r'\$temperature',str(temp),in_txt)
    in_txt = re.sub(r'\$x_min',str(x_min),in_txt)
    in_txt = re.sub(r'\$x_max',str(x_max),in_txt)
    in_txt = re.sub(r'\$y_min',str(y_min),in_txt)
    in_txt = re.sub(r'\$y_max',str(y_max),in_txt)
    in_txt = re.sub(r'\$z_min',str(z_min),in_txt)
    in_txt = re.sub(r'\$z_max',str(z_max),in_txt)
    in_txt = re.sub(r'\$grid_constant',str(grid_constant),in_txt)
    in_txt = re.sub(r'\$seed', str(seed), in_txt)

    #velfile = 'dump.melt_velocity_' + str(int(temp)) + "K_" + str(int(press/1e5)) + "GPa_PP"
    dumpfile = '/data/bh326/LAMMPS/free_expansion_airebo_m/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells)
    logfile = '/data/bh326/LAMMPS/free_expansion_airebo_m/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) +'.log'

    in_txt = re.sub(r'\$dumpfile',dumpfile,in_txt)
    in_txt = re.sub(r'\$logfile',logfile,in_txt)

    outfile = 'inputs_free_expansion/in.free_expansion_airebo_m_' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid)
    
    f = open(outfile, 'w')
    f.write(in_txt)
    f.close()

    slurm_txt = re.sub(r'\$temp',str(int(temp)),slurm_txt)
    slurm_txt = re.sub(r'\$n_cell',str(n_cells),slurm_txt)
    infile = re.sub('inputs_free_expansion/', '', outfile)
    slurm_txt = re.sub(r'\$infile',infile,slurm_txt)
    
    outfile_slurm  = 'inputs_free_expansion/lammps.free_expansion_airebo_m_' + str(temp) +'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid)

    f = open(outfile_slurm, 'w')
    f.write(slurm_txt)
    f.close()

grids = [3.35]
temperatures = [1000,2000,3000,4000]
for grid in grids:
    for temperature in temperatures:
        input_temp(in_txt, slurm_txt, temperature, grid, 8, t_init=temperature,seed=123456789)