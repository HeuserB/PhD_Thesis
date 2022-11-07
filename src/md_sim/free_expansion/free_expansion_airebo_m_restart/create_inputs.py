import numpy as np
import re

f = open("in.restart_free_ex",'r')
in_txt = f.read()
f.close()

f = open("lammps.SLURM_free_expansion_airebo_m_restart",'r')
slurm_txt = f.read()
f.close()

#liquidline = np.loadtxt('../../../../W_PhD_Data/Phase_data/Liquid_line.csv', delimiter=',')

def input_temp(in_txt, slurm_txt, temp, grid_constant, n_cells=8, seed=123456789):
    if (n_cells%2 != 0):
        print('Use an even number of grid cells!')
        return
    grid = re.sub('\.','_',str('%.2f' %grid_constant))
    temp = round(temp,1)

    old_dump = '/data/bh326/LAMMPS/free_expansion_airebo_m/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells)
    logfile = '/data/bh326/LAMMPS/free_expansion_airebo_m/grid_' +str(grid) +'/' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) +'_restart.log'

    in_txt      = re.sub(r'\$old_dump',old_dump ,in_txt)
    in_txt      = re.sub(r'\$logfile',logfile ,in_txt)
    new_dump    = old_dump + '_restart'
    in_txt      = re.sub(r'\$new_dump',new_dump,in_txt)

    outfile = 'inputs_free_expansion_restart/in.free_expansion_airebo_m_' + str(temp) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid) + '_restart'
    
    f = open(outfile, 'w')
    f.write(in_txt)
    f.close()

    slurm_txt = re.sub(r'\$temp',str(int(temp)),slurm_txt)
    slurm_txt = re.sub(r'\$n_cell',str(n_cells),slurm_txt)
    infile = re.sub('inputs_free_expansion_restart/', '', outfile)
    slurm_txt = re.sub(r'\$infile',infile,slurm_txt)
    
    outfile_slurm  = 'inputs_free_expansion_restart/lammps.free_expansion_airebo_m_' + str(temp) +'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid)

    f = open(outfile_slurm, 'w')
    f.write(slurm_txt)
    f.close()

grids = [3.50,3.35,3.4]
temperatures = [1000,1000,2000]
for grid, temperature in zip(grids, temperatures):
    input_temp(in_txt, slurm_txt,temp=temperature,grid_constant= grid,n_cells=8)