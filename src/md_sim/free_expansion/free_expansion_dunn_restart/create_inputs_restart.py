
import numpy as np
import re

equlibration_txt    = "# Equilibrate the diamond at a certain temperature\nfix 1 all nvt temp $temperature $temperature $(dt*100.0)\n\ntimestep	${timestep} # 0.0001 #when using metal units [ps]\nrun		$remaining_steps # equilibration timestes\nvariable E_equil equal etotal # store the Energy of the equilibrated system as a variable\n\nunfix 1\n"
expansion_txt       = "# Create a large box with fixed bounadary conditions\nchange_box all boundary f f f x final -200 200 y final -200 200 z final -200 200\n# ignore lost atoms along the way\nthermo_modify lost ignore flush yes\n\nfix 2 all nve\n\nrun $remaining_steps # expand\n"

def input_temp(in_txt, slurm_txt, temp, grid_constant, n_cells, last_step, t_init=17000.0,seed=123456789):
    print('Writing files now with n_cells:', n_cells, 'temp:', temp, 'grid_constant:', grid_constant)
    if (n_cells%2 != 0):
        print('Use an even number of grid cells!')
        return

    print('Last step was: ', last_step)
    if last_step < 20000:
        phase   =   "equilibration"
    else:
        phase   =   "expansion"
    print("Continuing with: ", phase, " phase.")

    if phase == "equilibration":
        remaining_steps     = 20000 - last_step
        in_txt              += re.sub(r'\$remaining_steps', str(int(remaining_steps)), equlibration_txt)
        in_txt              += re.sub(r'\$remaining_steps', '200000', expansion_txt)
    if phase == "expansion":
        remaining_steps     = 20000 + 200000 - last_step
        in_txt              += re.sub(r'\$remaining_steps', str(int(remaining_steps)), expansion_txt)

    grid = re.sub('\.','_',str('%.2f' %grid_constant))
    temp = round(temp,1)
    x_min, x_max, y_min, y_max, z_min, z_max = np.repeat(n_cells * grid_constant /2.,6)*np.array([-1.,1.,-1.,1.,-1.,1.])

    old_dump = '/data/bh326/LAMMPS/free_expansion_dunn/grid_' + str(grid) +'/' + str(int(temp)) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells)
    logfile = '/data/bh326/LAMMPS/free_expansion_dunn/grid_' + str(grid) +'/' + str(int(temp)) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) +'.log'

    in_txt      = re.sub(r'\$old_dump',old_dump ,in_txt)
    in_txt      = re.sub(r'\$logfile',logfile ,in_txt)
    new_dump    = old_dump + '_restart_step' + str(int(last_step))
    in_txt      = re.sub(r'\$new_dump',new_dump,in_txt)
    in_txt      = re.sub(r'\$last_step',str(last_step),in_txt)
    in_txt      = re.sub(r'\$t_init',str(t_init),in_txt)
    in_txt      = re.sub(r'\$temperature',str(temp),in_txt)
    in_txt      = re.sub(r'\$x_min',str(x_min),in_txt)
    in_txt      = re.sub(r'\$x_max',str(x_max),in_txt)
    in_txt      = re.sub(r'\$y_min',str(y_min),in_txt)
    in_txt      = re.sub(r'\$y_max',str(y_max),in_txt)
    in_txt      = re.sub(r'\$z_min',str(z_min),in_txt)
    in_txt      = re.sub(r'\$z_max',str(z_max),in_txt)
    in_txt      = re.sub(r'\$grid_constant',str(grid_constant),in_txt)
    in_txt      = re.sub(r'\$seed', str(seed), in_txt)

    outfile = 'restart_inputs/in.free_expansion_dunn_' + str(int(temp)) + 'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid) + '_restart'
    
    f = open(outfile, 'w')
    f.write(in_txt)
    f.close()

    slurm_txt = re.sub(r'\$temp',str(int(temp)),slurm_txt)
    slurm_txt = re.sub(r'\$n_cell',str(n_cells),slurm_txt)
    infile = re.sub('restart_inputs/', '', outfile)
    slurm_txt = re.sub(r'\$infile', infile, slurm_txt)
    
    outfile_slurm  = 'restart_inputs/lammps.free_expansion_dunn_' + str(int(temp)) +'K_' +str(n_cells) +'x' + str(n_cells) +'x' + str(n_cells) + 'grid_' + str(grid)

    f = open(outfile_slurm, 'w')
    f.write(slurm_txt)
    f.close()

#grids = [3.35,3.4,3.45,3.5]
#temperatures = [1000,2000,3000,4000]
#for grid in grids:
#    for temperature in temperatures:
#        input_temp(in_txt, slurm_txt, temperature, grid, 8, t_init=temperature,seed=123456789)