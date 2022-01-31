import numpy as np                                     # Matlab like syntax for linear algebra and functions
import re
import sys
import os
import imageio

NA = np.newaxis

''' 
Define the data directories and the RegEx for how the runs look like
'''

# define lineout dir
dir_lineout = '../lineouts/Single_shock_PET100/'

# define the naming pattern for a lineout file
pattern_run = re.compile(".*r[0-9]+.*.xy")

# define the pattern for a background file
pattern_bkg = re.compile(".*r[0-9]*.*ambient.*.xy")

# define a pattern that should be ignored
pattern_ign = re.compile('^(?!Q[0-9]).')

def list_runs(run_dir, pattern_run, pattern_bkg, pattern_ign):
    '''
    List all runs that are found a directory

    INPUT:
    run_dir: str ; directory where to look for runs
    run_patter: re.pattern ; pattern what a run file (eg. '*.tiff', '*.xy') looks like 
    pattern_bkg: re.pattern ; pattern what a cold run file (eg. '*.tiff', '*.xy') looks like (without laser shock)
    ign_patter: re.pattern ; pattern that is to be ignored during the scan

    OUTPUT:
    files_drive: list[str] ; filenames of the drive runs
    files_ambient: list[str]  ; filenames of the backgrounbd (cold) runs
    np.array(runs): list[int] ; integer numbers of the runs that where found

    '''
    files_drive = []
    files_ambient = []
    runs = []

    # iterate through all files in the directory
    for file in os.listdir(run_dir):
        # if a 'run' is found
        if re.match(pattern_run, file):
            if re.match(pattern_ign, file):
                run = re.split("_",file)[0]
                run = int(re.findall(r'\d+', run)[0])
                for file_amb in os.listdir(run_dir):
                    if re.match(pattern_ign, file_amb):
                        if re.match(pattern_bkg, file_amb):
                            # the background files are taken before the run so the number is 'run -1'
                            if re.search(str(int(run) - 1), file_amb):
                                files_drive.append(file)
                                files_ambient.append(file_amb)
                                runs.append(run)
    return files_drive, files_ambient, np.array(runs)

def load_run(run, files_drive, files_ambient, runs, run_dir):
    '''
    load a run file and the corresponding background as an numpy array

    INPUT:
    run: int ; the number of the run to load
    run_dir: str ; directory where to look for runs
    files_drive: list[str] ; filenames of the drive runs
    files_ambient: list[str]  ; filenames of the backgrounbd (cold) runs
    np.array(runs): list[int] ; integer numbers of the runs that where found

    OUTPUT:
    drive_data: np.array(N,2) ; lineout data of the run with [:,0]=two theta and [:,1]=intensity in arb. units
    ambient_data: np.array(N,2)  ; lineout data of the background with [:,0]=two theta and [:,1]=intensity in arb. units
    '''

    run_id = np.where(runs == run)[0][0]
    print(run_id)
    drive_data = np.loadtxt(os.path.join(run_dir,files_drive[run_id]))
    ambient_data = np.loadtxt(os.path.join(run_dir,files_ambient[run_id]))
    return ambient_data, drive_data

def list_runs(dir_data, pattern_runs):
    '''
    load a run file and the corresponding background as an numpy array

    INPUT:
    run_dir: str ; directory where to look for runs
    pattern_runs: pattern ; a pattern what a run folder with the *.tiff files looks like

    OUTPUT:
    run_data: np.array(N) ; integers of all available runs in the tiff
    '''
    runs = []
    # iterate through all files in the directory
    for file in os.listdir(dir_data):
        if re.match(pattern_runs, file):
            run_int = ''.join(re.findall(r'\d+', file))
            runs.append(int(run_int))
    return np.array(runs)

def load_tiff(run, quad, dir_data,pattern_folder="r$run_bkgCorrected",pattern_file="ePix10k_Quad$quad_r$run_0.tiff"):
    file = dir_data + "/" + re.sub(r'\$run',str(run),pattern_folder) + "/" + re.sub('\$quad',str(quad),re.sub(r'\$run',str(run),pattern_file))
    try:
        data = imageio.imread(file)
        return np.flipud(data)
    except:
        print(file)
        print("No such file!")
    