import numpy as np                                     # Matlab like syntax for linear algebra and functions
import re
import sys
import os
import imageio

sys.path.append("../xrd")
from combine_quads import merge_four_quads
quad_scale = np.loadtxt("../../.data_LW03/instprm/quad_scale_LW03.txt")

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
        data = imageio.imread(file);
        return np.flipud(data);
    except:
        print(file)
        print("No such file!")

def merge_quad23(run, dir_data, q2_ai, q3_ai,q2_mask,q3_mask):
    r_q2_tiff = load_tiff(run, 2, dir_data)
    r_q3_tiff = load_tiff(run, 3, dir_data)
    r_q2_result1d = q2_ai.integrate1d(r_q2_tiff,
                    npt=1000,
                    method='csr',
                    unit='2th_deg',
                    correctSolidAngle=True,
                    polarization_factor=0.99,
                    mask=q2_mask)
    r_q3_result1d = q3_ai.integrate1d(r_q3_tiff,
                    npt=1000,
                    method='csr',
                    unit='2th_deg',
                    correctSolidAngle=True,
                    polarization_factor=-0.99,
                    mask=q3_mask)

    merged_xy = merge_four_quads([np.array(r_q2_result1d),np.array(r_q3_result1d)], scale=[quad_scale[3],quad_scale[2]])
    np.savetxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy",merged_xy.T)

def write_SACLA_lineout(run, shot_id, dir_data, q_ai, mask, ref_id):
    pattern_file    =   "fpd_$run.tif"
    file            = dir_data + "/" + re.sub(r'\$run',str(shot_id),pattern_file)
    ref_file        = dir_data + "/" + re.sub(r'\$run',str(ref_id),pattern_file)
    try:
        data = imageio.imread(file)
        data = np.flipud(data)
        ref = imageio.imread(ref_file)
        ref = np.flipud(ref)
    except:
        print(file,ref_file)
        print("No such files!")

    r_result1d      = q_ai.integrate1d(data,
                        npt=1000,
                        method='csr',
                        unit='2th_deg',
                        correctSolidAngle=True,
                        polarization_factor=0.99,
                        mask=mask)
    ref_result1d    = q_ai.integrate1d(ref,
                        npt=1000,
                        method='csr',
                        unit='2th_deg',
                        correctSolidAngle=True,
                        polarization_factor=0.99,
                        mask=mask)

    script_dir = os.path.dirname(__file__)

    np.savetxt(os.path.join(script_dir, f"../../.data_SACLA/lineouts/r{run}.xy") ,np.array(r_result1d).T)
    np.savetxt(os.path.join(script_dir, f"../../.data_SACLA/lineouts/r{run}_ref.xy"),np.array(ref_result1d).T)

def load_SACLA(run, data_dir='../../.data_SACLA/lineouts/',theta_min=0.,theta_max=75.):
    try:
        run_data        = np.loadtxt(f"{data_dir}r{str(run)}.xy")
    except:
        print(f"Run data for run {run} does not exist in directory {data_dir}!")
    try:
        ref_data        = np.loadtxt(f"{data_dir}r{str(run)}_ref.xy")
    except:
        print(f"Reference data for run {run} does not exist in directory {data_dir}!")
    
    mask       = np.logical_and(run_data[:,0] > theta_min,run_data[:,0] < theta_max)

    return run_data[mask], ref_data[mask]