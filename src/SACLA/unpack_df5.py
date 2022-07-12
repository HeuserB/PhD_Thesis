from cProfile import run
from fileinput import filename
import numpy as np
import h5py
from PIL import Image
import os

from multiprocessing import Pool

srcdir          = '/mnt/TOSHIBA_EXT/HZDR_data/SACLA_hdf/hdf'
hfd5_files      = os.listdir(srcdir)
filenames       = [srcdir + '/' + name for name in hfd5_files]  


def unpack_SAXS(file, SAXS_directory  = "/mnt/TOSHIBA_EXT/HZDR_data/SACLA_2022/stevenson_2022a/SAXS"):
    f               = h5py.File(file , 'r')
    run_number      =   int(f['file_info']['run_number_list'][0])

    run             = f[list(f.keys())[1]]
    
    if 'detector_2d_assembled_1' not in list(run.keys()):
        print(f'No SAXS data for run {run_number}')
        return

    SAXS_detector   = run['detector_2d_assembled_1']
    SAXS_data       = np.array(SAXS_detector[list(SAXS_detector.keys())[1]].get('detector_data'))

    im              = Image.fromarray(SAXS_data)
    filename        = SAXS_directory + "/SAXS_" + str(run_number) + ".tif"

    im.save(filename)

if __name__ == "__main__":  # confirms that the code is under main function
    pool = Pool(os.cpu_count())
    pool.map(unpack_SAXS, filenames)