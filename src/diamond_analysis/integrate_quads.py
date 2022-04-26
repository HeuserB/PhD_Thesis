from xrdfit.spectrum_fitting import PeakParams, FitSpectrum
import re
import sys

import pyFAI
import numpy as np
import imageio
import matplotlib.pyplot as plt

from load_data import list_runs, load_tiff

def integrate1d(run, quad, dir, poni, mask):
    '''
    run : int ; number of the run to load
    quad : [0,1,2,3] ; the quad to be laoded
    dir : str ; datra directory
    poni : str ; name of the *.poni file from Dioptas
    mask : 2d array [int] ; mask to be used on the tiff everythin except 0 is ignored 

    returns: np.array([theta, I]) ; 2 by n array with the theta angle and the corresponding intensity
    '''
    data_tiff = load_tiff(run, quad, dir)
    ai = pyFAI.load(poni)
    result1d = ai.integrate1d(data_tiff,
                        npt=1000,
                        correctSolidAngle=True,
                        mask=mask)
    return np.array(result1d)

def integrate2d(run : int, quad : int, dir : str, poni : str, mask : np.array):
    '''
    run : int ; number of the run to load
    quad : [0,1,2,3] ; the quad to be laoded
    dir : str ; datra directory
    poni : str ; name of the *.poni file from Dioptas
    mask : 2d array [int] ; mask to be used on the tiff everythin except 0 is ignored 

    returns: np.array ; 2D array with the cake of the corresponding detector
    '''
    data_tiff = load_tiff(run, quad, dir)
    ai = pyFAI.load(poni)
    result2d = ai.integrate2d(run,
                        npt_rad=1000,
                        npt_azim=1000,
                        correctSolidAngle=True,
                        mask=mask)
    return np.array(result2d)