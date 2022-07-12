import sys

import re
from lmfit.models import GaussianModel, PolynomialModel, LorentzianModel, ExponentialModel, LinearModel
from lmfit import Parameters, fit_report, minimize, Model
from lmfit.model import save_modelresult, load_modelresult

import pandas as pd

import pyFAI
import numpy as np
import imageio
import matplotlib.pyplot as plt
sys.path.append("../xrd")

from load_data import load_tiff
from combine_quads import merge_four_quads
from peak_shapes import *

import logging
from stipcrawl import setup_logger
logger = setup_logger.logger
quad_scale = np.loadtxt("../../.data_LW03/instprm/quad_scale_LW03.txt")

wavelength = 1.3051 * 1e-10

# define the data directory
dir_data = "../../.data_LW03"

# the logbook file
logbook = "../../.data_LW03/logbook_LW03_subset.ods"

poni_file_q0 = "../../.data_LW03/poni/Quad0_MS.poni"
poni_file_q1 = "../../.data_LW03/poni/Quad1_MS.poni"
poni_file_q2 = "../../.data_LW03/poni/Quad2_MS.poni"
poni_file_q3 = "../../.data_LW03/poni/Quad3_MS.poni"

# the mask files for the quad detectors
mask_q0 = "../../.data_LW03/mask/Quad0_ext.mask"
mask_q1 = "../../.data_LW03/mask/Quad1_ext.mask"
mask_q2 = "../../.data_LW03/mask/Quad2.mask"
mask_q3 = "../../.data_LW03/mask/Quad3.mask"

# define the pattern for a background file
pattern_runs = re.compile("r[0-9]*_bkgCorrected")

# define the azimuthal intyegrator objects from pyFAI
q0_ai = pyFAI.load(poni_file_q0)
q1_ai = pyFAI.load(poni_file_q1)
q2_ai = pyFAI.load(poni_file_q2)
q3_ai = pyFAI.load(poni_file_q3)

# load the masks as np arrays

q0_mask = np.array(imageio.imread(mask_q0),dtype=int)
q1_mask = np.array(imageio.imread(mask_q1),dtype=int)
q2_mask = np.array(imageio.imread(mask_q2),dtype=int)
q3_mask = np.array(imageio.imread(mask_q3),dtype=int)

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
    return merged_xy
    #np.savetxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy",merged_xy.T)

def load_spectrum(run : int, dir_data : str = '/mnt/TOSHIBA_EXT/HZDR_data/Data/2020_LCLS-Kraus/Data/Detectors/ePix10k', q2_ai = q2_ai, q3_ai = q3_ai, q2_mask= q2_mask, q3_mask = q3_mask):
    merged_xy = merge_quad23(run, dir_data, q2_ai, q3_ai,q2_mask,q3_mask)
    return merged_xy