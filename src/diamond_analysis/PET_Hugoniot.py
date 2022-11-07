import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

Hug_PET             =    np.loadtxt('../../.data_LW03/Hugoniot_PET/hugoniot_mylar_dftmd.txt',skiprows=11)

us                  =    np.sqrt( 1 / (1 - 1/Hug_PET[:,1]) * (Hug_PET[:,4] + 1.5106100E-02) / (1.3800000) ) * 10.
up                  =    us * (1 - 1/Hug_PET[:,1])

interpolationP     = interpolate.interp1d(us, Hug_PET[:,4]*100.)
interpolationT     = interpolate.interp1d(us, Hug_PET[:,2])

def us_to_P(us):
    """
    Convert us [km/s] to P [GPa]
  
    From interpolated Mylar Hugoniot DFT data 

    Parameters:
    us (float): shock velocity in km/s
  
    Returns:
    float: Pressure in GPa
  
    """
    return interpolationP(us)

def us_to_T(us):
    """
    Convert us [km/s] to T [K]
  
    From interpolated Mylar Hugoniot DFT data 

    Parameters:
    us (float): shock velocity in km/s
  
    Returns:
    float: Temperature in K
  
    """
    return interpolationT(us)