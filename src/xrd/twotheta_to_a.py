import numpy as np

def twotheta_to_a(twotheta : float, twotheta_err : float, peak : np.array, wavelength=1.240343e-10):
    a           = wavelength / (2 * np.sin(np.deg2rad(twotheta/2.))) * np.sqrt(np.sum(peak**2)) * 1e10
    da          = wavelength / 2 * np.sqrt(np.sum(peak**2)) * np.abs(- (np.cos(np.deg2rad(twotheta/2.)))/(np.sin(np.deg2rad(twotheta/2.))**2) ) * twotheta_err  * 1e10 * 0.5
    return a, da