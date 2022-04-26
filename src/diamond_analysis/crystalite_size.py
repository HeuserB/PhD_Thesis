from pyexpat.errors import XML_ERROR_INVALID_TOKEN
import numpy as np

from lmfit.models import LinearModel
from lmfit import Parameters, fit_report, minimize, Model
from lmfit.model import save_modelresult, load_modelresult, ModelResult

import matplotlib.pyplot as plt

def crystalite_size(fit : ModelResult, wavelength : float, K : float):
    ### FWHM of Lorentzian peak is given by 2*sigma
    beta1, beta2                    = np.deg2rad(fit.params["C111_sigma"].value), np.deg2rad(fit.params["C220_sigma"].value)
    beta1_err, beta2_err            = np.deg2rad(fit.params["C111_sigma"].stderr), np.deg2rad(fit.params["C220_sigma"].stderr)
    y_values                        = np.array([np.log(beta1),np.log(beta2)])
    y_err                           = np.array([np.log(beta1_err),np.log(beta2_err)])
    
    theta1, theta2                  = fit.params["C111_center"].value/2., fit.params["C220_center"].value/2.
    theta1_err, theta2_err          = fit.params["C111_center"].stderr/2., fit.params["C220_center"].stderr/2.
    x_values                        = np.array([1/np.cos(theta1), 1/np.cos(theta2)])
    x_err                           = np.array([1/np.cos(theta1_err), 1/np.cos(theta2_err)])

    ### using lmfit to fit a linear model to the data points 
    line                            = LinearModel()
    pars = line.guess(y_values, x=x_values)
    pars.update(line.make_params())

    init = line.eval(pars, x=x_values)
    out = line.fit(y_values, pars, x=x_values)
    print(out.fit_report())

    fig, ax                         = plt.subplots()
    ax.scatter(x_values, y_values)