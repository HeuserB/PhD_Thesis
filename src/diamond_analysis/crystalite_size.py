from ipaddress import summarize_address_range
from pyexpat.errors import XML_ERROR_INVALID_TOKEN
import wave
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
    print(f'Thetas are: {theta1} and {theta2}')
    theta1_err, theta2_err          = fit.params["C111_center"].stderr/2., fit.params["C220_center"].stderr/2.
    print(f'Theta errors are: {theta1_err} and {theta2_err}')

    size_1                          = K * wavelength / (beta1 * np.cos(np.deg2rad(theta1)))
    size_2                          = K * wavelength / (beta1 * np.cos(np.deg2rad(theta2)))
    sizes                           = np.array([size_1,size_2]) / 2.
    print(f"Crystallite size of is: %.2f and %.2f " %(sizes[0]*1e9,sizes[1]*1e9))
    return sizes


def crystalite_size_advanced(fit : ModelResult, wavelength : float, K : float):
    ### FWHM of Lorentzian peak is given by 2*sigma
    beta1, beta2                    = np.deg2rad(fit.params["C111_sigma"].value), np.deg2rad(fit.params["C220_sigma"].value)
    beta1_err, beta2_err            = np.deg2rad(fit.params["C111_sigma"].stderr), np.deg2rad(fit.params["C220_sigma"].stderr)
    y_values                        = np.array([np.log(beta1),np.log(beta2)])
    y_err                           = np.array([np.log(beta1_err),np.log(beta2_err)])
    
    theta1, theta2                  = fit.params["C111_center"].value/2., fit.params["C220_center"].value/2.
    print(f'Thetas are: {theta1} and {theta2}')
    theta1_err, theta2_err          = fit.params["C111_center"].stderr/2., fit.params["C220_center"].stderr/2.
    print(f'Theta errors are: {theta1_err} and {theta2_err}')
    x_values                        = np.array([1/np.cos(np.radians(theta1)), 1/np.cos(np.radians(theta2))])
    x_err                           = np.array([1/np.cos(theta1_err), 1/np.cos(theta2_err)])

    ### using lmfit to fit a linear model to the data points 
    line                            = LinearModel()
    pars = line.guess(y_values, x=x_values)
    pars.update(line.make_params())

    init = line.eval(pars, x=x_values)
    out = line.fit(y_values, pars, x=x_values, weights=1/y_err, calc_covar=True)

    fig, ax                         = plt.subplots()
    ax.scatter(x_values, y_values)

    intercept                       = out.params["intercept"]
    size                            = K * wavelength * np.exp(- intercept) / 2.


    size_err_1                      = np.abs(K * wavelength * np.sin(theta1) * theta1_err / \
                                        (2.3548 * beta1 * np.cos(2 * theta1) + 1 ) 
                                        ) + \
                                      np.abs(K * wavelength * theta1_err / \
                                        (2.3548 * beta1**2 * np.cos(theta1) )  
                                        )
    size_err_2                      = np.abs(K * wavelength * np.sin(theta2) * theta2_err / \
                                        (2.3548 * beta2 * np.cos(2 * theta2) + 1 ) 
                                        ) + \
                                      np.abs(K * wavelength * theta1_err / \
                                        (2.3548 * beta2**2 * np.cos(theta2) )  
                                        )
    size_err_1                      = np.sqrt(np.sum(size_err_1**2 + size_err_2**2))

    print(f"Crystallite size of is: %.2f +/- %.2f " %(size*1e9,size_err_1*1e9))