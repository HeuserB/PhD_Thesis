import numpy as np
from load_data import *

from lmfit.models import GaussianModel, PolynomialModel, LorentzianModel, ExponentialModel, LinearModel
from lmfit import Parameters, fit_report, minimize, Model
from lmfit.model import save_modelresult, load_modelresult

NA = np.newaxis

def gauss(x, mu, sigma, C):
    return C / np.sqrt(2 * np.pi) / np.abs(sigma) * np.exp(-(x - mu)**2 / 2.0 / sigma**2)

def lingaussian(x,slope,intercept,amp_gauss, cen_gauss, wid_gauss):
#    "1-d gaussian: lingaussian(x, slope,intercept,amp, cen, wid)"
    return (slope*x+intercept+(amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2)))

def lin_gauss_lorentzian(x, slope,intercept,amp_gauss, cen_gauss, wid_gauss,amp_lor, cen_lor, wid_lor):
#    "1-d lorentzian: lorentzian(x, amp, cen, wid)"
    return ((slope*x+intercept)\
              +(amp_gauss/(np.sqrt(2*np.pi)*wid_gauss)) * np.exp(-(x-cen_gauss)**2 /(2*wid_gauss**2))\
              +(amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2))))

def lorentzian(x,amp_lor,cen_lor,wid_lor):
    return (amp_lor/np.pi * (wid_lor/((x-cen_lor)**2+wid_lor**2)))

def fit_diamond(data, plot=False):

    background = PolynomialModel(degree=7, prefix='bkg_')
    pars = background.guess(data[1], x=data[0])

    pars['bkg_c0'].set(value=np.random.normal(15000.,1000))
    pars['bkg_c1'].set(value=np.random.normal(-2500,1000))
    pars['bkg_c2'].set(value=np.random.normal(143,500))
    pars['bkg_c3'].set(value=np.random.normal(-4,5))
    pars['bkg_c4'].set(value=np.random.normal(0.06,0.005))
    pars['bkg_c5'].set(value=np.random.normal(-0.0005,5e-4))
    pars['bkg_c6'].set(value=np.random.normal(2.0e-06,5e-6))
    pars['bkg_c7'].set(value=np.random.normal(2.0e-06,5e-6),min=0.)

    peak1 = LorentzianModel(prefix='C111_')
    pars.update(peak1.make_params())

    pars['C111_center'].set(value=39, min=30, max=40)
    pars['C111_sigma'].set(value=np.random.normal(2.22,0.5))
    pars['C111_amplitude'].set(value=np.random.normal(3143,50))

    peak2 = LorentzianModel(prefix='C220_')
    pars.update(peak2.make_params())

    pars['C220_center'].set(value=63.9, min=60, max=67)
    pars['C220_sigma'].set(value=2.5, max=3.7)
    pars['C220_amplitude'].set(value=np.random.normal(1000,50), min=0.)

    PET_peak = LorentzianModel(prefix='PET_')
    pars.update(PET_peak.make_params())

    pars['PET_center'].set(value=27, min=26, max=29.)
    pars['PET_sigma'].set(value=1.5, max=2.5)
    pars['PET_amplitude'].set(value=np.random.normal(1000,50), min=0.)

    model = peak1 + peak2 + background + PET_peak

    init = model.eval(pars, x=data[0])
    out = model.fit(data[1], pars, x=data[0])
    #print(out.fit_report())

    # print(out.fit_report)
    if plot==True:

        fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
        axes[0].plot(data[0], data[1],label='XRD data')
        #axes[0].plot(data[0], init, '--', label='initial fit')
        axes[0].plot(data[0], out.best_fit, '-', label='best fit')
        axes[0].set_ylabel(r'$I$')
        axes[0].set_xlabel(r'$2\theta$')
        axes[0].legend()

        comps = out.eval_components(x=data[0])
        #print(out.params.valuesdict())

        axes[1].plot(data[0], data[1])
        axes[1].plot(data[0], comps['C111_'], '--', label='Peak 111')
        axes[1].plot(data[0], comps['C220_'], '--', label='Peak 220')
        axes[1].plot(data[0], comps['bkg_'], '--', label='Background')
        axes[1].plot(data[0], comps['PET_'], '--', label='PET peak')
        #axes[1].set_ylabel(r'$I [au]$')
        axes[1].set_xlabel(r'$2\theta$')
        axes[1].legend()

        plt.show() 

    else:
        fig = None
    return out, fig

def sampleChi2(runs=[182,186,188,190,192],rounds=10):
    ### run a routine to fit and only store if the Chi2 of the fit is better ###
    fits        = [[]] * len(runs)

    # repeat fitting with random inital parameters and overwrite if the fit is better 
    for i in range(rounds):
        for run_id, run in enumerate(runs):
                data                = np.loadtxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy").T
                fits[run_id]        = load_modelresult( f'../../.data_LW03/fits/r{run}__modelresult.sav')
                new_fit, fig        = fit_diamond(data)
                if new_fit.chisqr < fits[run_id].chisqr:
                    if new_fit.params["C111_center"].stderr is None:
                        print(f"Fit {i} for run {run} yielded no uncertanties and will therefore not be stored!")
                        continue
                    old_chi      = fits[run_id].chisqr
                    fits[run_id] = new_fit
                    save_modelresult(fits[run_id], f'../../.data_LW03/fits/r{run}__modelresult.sav')
                    print(f"Old Chi2 is {old_chi} new Chi2 is {new_fit.chisqr}")

def init_fits(runs=[182,186,188,190,192]):
    ### fit all runs and store the optimied models in files ###
    fits        = [[]] * len(runs)
    
    for run_id, run in enumerate(runs):
        data                    = np.loadtxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy").T
        fits[run_id], fig       = fit_diamond(data)
        #fig.savefig(f'../../.data_LW03/figures/r{run}__fit_plot.svg')
        save_modelresult(fits[run_id], f'../../.data_LW03/fits/r{run}__modelresult.sav')


'''
def fit_double_peak(data):
    twoTheta, I = data[0], data[1]
    fit_params          = Parameters()
    fit_params.add('scale', value=0.01)
    fit_params.add('intercept', value=200.)
    fit_params.add('amp_gauss', value=1.)
    fit_params.add('cen_gauss', value=30.9)
    fit_params.add('wid_gauss', value=0.9)
    fit_params.add('amp_lor', value=0.9)
    fit_params.add('cen_lor', value=37.)
    fit_params.add('wid_lor', value=0.9)
    fit_params.add('amp_lor2', value=66.)
    fit_params.add('cen_lor2', value=0.9)
    fit_params.add('wid_lor2', value=0.9)
    

    out                 = minimize(residual, fit_params, args=(twoTheta,), kws={'data': I});
    param_dict          = out.params.valuesdict()
    return out, param_dict


def residual(pars, x, data=None, eps=None):
    # unpack parameters: extract .value attribute for each parameter
    parvals     = pars.valuesdict()
    slope       = parvals['scale']
    intercept   = parvals['intercept']
    amp_gauss   = parvals['amp_gauss']
    cen_gauss   = parvals['cen_gauss']
    wid_gauss   = parvals['wid_gauss']
    amp_lor     = parvals['amp_lor'] 
    cen_lor     = parvals['cen_lor']
    wid_lor     = parvals['wid_lor']
    amp_lor2    = parvals['amp_lor2']
    cen_lor2    = parvals['cen_lor2']
    wid_lor2    = parvals['wid_lor2']

    model = lin_gauss_2lorentzian(x, slope,intercept,amp_gauss, cen_gauss, wid_gauss, amp_lor, cen_lor, wid_lor, amp_lor2, cen_lor2, wid_lor2)

    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data) / eps

'''

def fit_background(data):

    background = PolynomialModel(degree=7, prefix='bkg_')
    pars = background.guess(data[1], x=data[0])

    pars['bkg_c0'].set(value=np.random.normal(15000.,1000))
    pars['bkg_c1'].set(value=np.random.normal(-2500,1000))
    pars['bkg_c2'].set(value=np.random.normal(143,500))
    pars['bkg_c3'].set(value=np.random.normal(-4,5))
    pars['bkg_c4'].set(value=np.random.normal(0.06,0.005))
    pars['bkg_c5'].set(value=np.random.normal(-0.0005,5e-4))
    pars['bkg_c6'].set(value=np.random.normal(2.0e-06,5e-6))
    pars['bkg_c7'].set(value=np.random.normal(2.0e-06,5e-6),min=0.)

    model   = background

    init    = model.eval(pars, x=data[0])
    out     = model.fit(data[1], pars, x=data[0])
    #print(out.fit_report())

    return out, None