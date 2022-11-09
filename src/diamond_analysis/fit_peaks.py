import numpy as np
from load_data import *
import matplotlib.pyplot as plt

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

def custom_lorentzian(x, C111_amplitude, C111_center, C111_sigma, C220_center, C220_sigma):
    return (C111_amplitude/np.pi * (C111_sigma/((x-C111_center)**2+C111_sigma**2))) +\
        (C111_amplitude * 0.288 /np.pi * (C220_sigma/((x-C220_center)**2+C220_sigma**2)))

def fit_diamond(data, plot=False,C111_max=3.):

    background = PolynomialModel(degree=7, prefix='bkg_')
    pars = background.guess(data[1], x=data[0])

    pars['bkg_c0'].set(value=np.random.normal(15000.,2000))
    pars['bkg_c1'].set(value=np.random.normal(-2500,2000))
    pars['bkg_c2'].set(value=np.random.normal(143,500))
    pars['bkg_c3'].set(value=np.random.normal(-4,5))
    pars['bkg_c4'].set(value=np.random.normal(0.06,0.005))
    pars['bkg_c5'].set(value=np.random.normal(-0.0005,5e-4))
    pars['bkg_c6'].set(value=np.random.normal(2.0e-06,5e-6))
    pars['bkg_c7'].set(value=np.random.normal(2.0e-06,5e-6),min=0.)

    '''
    diamond_peaks   =    Model(custom_lorentzian, prefix='peaks_')
    pars.update(diamond_peaks.make_params())
    
    pars['peaks_C111_center'].set(value=38.443, min=30, max=40)
    pars['peaks_C111_sigma'].set(value=np.random.normal(0.22,0.5),min=0.,max=C111_max)
    pars['peaks_C111_amplitude'].set(value=np.random.normal(3143,50))

    pars['peaks_C220_center'].set(value=63.9, min=60, max=67)
    pars['peaks_C220_sigma'].set(value=2.5, min=0.,max=3.7)

    '''
    peak1 = LorentzianModel(prefix='C111_')
    pars.update(peak1.make_params())

    pars['C111_center'].set(value=38.45, min=30, max=40)
    pars['C111_sigma'].set(value=np.random.normal(0.22,0.5),min=0.,max=C111_max)
    pars['C111_amplitude'].set(value=np.random.normal(3143,50))

    peak2 = LorentzianModel(prefix='C220_')
    pars.update(peak2.make_params())

    pars['C220_center'].set(value=63.9, min=60, max=67)
    pars['C220_sigma'].set(value=2.5, max=3.7)
    #pars['C220_amplitude'].set(value=np.random.normal(1000,50), min=0.)
    pars['C220_amplitude'].set(value=np.random.normal(1000,50), min=0.)
    

    PET_peak = LorentzianModel(prefix='PET_')
    pars.update(PET_peak.make_params())

    pars['PET_center'].set(value=27, min=26, max=29.)
    pars['PET_sigma'].set(value=1.5, max=2.5)
    #pars['PET_amplitude'].set(value=np.random.normal(1000,50), min=0.)
    pars['PET_amplitude'].set(value=np.random.normal(1000,50), min=0.,max=0.1)

    model = peak1 + peak2 + background + PET_peak
    #model = diamond_peaks + background + PET_peak

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

def sampleChi2(runs=[182,186,188,190,192],rounds=10,C111_max=3.):
    ### run a routine to fit and only store if the Chi2 of the fit is better ###
    fits        = [[]] * len(runs)

    # repeat fitting with random inital parameters and overwrite if the fit is better 
    for i in range(rounds):
        for run_id, run in enumerate(runs):
                data                = np.loadtxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy").T
                fits[run_id]        = load_modelresult( f'../../.data_LW03/fits/r{run}__modelresult.sav')
                new_fit, fig        = fit_diamond(data, C111_max=C111_max)
                if new_fit.chisqr < fits[run_id].chisqr:
                    if new_fit.params["C111_center"].stderr is None:
                    #if new_fit.params["peaks_C111_center"].stderr is None:
                        print(f"Fit {i} for run {run} yielded no uncertanties and will therefore not be stored!")
                        continue
                    old_chi      = fits[run_id].chisqr
                    fits[run_id] = new_fit
                    save_modelresult(fits[run_id], f'../../.data_LW03/fits/r{run}__modelresult.sav')
                    print(f"Old Chi2 is {old_chi} new Chi2 is {new_fit.chisqr}")

def sampleChi2_SACLA(runs,rounds=10,C111_max=100.,modeltype="individual",include_PET=False,resample=False):
    ### run a routine to fit and only store if the Chi2 of the fit is better ###
    fits        = [[]] * len(runs)

    # repeat fitting with random inital parameters and overwrite if the fit is better 
    for i in range(rounds):
        for run_id, run in enumerate(runs):
                print(f'Fitting run {run}.')
                data,_                  = load_SACLA(run)
                data                    = data.T
                fits[run_id]            = load_modelresult( f'../../.data_SACLA/fits/r{run}__modelresult.sav')
                if resample == True:
                    fit = fits[run_id]
                    try:
                        beta1_err, beta2_err, beta3_err            = np.deg2rad(fit.params["C111_sigma"].stderr), \
                                                                    np.deg2rad(fit.params["C220_sigma"].stderr), \
                                                                    np.deg2rad(fit.params["C311_sigma"].stderr)
                        theta1_err, theta2_err,theta3_err          = fit.params["C111_center"].stderr/2., \
                                                                    fit.params["C220_center"].stderr/2.,\
                                                                    fit.params["C220_center"].stderr/2.
                        print(f'Skipping run {run} because it already has uncertainties stored.')
                        continue
                    except:
                        pass
                try:
                    new_fit, fig            = fit_diamond_SACLA(data, C111_max=C111_max,modeltype=modeltype,include_PET=include_PET)
                except:
                    print(f"Fit {i} failed")
                    continue
                if new_fit.chisqr < fits[run_id].chisqr:
                    if new_fit.params["C111_center"].stderr is None:
                    #if new_fit.params["peaks_C111_center"].stderr is None:
                        print(f"Fit {i} for run {run} yielded no uncertanties and will therefore not be stored!")
                        print(new_fit.fit_report())
                        continue
                    old_chi             = fits[run_id].chisqr
                    fits[run_id]        = new_fit
                    save_modelresult(fits[run_id], f'../../.data_SACLA/fits/r{run}__modelresult.sav')
                    print(f"Old Chi2 is {old_chi} new Chi2 is {new_fit.chisqr}")

def init_fits(runs=[182,186,188,190,192],C111_max=3.):
    ### fit all runs and store the optimied models in files ###
    fits        = [[]] * len(runs)
    
    for run_id, run in enumerate(runs):
        data                    = np.loadtxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy").T
        fits[run_id], fig       = fit_diamond(data,C111_max=C111_max)
        #fig.savefig(f'../../.data_LW03/figures/r{run}__fit_plot.svg')
        save_modelresult(fits[run_id], f'../../.data_LW03/fits/r{run}__modelresult.sav')

def init_fits_SACLA(runs,C111_max=100.,modeltype="individual",include_PET=False):
    ### fit all runs and store the optimied models in files ###
    fits        = [[]] * len(runs)
    
    for run_id, run in enumerate(runs):
        data,_                  = load_SACLA(run)
        data                    = data.T
        failed                  = True
        while failed == True:
            try:
                fits[run_id], fig       = fit_diamond_SACLA(data,C111_max=C111_max,modeltype=modeltype,include_PET=include_PET)
                failed = False 
            except:
                pass
        #fig.savefig(f'../../.data_LW03/figures/r{run}__fit_plot.svg')
        save_modelresult(fits[run_id], f'../../.data_SACLA/fits/r{run}__modelresult.sav')

def lorentzian(x, center, sigma):
    return 1/(np.pi) * ( (sigma) /( (x - center)**2 + (sigma)**2 ))

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

def fit_diamond_SACLA(data, plot=False,C111_max=100., modeltype="individual",include_PET=False):
    def background(x,a,b,c,d):
        return a * (x-d)**2 + b*(x-d) + c

    def liquid_peak(x,amp,sigma,mu):
        return amp / (np.sqrt(2*np.pi * sigma) ) * np.exp( - 1 *  (x - mu)**2 / (2 * sigma**2))

    bkg_model               = Model(background,prefix='bkg_')
    liquid_model            = Model(liquid_peak,prefix='liquid_')

    pars                    = bkg_model.make_params()
    pars.update(liquid_model.make_params())

    pars['bkg_a'].set(value=np.random.normal(0.21708,0.1))
    pars['bkg_b'].set(value=np.random.normal(-3.031587,1.5))
    pars['bkg_c'].set(value=np.random.normal(500.18140852,150))
    pars['bkg_d'].set(value=np.random.normal(39.,5.0))

    pars['liquid_amp'].set(value=np.random.normal(400,100.))
    pars['liquid_sigma'].set(value=np.random.normal(7.,1.))
    pars['liquid_mu'].set(value=np.random.normal(38.,2.))

    model                   = bkg_model + liquid_model

    if include_PET == True:
        print('Including PET peak')
        PET_peak = LorentzianModel(prefix='PET_')
        pars.update(PET_peak.make_params())

        pars['PET_center'].set(value=42.)
        pars['PET_sigma'].set(value=1.5)
        pars['PET_amplitude'].set(value=250)

        model += PET_peak

    if modeltype=='individual':
        peak1 = LorentzianModel(prefix='C111_')
        pars.update(peak1.make_params())

        peak2 = LorentzianModel(prefix='C220_')
        pars.update(peak2.make_params())

        peak3 = LorentzianModel(prefix='C311_')
        pars.update(peak3.make_params())

        pars['C111_center'].set(35.6, min=35.0482, max=40)
        pars['C220_center'].set(58.4, min=58.906, max=67)
        pars['C311_center'].set(71.9, min=70.4201, max=75)
        pars['C111_sigma'].set(2, min=0, max=100)
        pars['C220_sigma'].set(2, min=0, max=50)
        pars['C311_sigma'].set(7, min=0, max=50)

        model               += peak1 + peak2 + peak3
    
    if modeltype=='combined':
        def triple_peak(x, amp, C111_center, C111_sigma, C220_center, C220_sigma, C311_center, C311_sigma, scale220, scale311):
            return amp * ( lorentzian(x,C111_center,C111_sigma) + \
                scale220 * 0.46 * lorentzian(x,C220_center,C220_sigma) + \
                scale311 * 0.081 * lorentzian(x,C311_center,C311_sigma) )

        triple_peak         = Model(triple_peak)
        pars.update(triple_peak.make_params())

        pars['C111_center'].set(value=np.random.normal(34.3,1.), min=30, max=40)
        pars['C220_center'].set(value=np.random.normal(58,1), min=55, max=67)
        pars['C311_center'].set(value=np.random.normal(69.32,1.5), min=68, max=72)
        pars['C111_sigma'].set(value=np.random.normal(2,0.1), min=0, max=100)
        pars['C220_sigma'].set(value=np.random.normal(2,0.1), min=0, max=50)
        pars['C311_sigma'].set(value=np.random.normal(7,0.25), min=0, max=50)
        pars['amp'].set(value=np.random.normal(250,10),min=0.,max=10000.)
        pars['scale220'].set(value=np.random.normal(1.1,0.05), min=0.4, max=1.3)
        pars['scale311'].set(value=np.random.normal(1.,0.05), min=0.3, max=1.3)
        #pars['scale311'].set(value=np.random.normal(1.,0.05), min=0.7, max=1.3)

        model               += triple_peak
    
    out = model.fit(data[1], pars, x=data[0],weights=1/data[2])

    if plot==True:

        fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
        axes[0].plot(data[0], data[1],label='XRD data')
        axes[0].plot(data[0], out.best_fit, '-', label='best fit')
        axes[0].set_ylabel(r'$I$')
        axes[0].set_xlabel(r'$2\theta$')
        axes[0].legend()

        comps = out.eval_components(x=data[0])

        axes[1].plot(data[0], data[1])
        axes[1].plot(data[0], comps['C111_'], '--', label='Peak 111')
        axes[1].plot(data[0], comps['C220_'], '--', label='Peak 220')
        axes[1].plot(data[0], comps['bkg_'], '--', label='Background')
        axes[1].plot(data[0], comps['PET_'], '--', label='PET peak')
        axes[1].set_xlabel(r'$2\theta$')
        axes[1].legend()

        plt.show() 

    else:
        fig = None
    return out, fig