from scipy import interpolate
import numpy as np

data                    = np.loadtxt('../../.data/expansion_coefficient_diamond.txt')
a_to_T                  = interpolate.interp1d(data[:,0],data[:,1])

def T_from_a(a, lower_limit=3.5678):
    if  isinstance(a, np.ndarray):
        result          = np.zeros_like(a)
        result[np.where(a > lower_limit)] = a_to_T(a[np.where(a > lower_limit)])
    else:
        if a <= lower_limit:
            result = 0.
        else:
            result = a_to_T(a)
    return result

def twotheta_to_a(twotheta, twotheta_err, peak, wavelength=1.240343):
    a           = wavelength / (2 * np.sin(np.deg2rad(twotheta/2.))) * np.sqrt(np.sum(peak**2))
    da          = np.abs((wavelength / 2 * np.sqrt(np.sum(peak**2)) * (-1) * (np.cos(np.deg2rad(twotheta/2.)))/(np.sin(np.deg2rad(twotheta/2.))**2))**2) * (np.deg2rad(twotheta_err * 0.5))**2
    if  isinstance(twotheta, np.ndarray):
        return np.hstack([a[:,np.newaxis], np.sqrt(da)[:,np.newaxis]])
    else:
        return a, np.sqrt(da)

def T_from_peak(twotheta, twotheta_err, peak, wavelength=1.240343):
    a                   = wavelength / (2 * np.sin(np.deg2rad(twotheta/2.))) * np.sqrt(np.sum(peak**2))
    a_min               = wavelength / (2 * np.sin(np.deg2rad((twotheta+twotheta_err)/2.))) * np.sqrt(np.sum(peak**2))
    a_max               = wavelength / (2 * np.sin(np.deg2rad((twotheta-twotheta_err)/2.))) * np.sqrt(np.sum(peak**2))
    T, Tmin, Tmax       = T_from_a(a), T_from_a(a_min), T_from_a(a_max)
    V, Vmin, Vmax       = a**3, a_min**3, a_max**3
    return T, Tmin, Tmax, V, Vmin, Vmax

def fit_peak(shotfile="../../.data_SACLA/logbook/PET.pkl",run_filter=[424,436,422,434,420,432,418,430,412,426,414,617],thetarange=[31., 36.5],peak=np.array([1,1,1])):
    if ((peak == np.array([1,1,1])).all()):
        thetarange=[31., 36.5]
    elif ((peak == np.array([2,2,0])).all()):    
        thetarange=[55., 62.0]
    print("Fitting all runs and storing the results")
    from lmfit.model import Model, save_modelresult
    import pandas as pd
    import sys
    from peak_shapes import background, lorentzian
    
    sys.path.append("../diamond_analysis")
    from load_data import load_SACLA

    SACLA_shots         =   pd.read_pickle(shotfile)
    PET75_shots         =   SACLA_shots[SACLA_shots['target'] == 'PET75']
    PET75_shots         =   PET75_shots[PET75_shots['run'].isin(run_filter)]
    delays_SACLA            = np.array([PET75_shots.loc[PET75_shots['run'] == run]['delay'].values[0] for run in run_filter])
    SACLA_dict              = {}
    for delay, run in zip(delays_SACLA,run_filter):
        SACLA_dict.update({delay: run})

    for delay_id, delay in enumerate([8,9,10,11,12,13,14,15,16,17,18,19]):
        run                         = SACLA_dict[delay]
        run_data, ref_data          =   load_SACLA(run)

        thetamin, thetamax          = thetarange[0], thetarange[1]
        mask                        = np.logical_and(run_data[:,0] > thetamin, run_data[:,0] < thetamax)

        def background(x,a,b):
            return a * x + b

        def single_peak(x, amp, C_center, C_sigma):
            return lorentzian(x,amp,C_center,C_sigma)

        stderr       =    None
        while stderr==None:

            bkg_model                   = Model(background,prefix='bkg_')

            pars                        = bkg_model.make_params()

            pars['bkg_a'].set(value=np.random.normal(0.21708,0.1))
            pars['bkg_b'].set(value=np.random.normal(-3.031587,1.5))

            model                       = bkg_model
                
            peakmodel                   = Model(single_peak)
            pars.update(peakmodel.make_params())
            
            if ((peak == np.array([1,1,1])).all()):
                print("Fitting peak (111)!")
                pars['C_center'].set(value=np.random.normal(34.0,1.), min=30, max=40)
                pars['C_sigma'].set(value=np.random.normal(2,0.1), min=0, max=100)
                pars['amp'].set(value=np.random.normal(250,10),min=0.,max=10000.)
            
            elif (peak == np.array([2,2,0])).all():
                print("Fitting peak (220)!")
                pars['C_center'].set(value=np.random.normal(58.5,1.), min=57, max=60)
                pars['C_sigma'].set(value=np.random.normal(2,0.1), min=0, max=100)
                pars['amp'].set(value=np.random.normal(550,10),min=0.,max=10000.)

            else:
                print('Please enter a valid peak!')
                break
            model                       += peakmodel

            out                         = model.fit(run_data[:,1][mask], pars, x=run_data[:,0][mask],weights=1/run_data[:,2][mask])
            stderr                      = out.params["C_center"].stderr
            if stderr==None:
                continue
            elif np.isnan(float(stderr)):
                print(run)
                stderr=None



        str_peak    = ''
        for ikl in peak:
            str_peak += str(ikl)

        save_modelresult(out, f'../../.data_SACLA/fits/C{str_peak}_peak_r{run}_modelresult.sav')

if __name__=='__main__':
    fit_peak(peak=np.array([1,1,1]))
    fit_peak(peak=np.array([2,2,0]))