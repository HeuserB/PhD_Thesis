import matplotlib.pyplot as plt
from lmfit.model import load_modelresult
from load_data import load_SACLA
from fit_peaks import lorentzian
import numpy as np
import pandas as pd
import seaborn as sns

def plotBestFits_SACLA(runs, df):
    fits            = [[]] * len(runs)
    fig             = plt.figure(figsize=(12, 8))
    gs              = fig.add_gridspec(len(runs)*3, 2)

    for run_id, run in enumerate(runs):
        data,_              = load_SACLA(run)
        data                = data.T
        fits[run_id]        = load_modelresult( f'../../.data_SACLA/fits/r{run}__modelresult.sav')
        params              = fits[run_id].params
        t                   = df.loc[df['run'] == run]['delay'].values[0]

        #t_id                = (int(t-np.min(list(probingTime.values()))))
        t_id                = run_id
        comps               = fits[run_id].eval_components(x=data[0])
        
        ax1                 = fig.add_subplot(gs[t_id*3:t_id*3+2, 0:0+1])

        ax1.plot(data[0], data[1],label=f'XRD data run {run}')
        
        ax1.plot(data[0], fits[run_id].best_fit, '-', label='best fit')
        ax1.set_ylabel(r'$I$ [\textit{au}]',fontsize=20)

        if run_id == len(runs)-1:
            ax1.set_xlabel(r'$2\theta$ [\textit{deg}]',fontsize=20)
            
        ax1.legend()
        ax1.axes.yaxis.set_ticks([])
        
        ax2                 = fig.add_subplot(gs[t_id*3:t_id*3+3, 1:1+1])
        ax2.plot(data[0], data[1])

        y_111               = params['amp'] * lorentzian(data[0],params['C111_center'],params['C111_sigma'])
        y_220               = params['amp'] * lorentzian(data[0],params['C220_center'],params['C220_sigma']) * params['scale220'] * 0.46
        y_311               = params['amp'] * lorentzian(data[0],params['C311_center'],params['C311_sigma']) * params['scale311'] * 0.081

        try:
            ax2.plot(data[0], comps['PET_'], '--', label='PET peak')
        except:
            print()
            
        ax2.plot(data[0], y_111, '--', label='Peak 111')
        ax2.plot(data[0], y_220, '--', label='Peak 220')
        ax2.plot(data[0], y_311, '--', label='Peak 311')

        ax2.plot(data[0], comps['bkg_'] + comps['liquid_'], '--', label='Background + Liquid')
        ax1.axes.axvline(x=params['C111_center'].value, label='Diamond 111',linestyle='dashed')
        ax1.axes.axvline(x=params['C220_center'].value, label='Diamond 220',linestyle='dashed')
        ax1.axes.axvline(x=params['C311_center'].value, label='Diamond 311',linestyle='dashed')


        ax2.axes.yaxis.set_ticks([])
        
        if run_id == len(runs)-1:
            ax2.set_xlabel(r'$2\theta$ [\textit{deg}]',fontsize=20)
        ax2.legend()

        ax3                 = fig.add_subplot(gs[t_id*3+2:t_id*3+3, 0:0+1])
        residuals           = (fits[run_id].best_fit - data[1])/np.max(data[1])
        ax3.plot(data[0],residuals)
    fig.savefig(f'../../.data_SACLA/figures/bestFits.svg')
    plt.show() 

def plotBestFits(runs=[182,186,188,190,192,298,286]):
    fits            = [[]] * len(runs)
    fig             = plt.figure(figsize=(12, 40))
    gs              = fig.add_gridspec(len(runs)*3, 2)
    probingTime     = {182:9.,186:8.,188:7.,190:10.,192:11.,289:12,286:12}

    for run_id, run in enumerate(runs):
        data                = np.loadtxt(f"../../.data_LW03/lineouts/r{run}_Q23.xy").T
        fits[run_id]        = load_modelresult( f'../../.data_LW03/fits/r{run}__modelresult.sav')
        comps               = fits[run_id].params
        try:
            print(comps["C111_center"])
            print(comps["C111_sigma"])
            #print(comps)
        except:
            print(comps["peaks_C220_center"])
            print(comps["peaks_C220_sigma"])
        #print(comps["C111_amplitude"].value/comps["C220_amplitude"].value)

        t                   = int(probingTime[run])
        t_id                = (int(t-np.min(list(probingTime.values()))))
        comps               = fits[run_id].eval_components(x=data[0])
        #print(out.params.valuesdict())
        ax1                 = fig.add_subplot(gs[t_id*3:t_id*3+2, 0:0+1])

        ax1.plot(data[0], data[1],label=f'XRD data run {run}')
        #axes[0].plot(data[0], init, '--', label='initial fit')
        ax1.plot(data[0], fits[run_id].best_fit, '-', label='best fit')
        ax1.set_ylabel(r'$I$ [\textit{au}]',fontsize=20)
        if run_id == len(runs)-1:
            ax1.set_xlabel(r'$2\theta$ [\textit{deg}]',fontsize=20)
        #axes[t_id,0].plot(data[0], comps['bkg_'], '--', label='Background',color='r')
        ax1.legend()
        ax1.axes.yaxis.set_ticks([])
        
        ax2                 = fig.add_subplot(gs[t_id*3:t_id*3+3, 1:1+1])
        ax2.plot(data[0], data[1])
        try:
            ax2.plot(data[0], comps['peaks_'], '--', label='Peaks')
        except:
            ax2.plot(data[0], comps['C111_'], '--', label='Peak 111')
            ax2.plot(data[0], comps['C220_'], '--', label='Peak 220')
        ax2.plot(data[0], comps['bkg_'], '--', label='Background')
        ax2.plot(data[0], comps['PET_'], '--', label='PET peak')
        ax2.axes.yaxis.set_ticks([])
        #axes[1].set_ylabel(r'$I [au]$')
        if run_id == len(runs)-1:
            ax2.set_xlabel(r'$2\theta$ [\textit{deg}]',fontsize=20)
        ax2.legend()

        ax3                 = fig.add_subplot(gs[t_id*3+2:t_id*3+3, 0:0+1])
        residuals           = (fits[run_id].best_fit - data[1])/np.max(data[1])
        ax3.plot(data[0],residuals)
    fig.savefig(f'../../.data_LW03/figures/bestFits.svg')
    plt.show() 


def plot_waterfall(run_filter=[412,414,416,418,420,422,424,426,430,432,434,436,617,619], pickeldir="../../.data_SACLA/logbook/PET.pkl",target='PET75',outfile='../../.data_SACLA/figures/SACLA_waterfall.pdf'):
    all_shots           =   pd.read_pickle(pickeldir)
    subset_shots         =   all_shots[all_shots['target'] == target]
    subset_shots         =   subset_shots[subset_shots['run'].isin(run_filter)]

    timings             =   subset_shots['delay'].values
    timing_sort         =   np.sort(timings)

    fig, ax             =   plt.subplots(1,2,figsize=(17,7))

    for shot in range(len(subset_shots)):
        run                     =   subset_shots.iloc[shot]['run']
        E                       =   subset_shots.iloc[shot]['E_on_sample']
        delay                   =   subset_shots.iloc[shot]['delay']
        
        run_data, ref_data      =   load_SACLA(run)
        ### normalise Intensities
        run_data[:,1]           *=  1/max(run_data[:,1])
        ref_data[:,1]           *=  1/max(ref_data[:,1])

        offset                  =   int(np.where(timing_sort == delay)[0][0]) * 0.11

        ax[0].plot(run_data[:,0],run_data[:,1]-offset,label=f"run {run} at {E} J after {delay} ns")
        mask                    = np.logical_and(run_data[:,0] > 28.,run_data[:,0] < 44.)
        ax[1].plot(run_data[mask][:,0],run_data[mask][:,1]-offset,label=f"run {run} at {E} J after {delay} ns")

    ax[0].axvline(x=35.0540, label='Diamond 111',linestyle='dashed')
    ax[0].axvline(x=58.9159, label='Diamond 220',linestyle='dashed')
    ax[0].axvline(x=70.4332, label='Diamond 311',linestyle='dashed')

    ax[1].axvline(x=35.0540, label='Diamond 111',linestyle='dashed')

    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend([handles[idx] for idx in np.argsort(timings)],[labels[idx] for idx in np.argsort(timings)],loc='lower right')
    ax[0].set_xlabel(r"Scattering angle $2\Theta$")
    ax[0].set_ylabel(r"Intensity [au]")
    ax[1].set_xlabel(r"Scattering angle $2\Theta$")

    plt.savefig(outfile)

def plot_temperature_estimates():
    fig, ax_main                =   plt.subplots(1,2,figsize=(14,7))
    ax, ax_vol                  =   ax_main[0], ax_main[1]

    Tvals                   =   np.zeros([2,12])
    Terrors                 =   np.zeros([2,12])
    for peak, thetarange, color in zip([np.array([1,1,1]),np.array([2,2,0])],[[31., 36.5],[55., 62.0]],[sns.color_palette('dark')[0],sns.color_palette('dark')[1]]):
        str_peak    = ''
        for ikl in peak:
            str_peak += str(ikl)
        fit_pattern=f'../../.data_SACLA/fits/C{str_peak}_peak_rrun_modelresult.sav'
        for delay_id, delay in enumerate([8,9,10,11,12,13,14,15,16,17,18,19]):
            fig_spectrum, ax_spectrum       = plt.subplots(figsize=(7,7))

            run                         = SACLA_dict[delay]

            fitfile                     = re.sub('run',str(run),fit_pattern)
            fit                         = load_modelresult(fitfile)

            run_data, ref_data          =   load_SACLA(run)

            thetamin, thetamax          = thetarange[0], thetarange[1]
            mask                        = np.logical_and(run_data[:,0] > thetamin, run_data[:,0] < thetamax)

            try:
                T, Tmin, Tmax, V, Vmin, Vmax               = T_from_peak(fit.values['C_center'], fit.params['C_center'].stderr,peak)
            except:
                T, Tmin, Tmax, V, Vmin, Vmax                = T_from_peak(fit.values['C_center'], fit.values['C_center']*0.01 ,peak)
            if np.isnan(Tmax):
                Tmax                    = 2*T - Tmin
            natoms_uc, mass_atomic_C    =   8, 12.011
            mass_uc                     = natoms_uc * mass_atomic_C
            proportion                  = 1.66053907
            rho, rho_min, rho_max       = mass_uc/V*proportion, mass_uc/Vmin*proportion, mass_uc/Vmax*proportion

            ax.errorbar(delay, T, yerr=np.array([[T-Tmin,Tmax-T]]).T, fmt="o",c=color,markersize=5)
            print(rho_min)
            if rho_min>1.:
                ax_vol.errorbar(delay, rho, yerr=np.array([[rho-rho_min],[rho_max-rho]]), fmt="o",c=color,markersize=5)
            ax_spectrum.errorbar(run_data[::5,0], run_data[::5,1], yerr=run_data[::5,2], fmt="o",c=color,markersize=0.1)
            ax_spectrum.set_xlabel(r"Scattering angle $2\theta$ [$\circ$]",fontsize=30)
            ax_spectrum.set_ylabel(r"Intensity [au]",fontsize=30)
            fig_spectrum.tight_layout()
            title = f'../../../../../W_PhD_Articles/Heuser_ND_recovery/figures/lineout_{str(run)}.pdf'
            plt.savefig(title, format='pdf', dpi=1200)
            plt.close(fig_spectrum)

    ax.set_xlabel(r"Probe laser delay [ns]",fontsize=30)
    ax.set_ylabel(r"Diamond temperature [K]",fontsize=30)
    ax_vol.set_xlabel(r"Probe laser delay [ns]",fontsize=30)
    ax_vol.set_ylabel(r"Diamond density [g/cc]",fontsize=30)
    #ax.axhline(y=weighted_mean,xmin=9.,xmax=19.,c="r",linewidth=10)
    #ax.plot([9,19],[2700,2700],c='r',linestyle='dashed',label='After breakout mean T')
    ax.tick_params(axis='both', which= 'major', labelsize=25)
    ax_vol.tick_params(axis='both', which= 'major', labelsize=25)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    legend_elements = [ax.errorbar([], [], yerr=[0.01], fmt="o",c=sns.color_palette('dark')[0],markersize=5,label='(111) peak')]
    legend_elements = [ax.errorbar([], [], yerr=[0.01], fmt="o",c=sns.color_palette('dark')[1],markersize=5,label='(220) peak')]

    #legend_elements_vol = [ax_vol.errorbar([], [], yerr=[0.01], fmt="o",c=sns.color_palette('dark')[0],markersize=5,label='(111) peak')]
    #legend_elements_vol = [ax_vol.errorbar([], [], yerr=[0.01], fmt="o",c=sns.color_palette('dark')[1],markersize=5,label='(220) peak')]
    
    ax.legend(fontsize=20)
    #ax_vol.legend(fontsize=20)
    sns.set(style="ticks")
    sns.set_style("whitegrid")
    fig.tight_layout()
    plt.savefig('../../../../../W_PhD_Articles/Heuser_ND_recovery/figures/TempDiamond.pdf', format='pdf', dpi=1200)

if __name__=='__main__':
    run_filter          =   [412,414,416,418,420,422,424,426,430,432,434,436,617,619]
    SACLA_shots         =   pd.read_pickle("../../.data_SACLA/logbook/PET.pkl")
    PET75_shots         =   SACLA_shots[SACLA_shots['target'] == 'PET75']
    PET75_shots         =   PET75_shots[PET75_shots['run'].isin(run_filter)]
    runs_SACLA              = [424,436,422,434,420,432,418,430,412,426,414,617]
    delays_SACLA            = np.array([PET75_shots.loc[PET75_shots['run'] == run]['delay'].values[0] for run in runs_SACLA])
    SACLA_dict              = {}
    for delay, run in zip(delays_SACLA,runs_SACLA):
        SACLA_dict.update({delay: run})

    plot_waterfall()

    import sys
    import re
    sys.path.append('../xrd')
    from expansion_coefficient import T_from_peak
    plot_temperature_estimates()