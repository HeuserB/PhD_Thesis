import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import logging
import sys
import pandas as pd
import seaborn as sns



def load_hf5(h5file):
    f   = h5.File(h5file,'r')
    equilibration_thermo    = pd.DataFrame()
    expansion_thermo    = pd.DataFrame()
    data_frames    = [equilibration_thermo, expansion_thermo]
    groupnames  = ['ThermoData/equilibration','ThermoData/expansion']

    for purpose in [0,1]:
        for key in f[groupnames[purpose]].keys():
            tmp =   str(groupnames[purpose] + '/' + key)
            data_frames[purpose][key]    = f[tmp][:]
    f.close()
    return data_frames

def running_average(x, avg_size):
    # calculating the average of an array over avg_size chunks
    avg     =   np.average(np.array(x[:-(len(x)%avg_size)]).reshape(-1,avg_size), axis=1)
    avg     =   np.hstack([avg,np.array(np.sum(x[-(len(x)%avg_size):])/(len(x)%avg_size)).reshape(-1)])
    return avg


def plot_thermo(filename_base):
    label   = {0:50,1:100,2:200}
    files   = [filename_base + str(label[i]) + 'ps.h5' for i in range(3)]  
    logger.info(f'Loading the files :{files}')
    data_frames = [load_hf5(file) for file in files]

    fig, axes   =   plt.subplots(2,2,figsize=(20,10))
    for i in range(3):
        average_step         = running_average(data_frames[i][0]['Step'], 50)
        average_T0           = running_average(data_frames[i][0]['Temp'], 50)
        average_vol         = running_average(data_frames[i][1]['Volume'], 50)
        average_T           = running_average(data_frames[i][1]['Temp'], 50)
        axes[0,1].plot(average_vol, average_T,label=f"{label[i]}ps")
        axes[0,0].plot(average_step, average_T0,label=f"{label[i]}ps",c='k')

        average_P           = running_average(data_frames[i][1]['Press'], 50) / 10000
        average_P0           = running_average(data_frames[i][0]['Press'], 50) / 10000
        axes[1,1].plot(average_vol, average_P,label=f"{label[i]}ps")
        axes[1,0].plot(average_step, average_P0,label=f"{label[i]}ps",c='k')

    axes[0,1].get_shared_x_axes().join(axes[0,1], axes[1,1])
    axes[0,1].get_shared_x_axes().join(axes[0,0], axes[1,0])
    axes[0,0].get_shared_y_axes().join(axes[0,0], axes[0,1])
    axes[1,0].get_shared_y_axes().join(axes[1,0], axes[1,1])
    axes[0,1].legend(fontsize=20)
    axes[1,1].legend(fontsize=20)
    axes[0,0].spines['bottom'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[0,0].spines['right'].set_visible(False)
    axes[0,1].spines['bottom'].set_visible(False)
    axes[0,1].spines['top'].set_visible(False)
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].spines['left'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['left'].set_visible(False)
    #axes[0,0].set_xticks([])
    #axes[0,1].set_xticks([])
    #axes[0,1].set_yticks([])
    #axes[1,1].set_yticks([])
    #axes[0,0].get_xaxis().set_visible(False)
    axes[0,0].tick_params(axis='both', which= 'major', labelsize=20)
    axes[0,1].tick_params(axis='both', which= 'major', labelsize=20)
    axes[1,0].tick_params(axis='both', which= 'major', labelsize=20)
    axes[1,1].tick_params(axis='both', which= 'major', labelsize=20)

    axes[0,0].set_ylabel(r'Temperature [K]',fontsize=20)
    axes[1,0].set_ylabel(r'Pressure [GPa]',fontsize=20)
    axes[1,1].set_xlabel(r'Volume [$\AA^3$]',fontsize=20)
    axes[1,0].set_xlabel(r'Steps',fontsize=20)

    axes[0,0].grid(b=True,linestyle='--')
    axes[0,1].grid(b=True,linestyle='--')
    axes[1,0].grid(b=True,linestyle='--')
    axes[1,1].grid(b=True,linestyle='--')
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Creating an object
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO)
    logger.setLevel(logging.INFO)

    if len(sys.argv) == 1:
        logger.error('Please provide filenames!')

    data_frames = plot_thermo(sys.argv[1])