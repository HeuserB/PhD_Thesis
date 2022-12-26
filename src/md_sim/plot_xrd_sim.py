import numpy as np
import matplotlib.pyplot as plt
import sys
import logging 
import re

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)

def load_xrd_file(filename):
    try:
        with open(filename, 'r') as file:
            content = file.readlines()
            print(content[1])
            print(re.split('\s',content[3]))
            TimeStep = int(re.split('\s',content[3])[0])
            nBins = int(re.split('\s',content[3])[1])
            counts = float(re.split('\s',content[3])[2])
            missing_counts = float(re.split('\s',content[3])[3])
            min = float(re.split('\s',content[3])[4])
            min = float(re.split('\s',content[3])[5])
            #startline = 4
            timesteps_total = int((len(content) - 3 ) / (1 + nBins))
            hist = np.empty([timesteps_total,4,nBins])
            for timestep in range(timesteps_total):
                startline = 4 + timestep * (1 + nBins)
                for lineid, linenumber in enumerate(range(startline,startline+nBins)):
                    for id in range(4):
                        hist[timestep,id,lineid] = float(re.split('\s',content[linenumber])[id])
        logger.info(f'Loaded a total of {timesteps_total + 1} timesteps from file {filename}.')
        return hist 
    except:
        logger.error('Error loading file {filename}!')
        return None

if __name__=='__main__':
    hist = load_xrd_file(sys.argv[1])
    specturm = np.loadtxt(sys.argv[2])
    fig, ax = plt.subplots()
    for i in range(len(hist)):
        ax.plot(hist[i,1,:],hist[i,2,:]/np.max(hist[i,2,:]), label=f'Snapshot {i}')
    ax.plot(specturm[:,0], specturm[:,1]/np.max(specturm[:,1]))
    ax.set_xlabel(r'$2\theta$')
    ax.set_ylabel(r'Simulated intensity')
    ax.legend()
    plt.show()