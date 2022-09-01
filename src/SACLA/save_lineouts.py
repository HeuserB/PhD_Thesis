import pyFAI
import pandas as pd
import sys
import numpy as np
import imageio

sys.path.append('../diamond_analysis')
from load_data import write_SACLA_lineout


mask_SACLA  = "/home/benjamin//Nextcloud/Work/W_PhD/W_PhD_Git/PhD_Git/PhD_Thesis/.data_SACLA/mask/SACLA_extensive.mask"
SACLA_mask  = np.array(imageio.imread(mask_SACLA),dtype=int);

poni_file_SACLA = "/home/benjamin//Nextcloud/Work/W_PhD/W_PhD_Git/PhD_Git/PhD_Thesis/.data_SACLA/poni/Q_SACLA_Au.poni"
SACLA_ai = pyFAI.load(poni_file_SACLA);

data_dir_SACLA      =   "/home/benjamin/Nextcloud/Work/W_PhD/W_PhD_Analysis/SACLA_2022/stevenson_2022a/fpd"
SACLA_shots         =   pd.read_pickle("/home/benjamin//Nextcloud/Work/W_PhD/W_PhD_Git/PhD_Git/PhD_Thesis/.data_SACLA/logbook/PET.pkl")

for run_id in range(len(SACLA_shots)):
    run         =   SACLA_shots.iloc[run_id]['run']
    shot_id     =   SACLA_shots.iloc[run_id]['shot_id']
    ref_id      =   SACLA_shots.iloc[run_id]['ref_id']
    write_SACLA_lineout(run, shot_id, data_dir_SACLA, SACLA_ai, SACLA_mask,ref_id)