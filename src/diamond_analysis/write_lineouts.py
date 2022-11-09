import pandas as pd
import numpy as np

from load_data import write_SACLA_lineout

if __name__=='main':
    run_filter          =   [412,414,416,418,420,422,424,426,430,432,434,436,617,619]
    print(f"Writing SACLA shots {run_filter} ")

    SACLA_shots         =   pd.read_pickle("../../.data_SACLA/logbook/PET.pkl")
    PET75_shots         =   SACLA_shots[SACLA_shots['target'] == 'PET75']
    PET75_shots         =   PET75_shots[PET75_shots['run'].isin(run_filter)]

    for shot in run_filter:
        write_SACLA_lineout(shot, PET75_shots[PET75_shots['run'] == shot]['shot_id'].values[0] , "/media/benjamin/BenDrive/Data/SACLA_2022/stevenson_2022a/fpd/", SACLA_ai, SACLA_mask, PET75_shots[PET75_shots['run'] == shot]['ref_id'].values[0])