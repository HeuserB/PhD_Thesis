import pandas as pd
import numpy as np

SACLA_target_dict   = {0:'PET75',1:'AL/PET100',2:'PET1000',3:'PET38/Al6.5',4:'Al10/PET75',5:'PET100/Al'}
SACLA_data          = \
    np.array([[23,25,27,29,31,33,35,76,78,80,82,146,148,150,152,155,157,159,163,165,167,169,171,173,182,185,188,192,195,199,204,211,213,215,217,222,224,226,307,309,311,313,316,320,323,325,327,329,331,333,335,358,360,362,365,370,372,375,412,414,416,418,420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,611,613,615,617,619,621,623,625,627,629,631,631],\
    [1162359,1162361,1162363,1162365,1162367,1162369,1162371,1162415,1162417,1162419,1162421,1162485,1162487,1162489,1162491,1162494,1162496,1162498,1162502,1162504,1162506,1162508,1162510,1162512,1162521,1162524,1162527,1162531,1162534,1162537,1162540,1162547,1162549,1162551,1162553,1162558,1162560,1162562,1162645,1162647,1162649,1162651,1162654,1162658,1162661,1162663,1162665,1162667,1162669,1162671,1162673,1162697,1162699,1162701,1162704,1162709,1162711,1162714,1162751,1162753,1162755,1162757,1162759,1162761,1162763,1162765,1162767,1162769,1162771,1162773,1162775,1162777,1162779,1162781,1162783,1162785,1162787,1162789,1162791,1162948,1162950,1162952,1162954,1162956,1162958,1162960,1162962,1162964,1162966,1162968,1162970],\
    [1162358,1162360,1162362,1162364,1162366,1162368,1162370,1162414,1162416,1162418,1162420,1162484,1162486,1162488,1162490,1162493,1162495,1162497,1162501,1162503,1162505,1162507,1162509,1162511,1162520,1162523,1162526,1162530,1162533,1162536,1162539,1162546,1162548,1162550,1162552,1162557,1162559,1162561,1162644,1162646,1162648,1162650,1162653,1162657,1162660,1162662,1162664,1162666,1162668,1162670,1162672,1162696,1162698,1162700,1162703,1162708,1162710,1162713,1162750,1162752,1162754,1162756,1162758,1162760,1162762,1162764,1162766,1162768,1162770,1162772,1162774,1162776,1162778,1162780,1162782,1162784,1162786,1162788,1162790,1162947,1162949,1162951,1162953,1162955,1162957,1162959,1162961,1162963,1162965,1162967,1162969],\
    [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,1162522,1162525,1162528,None,1162535,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None],\
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2,0.2,0.2,0.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
    [2,5,4,3,3,4,5,5,5,5,6,6,6,10,8,8,12,7,4,5,4,3,6,10,3,3,3,4,7,5,5,3,4,5,2,2,2,4,5,8,10,11,13,15,12,13,11,10,14,15,13,12,14,14,12,16,16,12,11,9,7,13,15,17,19,10,8,12,14,16,18,11,12,15,18,12,9,6,4,12,14,10,8,6,4,3,8,8,8,8,8],\
    [1.24E+01,9.68E+00,1.20E+01,1.32E+01,8.34E+00,9.30E+00,1.07E+01,3.41E+00,1.95E+00,4.45E+00,7.69E+00,2.91E+00,5.90E+00,4.00E+00,3.98E+00,3.96E+00,3.40E+00,2.79E+00,2.83E+00,2.39E+00,1.95E+00,2.41E+00,2.31E+00,2.32E+00,1.55E+01,1.32E+01,1.74E+01,1.60E+01,1.62E+01,1.74E+01,2.01E+01,2.45E+00,2.01E+00,2.03E+00,2.51E+00,1.87E+00,3.03E+00,2.88E+00,6.80E+00,6.05E+00,7.16E+00,7.54E+00,6.59E+00,4.88E+00,5.60E+00,6.23E+00,6.21E+00,5.26E+00,7.05E+00,7.06E+00,6.32E+00,5.21E+00,7.46E+00,6.24E+00,6.16E+00,6.16E+00,4.33E+00,5.10E+00,6.13E+00,8.12E+00,5.38E+00,6.08E+00,5.12E+00,6.23E+00,5.90E+00,5.54E+00,5.73E+00,5.07E+00,4.41E+00,4.06E+00,4.07E+00,3.36E+00,4.74E+00,6.35E+00,7.02E+00,4.84E+00,6.35E+00,7.36E+00,7.42E+00,6.10E+00,6.99E+00,7.51E+00,7.38E+00,6.46E+00,7.01E+00,6.93E+00,4.98E+00,5.28E+00,2.58E+00,3.29E+00,4.05E+00],\
    [0.8,0.8,0.8,0.8,0.6,0.6,0.6,0.2,0.15,0.25,0.4,0.2,0.3,0.3,0.25,0.2,0.2,0.15,0.15,0.13,0.13,0.13,0.13,0.13,0.13,0.2,0.5,1,0.2,1,0.13,0.13,0.13,0.13,0.13,0.13,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.42,0.4,0.4,0.4,0.4,0.4,0.3,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.3,0.3,0.3,0.3,0.3,0.25,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.35,0.25,0.25,0.28],\
    [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,0,0,0,0,0,0,4,4,4,4,4,4,4,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0]]
    )
#for i in SACLA_data:
#   print(len(i))    

SACLA_shots         = pd.DataFrame(SACLA_data.T,columns=['run','shot_id','ref_id','post_id','Xray_attenuation',\
   'delay','E_on_sample','laser_transmission','target'])

datatypes           =  [int, int, int, int, float, int, float, float, int]


[SACLA_shots.__setitem__(key,SACLA_shots[key].convert_dtypes(datatypes[id])) for id, key in enumerate(SACLA_shots.keys())];
SACLA_shots['target'] = [SACLA_target_dict[target_id] for target_id in SACLA_shots['target'].values]
#for id, key in enumerate(SACLA_shots.keys()): 
#    SACLA_shots[key] = SACLA_shots[key].convert_dtypes(datatypes[id])
SACLA_shots.to_pickle("../../.data_SACLA/logbook/PET.pkl")