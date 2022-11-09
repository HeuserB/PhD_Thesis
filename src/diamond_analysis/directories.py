# define the data directory
dir_data = "../../.data_LW03"
import re

# the logbook file
logbook = "../../.data_LW03/logbook_LW03_subset.ods"

# the calibration .poni files created by Dioptas
poni_file_q0    = "../../.data_LW03/poni/Quad0_CeO2.poni"
poni_file_q1    = "../../.data_LW03/poni/Quad1_CeO2.poni"
poni_file_q2    = "../../.data_LW03/poni/Quad2_CeO2.poni"
poni_file_q3    = "../../.data_LW03/poni/Quad3_CeO2.poni"

poni_file_q0    = "../../.data_LW03/poni/Quad0_MS.poni"
poni_file_q1    = "../../.data_LW03/poni/Quad1_MS.poni"
poni_file_q2    = "../../.data_LW03/poni/Quad2_MS.poni"
poni_file_q3    = "../../.data_LW03/poni/Quad3_MS.poni"
poni_file_SACLA = "../../.data_SACLA/poni/Q_SACLA.poni"

# the mask files for the quad detectors
mask_q0     = "../../.data_LW03/mask/Quad0_ext.mask"
mask_q1     = "../../.data_LW03/mask/Quad1_ext.mask"
mask_q2     = "../../.data_LW03/mask/Quad2.mask"
mask_q3     = "../../.data_LW03/mask/Quad3.mask"
mask_SACLA  = "../../.data_SACLA/mask/SACLA_extensive.mask"

# define the pattern for a background file
pattern_runs = re.compile("r[0-9]*_bkgCorrected")


inst_file = '../../.data_LW03/instprm/LW03.instprm'
cif_file = '../../.data_LW03/cif/diamond.cif'