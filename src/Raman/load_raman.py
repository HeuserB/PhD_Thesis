import numpy as np
import re
import os

pattern_spectrum            = 'Spektrum[0-9]*.*'
pattern_X                   = 'X[-]*[0-9]+,[0-9]+'
pattern_Y                   = 'Y[-]*[0-9]+,[0-9]+'

single_spectrum_pattern     =  "[0-9]*[\s]*[_]*[a-zA-Z]*[\s]*[_]*[0-9]*[\s]*[_]*[0-9]*_[0-9]*_[0-9]*.csv"
run_pattern                 =  "[0-9]*_[0-9]*_[0-9]*.csv"

def load_single_raman(run_id, root_dir):
    for root, dirs, files in os.walk(root_dir):
        for filename in files:
            if re.match(single_spectrum_pattern,filename) != None:
                ID              =   re.search(run_pattern,filename)[0][:-4]
                if ID == run_id:
                    file            =   os.path.join(root,filename)
                    shift, I        =   np.loadtxt(file,skiprows=1).T
                    return shift, I
    print(f"Run not found!")
    return None, None


def load_raman(file,scan=False):
    with open(file) as f:
        if scan == True:
            txt = "".join(f.readlines())
        else:
            txt = "".join(f.readlines()[1:])
    txt = re.sub(',','.',txt.strip())
    if scan == True:
        #return txt
        return(np.fromstring(txt.strip(),sep='\n').reshape(-1))
    return(np.fromstring(txt,sep=' ').reshape(-1,2))