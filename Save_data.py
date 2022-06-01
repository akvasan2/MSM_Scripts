import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import mdshare
import pyemma
from pyemma.util.contexts import settings
import matplotlib.pyplot as plt
import numpy as np
import mdshare
import pyemma
import mdtraj
import pickle
##################### Define parameters ##################
ic = 5 ###### num IC for TICA
tica_lag =100 ####### lag time for tica= 100 steps or 1 ns
msm_lag = 200 ####### lag time for msm = 200 steps or 2 ns
num_clust = 1000 ##### num microstates

######################Loading data#################
###################################################################
data = []
pwd = '../Distance_for_MSM'
pwd_d120add = '../Distance_for_MSM_additional_D120'

for i in range(0,10):
    # i = 1
    if (i!=1 and i!=4 and i!=5 and i!=6 and i!=7):
        line_PROA = np.loadtxt(f'{pwd}/Run_{i}-PROA.dat')
        line_PROA = np.hstack((line_PROA,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROA.dat')))

    if (i!=1 and i!=3):
        line_PROB = np.loadtxt(f'{pwd}/Run_{i}-PROB.dat')
        line_PROB = np.hstack((line_PROB,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROB.dat')))

    if (i!=0 and i!=2 and i!=4 and i!=5 and i!=6):
        line_PROC = np.loadtxt(f'{pwd}/Run_{i}-PROC.dat')
        line_PROC = np.hstack((line_PROC,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROC.dat')))
    
    if (i!=1 and i!=4 and i!=5 and i!=6 and i!=7):
        data.append(line_PROA)

    if (i!=1 and i!=3):
        data.append(line_PROB)

    if (i!=0 and i!=2 and i!=4 and i!=5 and i!=6):
        data.append(line_PROC)

    print(line_PROA.shape)

data_concatenated = np.concatenate(data)

with open('data.pickle', 'wb') as f:
    pickle.dump(data,f)

with open('data_concat.pickle', 'wb') as f:
    pickle.dump(data_concatenated,f)


#
#######################Converting data to tica #################
####################################################################
#
tica = pyemma.coordinates.tica(data,lag=tica_lag,dim=ic)

tica_output = tica.get_output()

tica_concatenated = np.concatenate(tica_output)

with open('tica.pickle', 'wb') as f:
    pickle.dump(tica,f)

with open('tica_output.pickle', 'wb') as f:
    pickle.dump(tica_output,f)

with open('tica_concat.pickle', 'wb') as f:
    pickle.dump(tica_concatenated,f)


