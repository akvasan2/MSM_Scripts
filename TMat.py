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

##################### Define parameters ##################
ic = 5 ###### num IC for TICA
tica_lag =100 ####### lag time for tica= 100 steps or 1 ns
msm_lag = 200 ####### lag time for msm = 200 steps or 2 ns
num_clust = 1000  ##### num microstates

######################Loading data#################
###################################################################
data = []
pwd = '/Scr/nandan/OmpF/Ompf-Trimer/equilibration_150mM_NaCl/Analysis/Distance_for_MSM'

for i in range(1,11): # 1,11
    # i = 1
    line_PROA = np.loadtxt(f'{pwd}/Run_{i}-PROA.dat')
    line_PROB = np.loadtxt(f'{pwd}/Run_{i}-PROB.dat')
    line_PROC = np.loadtxt(f'{pwd}/Run_{i}-PROC.dat')
    
    data.append(line_PROA)
    data.append(line_PROB)
    data.append(line_PROC)

data_concatenated = np.concatenate(data)

print(data_concatenated.shape)

######################Converting data to tica #################
###################################################################

tica = pyemma.coordinates.tica(data,lag=tica_lag,dim=ic)

tica_output = tica.get_output()

tica_concatenated = np.concatenate(tica_output)

######################Clustering Dataset#################
###################################################################

cluster_dtrajs = []
pwd = 'dtrajs_files'
for i in range(0,30):
    # i = 1
    line_cluster = np.loadtxt(f'{pwd}/cluster_{i}.dtraj')
    
    cluster_dtrajs.append(line_cluster)

dtrajs_concatenated_float = np.concatenate(cluster_dtrajs)
dtrajs_concatenated = [int(dc) for dc in dtrajs_concatenated_float]

cluster = pyemma.load('distance.pyemma', model_name='distance_ompf_cluster')

#dtrajs_concatenated = np.concatenate(cluster.dtrajs)

msm = pyemma.load('distance.pyemma', model_name='distance_ompf_msm')
bayesian_msm = pyemma.load('distance.pyemma', model_name='distance_ompf_bayesian_msm')

######################Creating MSM and running CK Test#################
###################################################################

np.savetxt('Count_matrix.dat',msm.count_matrix_active)
np.savetxt('Transition_matrix.dat',msm.transition_matrix)



print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))

############################Plot stationary distribution/eigenvectors#######################################
###################################################################

