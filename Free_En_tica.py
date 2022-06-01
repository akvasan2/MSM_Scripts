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
import matplotlib
import pickle

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"

######################Loading data#################
###################################################################

with open('data.pickle','rb') as f:
    data = pickle.load(f)

with open('data_concat.pickle','rb') as f:
    data_concatenated = pickle.load(f)

with open('tica.pickle','rb') as f:
    tica = pickle.load(f)

with open('tica_output.pickle','rb') as f:
    tica_output = pickle.load(f)

with open('tica_concat.pickle','rb') as f:
    tica_concatenated = pickle.load(f)

######################Converting data to tica #################
###################################################################

#tica = pyemma.coordinates.tica(data,lag=100,dim=5)
#tica_output = tica.get_output()
#tica_concatenated = np.concatenate(tica_output)

######################Load Clusters, MSMs#################
###################################################################
#cluster_dtrajs = []
#pwd = 'dtrajs_files'
#for i in range(0,30):
#    # i = 1
#    line_cluster = np.loadtxt(f'{pwd}/cluster_{i}.dtraj')
#    
#    cluster_dtrajs.append(line_cluster)
#
#dtrajs_concatenated_float = np.concatenate(cluster_dtrajs)
#dtrajs_concatenated = [int(dc) for dc in dtrajs_concatenated_float]

cluster = pyemma.load('distance.pyemma', model_name='distance_ompf_cluster')

msm = pyemma.load('distance.pyemma', model_name='distance_ompf_msm')

bayesian_msm = pyemma.load('distance.pyemma', model_name='distance_ompf_bayesian_msm')

#############################Plot stationary distribution/eigenvectors#######################################
####################################################################

pyemma.plots.plot_free_energy(
        *tica_concatenated[:, :2].T,
    weights=np.concatenate(msm.trajectory_weights()),
    vmax = 6,
    kT=0.62,
    ax=None,
    legacy=False)
#plt.xticks(np.arange(-1.5,3,step=0.5))
#plt.yticks(np.arange(-1.5,5,step=0.5))

#plt.grid()
plt.savefig('freeen_tica.grid.G119D.png',bbox_inches='tight',dpi=300)
plt.close()

