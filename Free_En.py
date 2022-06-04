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
import matplotlib
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


######################Load Clusters, MSMs#################
###################################################################
cluster_dtrajs = []
pwd = 'dtrajs_files'
for i in range(0,18):
    line_cluster = np.loadtxt(f'{pwd}/cluster_{i}.dtraj')
    
    cluster_dtrajs.append(line_cluster)

dtrajs_concatenated_float = np.concatenate(cluster_dtrajs)
dtrajs_concatenated = [int(dc) for dc in dtrajs_concatenated_float]

cluster = pyemma.load('distance.pyemma', model_name='distance_ompf_cluster')

msm = pyemma.load('distance.pyemma', model_name='distance_ompf_msm')
bayesian_msm = pyemma.load('distance.pyemma', model_name='distance_ompf_bayesian_msm')

eigvec = msm.eigenvectors_right()
#############################Plot stationary distribution/eigenvectors#######################################
####################################################################
pos_corr_feat = 9
pos_label = 'E117-Y22'
neg_corr_feat = 26#6
neg_label ='D119-R132' #'D121-R82'
data_proj = np.vstack((data_concatenated[:,pos_corr_feat],data_concatenated[:,neg_corr_feat])).T
pyemma.plots.plot_free_energy(                        
    *data_proj.T,
    weights=np.concatenate(msm.trajectory_weights()),
    kT=0.6,
    vmin=0,
    ax=None,)
    #legacy=False)
plt.xlim(0,20)
plt.ylim(0,20)
plt.yticks(np.arange(0, 21, step=5))
plt.xticks(np.arange(0, 21, step=5))
plt.xlabel(pos_label)
plt.ylabel(neg_label)
plt.savefig(f'Images/Free_en{pos_label}_{neg_label}.G120D.new.png',bbox_inches='tight',dpi=300)        
plt.close()

pyemma.plots.plot_free_energy(                        
    *data_proj.T,
#    weights=np.concatenate(msm.trajectory_weights()),
    kT=0.6,
    vmin=0,
    ax=None,)
    #legacy=False)
plt.xlim(0,20)
plt.ylim(0,20)
plt.yticks(np.arange(0, 21, step=5))
plt.xticks(np.arange(0, 21, step=5))
plt.xlabel(pos_label)
plt.ylabel(neg_label)
plt.savefig(f'Images/Free_en{pos_label}_{neg_label}_unweighted.G120D.new.png',bbox_inches='tight',dpi=300)        
plt.close()



































































