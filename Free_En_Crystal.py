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
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"

######################Loading data#################
###################################################################

data = []
pwd = '/Scr/nandan/OmpF/Ompf-Trimer/equilibration_150mM_NaCl/Analysis/Distance_for_MSM'

for i in range(1,11):
    # i = 1
    line_PROA = np.loadtxt(f'{pwd}/Run_{i}-PROA.dat')
    line_PROB = np.loadtxt(f'{pwd}/Run_{i}-PROB.dat')
    line_PROC = np.loadtxt(f'{pwd}/Run_{i}-PROC.dat')
    
    data.append(line_PROA)
    data.append(line_PROB)
    data.append(line_PROC)

data_concatenated = np.concatenate(data)

######################Converting data to tica #################
###################################################################

tica = pyemma.coordinates.tica(data,lag=100,dim=5)
tica_output = tica.get_output()
tica_concatenated = np.concatenate(tica_output)

######################Load Clusters, MSMs#################
###################################################################
cluster_dtrajs = []
pwd = 'dtrajs_files'
for i in range(0,30):
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
neg_corr_feat = [1]#,15]
y = 1
neg_label = ['D121-R132','D121-R82']#,'E117-Y302']
it = 0

it+=1
it = 0
y_int = int(y)
data_proj = np.vstack((data_concatenated[:,pos_corr_feat],data_concatenated[:,y_int])).T
pyemma.plots.plot_free_energy(                        
    *data_proj.T,
    weights=np.concatenate(msm.trajectory_weights()),
    kT=0.6,
    vmin=0,
    #vmax=6,
    ax=None,
    legacy=False)

data_cryst = np.loadtxt('Crystal_Structure_Data/Distance.dat')
plt.scatter(data_cryst[:,1],data_cryst[:,0],c='black',s=20)

plt.xlim(0,20)
plt.ylim(0,20)
plt.yticks(np.arange(0, 21, step=5))
plt.xticks(np.arange(0, 21, step=5))
plt.xlabel(pos_label)
plt.ylabel(neg_label[it])
plt.savefig(f'Images/FreeEn_Features/Free_en{pos_label}_{neg_label[it]}_cryst.png',bbox_inches='tight',dpi=300)        
plt.close()
