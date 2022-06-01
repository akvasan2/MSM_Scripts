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
    if (i!=6 and i!=7):
        line_PROA = np.loadtxt(f'{pwd}/Run_{i}-PROA.dat')
        line_PROA = np.hstack((line_PROA,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROA.dat')))

    line_PROB = np.loadtxt(f'{pwd}/Run_{i}-PROB.dat')
    line_PROB = np.hstack((line_PROB,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROB.dat')))

    if (i!=0 and i!=8):
        line_PROC = np.loadtxt(f'{pwd}/Run_{i}-PROC.dat')
        line_PROC = np.hstack((line_PROC,np.loadtxt(f'{pwd_d120add}/Run_{i}-PROC.dat')))
    
    if (i!=6 and i!=7):
        data.append(line_PROA)

    data.append(line_PROB)

    if (i!=0 and i!=8):
        data.append(line_PROC)

#    print(line_PROA.shape)

data_concatenated = np.concatenate(data)

tica = pyemma.coordinates.tica(data,lag=tica_lag,dim=ic)

tica_output = tica.get_output()

tica_concatenated = np.concatenate(tica_output)


with open('data.pickle', 'wb') as f:
    pickle.dump(data,f)

with open('data_concat.pickle', 'wb') as f:
    pickle.dump(data_concatenated,f)

with open('tica.pickle', 'wb') as f:
    pickle.dump(tica,f)

with open('tica_output.pickle', 'wb') as f:
    pickle.dump(tica_output,f)

with open('tica_concat.pickle', 'wb') as f:
    pickle.dump(tica_concatenated,f)


#######################Clustering Dataset#################
####################################################################

cluster = pyemma.coordinates.cluster_kmeans(tica, k=num_clust, max_iter=50, stride=10)

dtrajs_concatenated = np.concatenate(cluster.dtrajs) ## dtrajs is the cluster at each frame of each traject.  dtrajs_concatenated concatenates each 

cluster.save('distance.pyemma', model_name='distance_ompf_cluster', overwrite=True) ## distance.pyemma  can contain clustering, msm, and bayesian msm information
cluster.save_dtrajs(prefix="cluster",output_dir="dtrajs_files",extension='.dtraj') ## need to save cluster dtrajs independently

######################### If you already have clusters stored, you can just load them here by uncommenting these lines####################

#cluster = pyemma.load('distance.pyemma', model_name='distance_ompf_cluster')
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

######################Creating MSM #################
###################################################################

msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=msm_lag, dt_traj='10 ps')
bayesian_msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=msm_lag, dt_traj='10 ps', conf=0.95)

###################################################################
############################ Store MSM #######################################
###################################################################

msm.save('distance.pyemma', model_name='distance_ompf_msm', overwrite=True)
bayesian_msm.save('distance.pyemma', model_name='distance_ompf_bayesian_msm', overwrite=True)

print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))

############ IF you already have msm you can load them here ################

msm = pyemma.load('distance.pyemma', model_name='distance_ompf_msm')
bayesian_msm = pyemma.load('distance.pyemma', model_name='distance_ompf_bayesian_msm')

###########################If you want to plot stationary distribution/eigenvectors uncomment these lines#######################################
##################################################################

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
pyemma.plots.plot_contour(
        *tica_concatenated[:, :2].T,
    msm.pi[dtrajs_concatenated],
    ax=axes[0],
    mask=True,
    cbar_label='stationary distribution')
pyemma.plots.plot_free_energy(
        *tica_concatenated[:, :2].T,
    weights=np.concatenate(msm.trajectory_weights()),
    vmax = 6,
    kT=0.6,
    ax=axes[1],
    legacy=False)
for ax in axes.flat:
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
axes[0].set_title('Stationary distribution', fontweight='bold')
axes[1].set_title('Reweighted free energy surface', fontweight='bold')
fig.tight_layout()
plt.ylim(-4,4)
plt.xlim(-4,4)
plt.savefig('stat_reweight_freeen.png')
plt.close()

eigvec = msm.eigenvectors_right()


print('first eigenvector is one: {} (min={}, max={})'.format(
    np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
            *tica_concatenated[:, :2].T, eigvec[dtrajs_concatenated, i + 1], ax=ax, cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(i + 2), mask=True)
    ax.set_xlabel('IC1')
    ax.set_aspect('equal')
plt.xlim(-4,4)
plt.ylim(-4,4)
axes[0].set_ylabel('IC2')
fig.tight_layout()
plt.savefig('Eigvecs_tica.png')

