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
import sys
######################Loading data#################
###################################################################

data = []
pwd = '../Distance_for_MSM'

for i in range(1,10):
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

for ic in range(2,8):
    
    tica = pyemma.coordinates.tica(data,lag=100,dim=ic)
    print(tica.dimension())
    tica_output = tica.get_output()
    tica_concatenated = np.concatenate(tica_output)
    
    ######################Clustering Dataset#################
    ###################################################################
    
    n_clustercenters = [100,500,1000,2000]
    
    scores = np.zeros((len(n_clustercenters), 5))
    for n, k in enumerate(n_clustercenters):
        print (n)
        print(k)
        for m in range(5):
            with pyemma.util.contexts.settings(show_progress_bars=False):
                _cl = pyemma.coordinates.cluster_kmeans(
                    tica_output, k=k, max_iter=50, stride=50)
                _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 100)
                scores[n, m] = _msm.score_cv(
                    _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))
        np.savetxt(f'VAMP_score_data/vamp_scores_{ic}.dat',scores) 


