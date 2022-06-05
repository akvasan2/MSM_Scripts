import matplotlib.pyplot as plt
import numpy as np
import mdshare
import pyemma
import mdtraj


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

#######################Converting data to tica #################
####################################################################

tica = pyemma.coordinates.tica(data,lag=100,dim=5)

tica_output = tica.get_output()

tica_concatenated = np.concatenate(tica_output)

######################Calculating ITS#################
###################################################################

for i, k in enumerate([1000]):

    cluster = pyemma.coordinates.cluster_kmeans(tica, k=k, max_iter=50, stride=10)
    ITS = pyemma.msm.its(cluster.dtrajs, lags=[1,2,5,10,20,50, 100, 200,300,400,500,600,700,800,900,1000], nits=10, errors = 'bayes')

np.savetxt(f'ITS_bayesian.tscale.dat',ITS.timescales)
np.savetxt(f'ITS_bayesian.mean.dat',ITS.sample_mean)
np.savetxt(f'ITS_bayesian.std.dat',ITS.sample_std)
