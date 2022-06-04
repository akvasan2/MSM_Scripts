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
import seaborn as sns

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"

######################Loading data#################
###################################################################

data = []
pwd = '../Distance_for_MSM/'

for i in range(1,10):
    line_PROA = np.loadtxt(f'{pwd}/Run_{i}-PROA.dat')
    line_PROB = np.loadtxt(f'{pwd}/Run_{i}-PROB.dat')
    line_PROC = np.loadtxt(f'{pwd}/Run_{i}-PROC.dat')
    
    data.append(line_PROA)
    data.append(line_PROB)
    data.append(line_PROC)

data_concatenated = np.concatenate(data)

######################Converting data to tica #################
###################################################################

tica = pyemma.coordinates.tica(data,lag=100,var_cutoff=1)
tica_output = tica.get_output()
tica_concatenated = np.concatenate(tica.get_output())
tica_eigvals = tica.eigenvalues

###### Plotting each histogram of tica ########

for i in range(0,tica_concatenated.shape[1]):
    sns.kdeplot(tica_concatenated[:,i],shade=True,alpha=0.4)#,label=f'State {i+1}')
#    plt.hist(tica_concatenated[:,i],density=True,bins=15)
    plt.title(f'IC {i+1}')
    plt.ylim(0,1.2)
    plt.yticks(np.arange(0,1.1,0.5))
    plt.xlim(-4,4)
    plt.savefig(f'Images/IC_{i}_hist.png',bbox_inches='tight')
    plt.close()

