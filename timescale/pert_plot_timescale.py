import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"] = 18
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.dpi"] = 200

import cantera as ct
gas0 = ct.Solution('mechanism.yaml')
species = gas0.species()
reactions = gas0.reactions()


import glob
nreactions = len(glob.glob("pert_timescale*"))

from random import randint
colors = []
for i in range(nreactions):
    colors.append('#%06X' % randint(0, 0xFFFFFF))

import os
for ii, filename in enumerate(sorted(glob.glob("pert_timescale*dat"))):

    if os.path.exists(filename[:-4] + '.png'):
        continue

    data = np.genfromtxt('ct_timescale.dat')
    temp = data[:,0]
    timescale = np.abs(data[:,:])
    nspecies = timescale.shape[1]
    for i in range(0,nspecies):
        plt.scatter(temp, timescale[:,i], marker='o',color='black')

    print(filename)
    data = np.genfromtxt(filename)
    temp = data[:,0]
    timescale = np.abs(data[:,:])
    for i in range(0,nspecies-1):
        plt.scatter(temp, timescale[:,i], marker='.',color=colors[ii])
    plt.scatter(temp, timescale[:,-1], marker='.',color=colors[ii], label=reactions[ii].equation)

    plt.plot([900,3100],[1e-7,1e-7], ':',color='black')
    plt.plot([900,3100],[1e-6,1e-6],'--',color='black')

    plt.ylim([1e-9,0.01])
    plt.xlim([900,3100])

    plt.legend()
    plt.yscale('log')
    plt.ylabel('Timescale')
    plt.xlabel('Temperature')
    plt.savefig(filename[:-4] + '.png',dpi=200)
#    plt.show()
    plt.close()

