import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"] = 18
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.dpi"] = 200

filename = 'timescale'
data = np.genfromtxt(filename + '.dat')

temp = data[:,0]
timescale = np.abs(data[:,:])

nspecies = timescale.shape[1]
for i in range(0,nspecies):
    plt.scatter(temp, timescale[:,i], marker='*')

#species = ['H2', 'H', 'O2', 'O', 'OH', 'HO2', 'H2O2', 'H2O', 'N2']
#species = ['C2H4', 'O2', 'CO2', 'CO', 'H2O', 'H2', 'N2']

#plt.plot(temp, timescale[:,0], color='black', label=species[0], marker='*')
#plt.plot(temp, timescale[:,1], color='red', label=species[1])
#plt.plot(temp, timescale[:,2], color='green', label=species[2])
#plt.plot(temp, timescale[:,3], color='blue', label=species[3])
#plt.plot(temp, timescale[:,4], color='orange', label=species[4])
#plt.plot(temp, timescale[:,5], color='purple', label=species[5])
#plt.plot(temp, timescale[:,6], color='magenta', label=species[6])
#plt.plot(temp, timescale[:,7], color='cyan', label=species[7])
#plt.plot(temp, timescale[:,8], color='yellow', label=species[8])

plt.plot([900,3100],[1e-7,1e-7], ':',color='black')
plt.plot([900,3100],[5e-7,5e-7],'--',color='black')

#plt.ylim([1e-15,0.001])
plt.ylim([1e-9,0.01])
plt.xlim([900,3100.0])

#plt.legend()
plt.yscale('log')
plt.ylabel('Timescale')
plt.xlabel('Temperature')
plt.savefig(filename + '.png',dpi=200)
plt.show()

