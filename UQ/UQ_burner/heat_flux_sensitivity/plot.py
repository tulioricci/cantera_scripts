import numpy as np
import matplotlib.pyplot as plt

phi_array = np.array([0.55,0.70,0.85,1.00,1.15,1.30])

data = np.zeros((6,5))

f = open("reaction_0000.dat")

ii = 0
for phi in phi_array:

    jj = 0
    for pertub in [0.5, 0.9, 1.1, 1.5, 2.0]:

        line1 = f.readline()
        line2 = f.readline()
        print(line1, float(line2))
        
        data[ii,jj] =  float(line2)
        
        jj += 1 
    
    ii += 1
    
baseline = np.array([    
-4.422205389872788, 
-4.821337964752764,
-5.203500865090126,
-5.702739876609651,
-5.330154832546809,
-5.550276829711381])


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
    
plt.plot(phi_array, -baseline, color="black")
plt.plot(phi_array, -data[:,0], color="purple", label="0.5")
plt.plot(phi_array, -data[:,1], color="firebrick", label="0.9")
plt.plot(phi_array, -data[:,2], color="red", label="1.1")
plt.plot(phi_array, -data[:,3], color="orange", label="1.5")
plt.plot(phi_array, -data[:,4], color="magenta", label="2.0")
plt.legend()
plt.ylabel("Heat flux")
plt.xlabel("phi")
plt.savefig("sensitivity_reaction0000.png",dpi=200)
plt.close()
