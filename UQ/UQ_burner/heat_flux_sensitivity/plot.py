import numpy as np
import matplotlib.pyplot as plt

phi_array = [0.55,0.70,0.85,1.00,1.15,1.30]
perturb_array = [0.5, 2.0]
nrxn = 78

nphi = len(phi_array)
nperturb = len(perturb_array)

baseline = np.array([    
-4.422205389872788, 
-4.821337964752764,
-5.203500865090126,
-5.702739876609651,
-5.330154832546809,
-5.550276829711381])

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"


for reaction in range(0,nrxn):
    heatFlux = np.zeros((nphi,nperturb))

    jj = 0
    for perturb in perturb_array:
        ii = 0
        filename = ("linearParameter_A/reaction_"
                     + str('%04d' % reaction)
                     + "_" + str('%0.2f' % perturb)
                     + ".dat")
        try:
            f = open(filename)
            for phi in phi_array:

                line1 = f.readline()
                line1 = line1.split(" ",4)
                flux = float(line1[1])
                heatFlux[ii,jj] =  flux
                
                if flux > -1.0:
                    heatFlux[ii,jj] = np.nan
                    print(phi, reaction, perturb, flux)
                
                ii += 1
            f.close()            
            jj += 1
        except:
            heatFlux[:,:] = np.nan

    plt.plot(phi_array, -heatFlux[:,0], color="purple")
    plt.plot(phi_array, -heatFlux[:,1], color="magenta")

plt.plot(phi_array, -baseline, color="black")
plt.legend()
plt.ylabel("Heat flux")
plt.xlabel("phi")
plt.ylim(4.0, 6.5)
plt.savefig("sensitivity_linearParam_A.png",dpi=200)
plt.close()
