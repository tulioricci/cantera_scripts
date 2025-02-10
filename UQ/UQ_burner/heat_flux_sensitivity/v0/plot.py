import numpy as np
import matplotlib.pyplot as plt

color_list = ["red", "green", "blue", "orange", "magenta"]
phi_array = [0.55, 0.70, 0.85, 1.00, 1.15, 1.30]
perturb_array = [0.9, 1.1]

nImportantReactions = 3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

baseline = np.array([    
-4.422205389872788, 
-4.821337964752764,
-5.203500865090126,
-5.702739876609651,
-5.330154832546809,
-5.550276829711381])

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

nphi = len(phi_array)
nperturb = len(perturb_array)
nreactions = 77 # FIXME

reactionList = []

kk = 0
heatFlux = np.zeros((nphi,nperturb,nreactions))
for reaction in range(0,nreactions):

    jj = 0
    for perturb in perturb_array:
        filename = ("linearParameter_A/reaction_"
                     + str('%04d' % reaction)
                     + "_" + str('%0.2f' % perturb)
                     + ".dat")
        ii = 0
        try:
            f = open(filename)
            for phi in phi_array:

                readLine = f.readline()
                data = readLine.split(" ",4)
                flux = -float(data[1])
                heatFlux[ii,jj,kk] = flux
                
                if flux < 1.0:
                    heatFlux[ii,jj,kk] = np.nan
                
                ii += 1
            f.close()
        except:
            heatFlux[ii,jj,kk] = np.nan
        
        reactionList.append(data[-1].replace("\n",""))
        
        jj += 1

        k = kk % len(color_list)
        plt.plot(np.array(phi_array), heatFlux[:,0,kk], "--", color=color_list[k])
        plt.plot(np.array(phi_array), heatFlux[:,1,kk], color=color_list[k])
    kk = kk + 1

plt.plot(phi_array, -baseline, color="black")
plt.legend()
plt.ylabel("Heat flux")
plt.xlabel("phi")
plt.ylim(4.0, 6.5)
plt.savefig("sensitivity_linearParam_A.png",dpi=200)
plt.close()

important_reactions = []
for ii in range(len(phi_array)):
    important_reactions.append([])
    
    print("baseline = ", -baseline[ii])
    for jj in range(len(perturb_array)):
        idx = np.flip(np.argsort(heatFlux[ii,jj,:]))
        
        for kk in range(nImportantReactions):
            print(idx[kk], perturb_array[jj], reactionList[idx[kk]], heatFlux[ii,jj,idx[kk]])
            important_reactions[ii].append(idx[kk])
            
    print()
