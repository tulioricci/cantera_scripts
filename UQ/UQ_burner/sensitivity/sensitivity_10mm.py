import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import *

#plt.rc('text', usetex=True)
#plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"
plt.rcParams.update({'axes.labelsize': 20})

import cantera as ct
import pandas as pd
import sys

plt.close("all")

volume_per_min = 25
mdot_rate = 0.1693/25*volume_per_min

width = 0.010 # m

phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30, 2.0])

color_list = ["red", "green", "blue", "black", "orange", "magenta"]

####################################

surf_temp = 300.0
burner_temp = 300.0

mech = "uiuc_20sp.yaml"

transp_model = 'mixture-averaged'
#transp_model = 'multicomponent'
#enable_soret = True

os.system("mkdir -p csv_baseline")

##############

ratio = 2
slope = 0.1
curve = 0.2
loglevel = 0

if transp_model == 'mixture-averaged':
   enable_soret = False

r_int = 2.38*25.4/2000
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

gas = ct.Solution(mech)
gas.transport_model = transp_model

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kk = 0
for phi in phi_array:
    air = "O2:0.21,N2:0.79"
    fuel = "C2H4:1"
    gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

    gas.TP = burner_temp, 101325.0
    temp, rho_int, x = gas.TDX

    mass_unb = volume_per_min/np.sum(x[:])
    u_int = mass_unb*lmin_to_m3s/A_int
    rhoU_int = rho_int*u_int
    mdot_int = rho_int*u_int*A_int

    sim = ct.ImpingingJet(gas=gas, width=width)
    sim.inlet.mdot = mdot_rate
    sim.surface.T = surf_temp
    sim.soret_enabled = enable_soret
    sim.set_initial_guess(products='equil')
    sim.radiation_enabled = False

    sim.set_refine_criteria(ratio=4, slope=0.9, curve=0.9, prune=0.1)
    sim.solve(loglevel, refine_grid=True, auto=True)

    sim.set_refine_criteria(ratio=2, slope=0.5, curve=0.5, prune=0.1)
    sim.solve(loglevel, refine_grid=True, auto=True)

    sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=0.)
    sim.solve(loglevel, refine_grid=True, auto=False)

    # Create a dataframe to store sensitivity-analysis data
    sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))  

    # write the velocity, temperature, and mole fractions to a CSV file
    sim.write_csv('./csv_baseline/stagnation_flame_uiuc20sp' +
        '_m' + str('%4.2f' % volume_per_min) + 
        '_phi' + str('%4.2f' % phi) +
        '.csv', species="X", quiet=False)

#    Su0 = np.max(sim.T)
#    #dT = sim.T[-2] - sim.T[-1]
#    #dx = sim.grid[-2] - sim.grid[-1]
#    #Su0 = sim.thermal_conductivity[-1]*dT/dx/10000.0

#    dk = 1e-2

#    # Create an empty column to store the sensitivities data
#    sensitivities["baseCase"] = ""


#    for m in range(gas.n_reactions):
#        print(m, end="\r")
#        gas.set_multiplier(1.0) # reset all multipliers
#        gas.set_multiplier(1+dk, m) # perturb reaction m   
#        
#        # Always force loglevel=0 for this
#        # Make sure the grid is not refined, otherwise it won't strictly 
#        # be a small perturbation analysis
#        sim.solve(loglevel, refine_grid=False, auto=False)
#        
#        # The new flame speed
#        Su = np.max(sim.T)
#        #dT = sim.T[-2] - sim.T[-1]
#        #dx = sim.grid[-2] - sim.grid[-1]
#        #Su = sim.thermal_conductivity[-1]*dT/dx/10000.0
#        
#        sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)

#    # This step is essential, otherwise the mechanism will have been altered
#    gas.set_multiplier(1.0)

#    # Reaction mechanisms can contains thousands of elementary steps. Choose a
#    # threshold to see only the top few
#    threshold = 0.002

#    firstColumn = sensitivities.columns[0]

#    # For plotting, collect only those steps that are above the threshold
#    # Otherwise, the y-axis gets crowded and illegible
#    #sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
#    sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
#    sensitivitiesSubset = sensitivitiesSubset[:20]
#    indicesMeetingThreshold = sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
#    sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="",legend=None)
#    #print(sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False))

#    plt.gca().invert_yaxis()

#    plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');
#    plt.tight_layout()
#    #plt.xlim(-0.15,0.10)
#    plt.savefig('sensitivity_' + str('%03d' % (phi*100)) + '.png',dpi=300)
#    plt.close()

#    kk += 1
#    print('m =', u_int*A_int/lmin_to_m3s,', phi =',phi,
#          ', Tb =', sim.T[0], sim.T[-1], sim.velocity[0])

"""########################################################################"""

#cantera_uiuc_20sp = np.array(cantera_uiuc_20sp)

#plt.rcParams["axes.labelsize"] = 16
#plt.rcParams["xtick.labelsize"] = 14
#plt.rcParams["ytick.labelsize"] = 14
#plt.rcParams["figure.dpi"] = 200

#plt.close('all')
#fig = plt.figure(1,figsize=[6.4,5.2])
#ax1 = fig.add_subplot(111)
#ax2 = plt.twiny()
#ax3 = plt.twinx()
#ax1.set_position([0.15,0.12,0.8,0.84])
#ax2.set_position([0.15,0.12,0.8,0.84])
#ax3.set_position([0.15,0.12,0.8,0.84])
#ax1.tick_params(axis='both', labelsize=16, bottom=True, top=False, left=True, right=False)
#ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
#ax2.xaxis.set_minor_locator(AutoMinorLocator())
#ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
#ax3.yaxis.set_minor_locator(AutoMinorLocator())

#ax1.plot(splm170_exp[:,0], splm170_exp[:,1], lw=3.0, color='black', marker="o", markerfacecolor="white", markersize=8, label='Experiment')
#ax1.plot(cantera_uiuc_20sp[:,0], -cantera_uiuc_20sp[:,4], '-.', lw=3.0, color='forestgreen', marker='P', markersize=8, label='Cantera, 20sp (1D)')
#ax1.set_xlim(0.48, 1.32)
#ax2.set_xlim(0.48, 1.32)
#ax1.set_ylim(2.5, 5.0)
#ax3.set_ylim(2.5, 5.0)
#ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
#ax1.set_ylabel(r'$\bf{Heat \;  flux \; (W/cm^2)}$')
#ax1.legend(loc='upper left', bbox_to_anchor=(-0.01,1.01), fontsize=13)
##ax3.legend(loc='upper left', Bbox_to_anchor=(-0.01,0.79), fontsize=13)
#legend = ax1.get_legend()
#kk = 0
#for text in legend.get_texts():
#    text.set_color(ax1.get_lines()[kk].get_color())
#    kk = kk + 1
#    if kk == 7:
#        break
#plt.savefig('burner_heat_flux_m' + str('%4.2f' % volume_per_min) + '.png')
#plt.show()

