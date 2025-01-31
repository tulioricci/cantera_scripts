import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import *
import cantera as ct
import sys
from statistics import median

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

skip = False

species_output = "X"

volume_per_min = 25

burner_temp = 300.0
surf_temp = 300.0

width = 0.010  # m

use_radiation = True

#transp_model = 'unity-Lewis-number'
#enable_soret = False

transp_model = 'mixture-averaged'
enable_soret = False

#transp_model = 'multicomponent'
##enable_soret = False
#enable_soret = True

####################################

path = "csv_10mm_" + species_output + "_" + transp_model
if enable_soret:
    path = "csv_10mm_" + species_output + "_" + transp_model + "_soret"

if use_radiation:
    path = path + "_radiation"

os.system("mkdir -p " + path)

phi_array = np.array([0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 2.0])

ratio = 2
slope = 0.05
curve = 0.05
loglevel = 0

r_int = 2.38*25.4/2000 #radius, actually
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

####################################

for mechanism in ["uiuc_20sp"]:

    mech_data = []
    for phi in phi_array:

        filename = (
            path + '/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        print("")
        gas = ct.Solution(mechanism + ".yaml")
        gas.transport_model = transp_model

        x = np.zeros(gas.n_species,)
        x[gas.species_index('N2')] = 1.0
        gas.TPX = 273.15, 101325.0, x
        rho_ref = gas.density

        air = "O2:0.21,N2:0.79"
        fuel = "C2H4:1"
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

        gas.TP = burner_temp, 101325.0
        _, rho_int, x = gas.TDX

        u_ref = volume_per_min*lmin_to_m3s/A_int
        rhoU_ref = rho_ref*u_ref

        sim = ct.ImpingingJet(gas=gas, width=width)
        sim.inlet.T = burner_temp
        sim.surface.T = surf_temp
        sim.soret_enabled = enable_soret
        sim.set_initial_guess(products='equil')
        sim.radiation_enabled = use_radiation

        flag = False
        while flag is False:

            sim.inlet.mdot = rhoU_ref
            
            print(rhoU_ref)
            sim.set_refine_criteria(ratio=4, slope=0.25, curve=0.25,
                                    prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            os.system("date")

            sim.set_refine_criteria(ratio=3, slope=0.15, curve=0.15,
                                    prune=0.05)
            sim.solve(loglevel, refine_grid=True, auto=True)
            os.system("date")

            sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve,
                                    prune=0.025)
            sim.solve(loglevel, refine_grid=True, auto=True)
            os.system("date")

            sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve,
                                    prune=0.025)
            sim.solve(loglevel, refine_grid=True, auto=True)
            os.system("date")
           
            flag = True

        dT = sim.T[-2] - sim.T[-1]
        dx = sim.grid[-2] - sim.grid[-1]
        kappa = sim.thermal_conductivity
        flux = kappa[-1]*dT/dx/10000.0
        maxT = max(sim.T)

        # write the velocity, temperature, and mole fractions to a CSV file
        #sim.save(csv_file, basis="mole")
         
        print(maxT, flux)
        data = [maxT, flux]          
        
        result_file = ('./flux/stagnation_flame_uiuc20sp' +
                       '_m' + str('%4.2f' % volume_per_min) + 
                       '_phi' + str('%4.2f' % phi) +
                       '_baseline.csv')         
        np.savetxt(result_file, data, fmt="%s")
    
"""########################################################################"""

