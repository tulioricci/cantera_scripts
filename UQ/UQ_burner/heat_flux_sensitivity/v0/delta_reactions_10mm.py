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

width = 0.010 # m

phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])

####################################

surf_temp = 300.0
burner_temp = 300.0

mech = "wang99_51sp"

transp_model = 'mixture-averaged'
#transp_model = 'multicomponent'
#enable_soret = True

use_radiation = True

os.system("mkdir -p csv")
os.system("mkdir -p flux")

##############

ratio = 2
slope = 0.05
curve = 0.05
loglevel = 0

if transp_model == 'mixture-averaged':
   enable_soret = False

r_int = 2.38*25.4/2000
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

gas0 = ct.Solution(mech + ".yaml")
species = gas0.species()
reactions = gas0.reactions()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#selected_reacts = np.arange(0,10,1)
#p_or_m = ["-","+","+","-","+","-","-","+","-","+"]

#for i in range(0,10):
#    print(reactions[i], p_or_m[i])

#from itertools import chain, combinations

#def powerset(iterable):
#    s = list(iterable)  # allows duplicate elements
#    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

#list_of_combos = []
#list_of_combos.append([])
#list_of_combos.append((0, 1, 2, 3, 4, 5, 6, 7, 8, 9))


flux = []
maxT = []
idx = []
#for i, combo in enumerate(powerset(selected_reacts), 1):
#for i, combo in enumerate(list_of_combos):
    #print('combo #{}: {}'.format(i, combo))

for ireact, reaction in enumerate(reactions):
    #print(ireact, reactions[ireact])

    coeff = 1.1
    
    custom_reactions = [r for r in reactions]
    rxn_type = custom_reactions[ireact].reaction_type

    if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
        A_parameter = reactions[ireact].rate.input_data['rate-constant']['A']
    if rxn_type == "falloff-Troe":
        A_parameter = reactions[ireact].rate.input_data['low-P-rate-constant']['A']
    if np.abs(A_parameter) > 0.0:
        print(ireact, reactions[ireact].equation, A_parameter)
    else:
        continue

    if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(coeff*reactions[ireact].rate.input_data['rate-constant']['A'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
            third_body=custom_reactions[ireact].third_body)
            
    if rxn_type == "falloff-Troe":
    
        low = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['low-P-rate-constant']['A'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['b'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['Ea'])
    
        high = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['high-P-rate-constant']['A'],
                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['b'],
                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['Ea']) 
        falloff_coeffs = np.array([
                reactions[ireact].rate.input_data["Troe"]["A"],
                reactions[ireact].rate.input_data["Troe"]["T3"],
                reactions[ireact].rate.input_data["Troe"]["T1"],
                reactions[ireact].rate.input_data["Troe"]["T2"]])
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.TroeRate(high=high, low=low, falloff_coeffs=falloff_coeffs)
            )

    gas = ct.Solution(thermo='ideal-gas', kinetics='gas', 
                       species=species, reactions=custom_reactions)
    gas.transport_model = transp_model

    x = np.zeros(gas.n_species,)
    x[gas.species_index('N2')] = 1.0
    gas.TPX = 273.15, 101325.0, x
    rho_ref = gas.density
    
    for phi in phi_array:
        print("phi=", phi)

        result_file = ('./flux/stagnation_flame_' + mech +
                       '_m' + str('%4.2f' % volume_per_min) + 
                       '_phi' + str('%4.2f' % phi) +
                       '_S' + str('%04d' % ireact) + 
                       '_coeff' + str('%4.2f' % coeff) +                          
                       '.csv') 

        csv_file = ('./csv/stagnation_flame_' + mech +
                   '_m' + str('%4.2f' % volume_per_min) + 
                   '_phi' + str('%4.2f' % phi) +
                   '_S' + str('%04d' % ireact) +
                   '_coeff' + str('%4.2f' % coeff) +
                   '.csv')

        if os.path.exists(result_file) is False:

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

            try:
                sim.inlet.mdot = rhoU_ref
                
                sim.set_refine_criteria(ratio=4, slope=0.25, curve=0.25,
                                        prune=0.1)
                sim.solve(loglevel, refine_grid=True, auto=True)

                sim.set_refine_criteria(ratio=3, slope=0.15, curve=0.15,
                                        prune=0.05)
                sim.solve(loglevel, refine_grid=True, auto=True)

                sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve,
                                        prune=0.025)
                sim.solve(loglevel, refine_grid=True, auto=True)

                sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve,
                                        prune=0.025)
                sim.solve(loglevel, refine_grid=True, auto=True)

                dT = sim.T[-2] - sim.T[-1]
                dx = sim.grid[-2] - sim.grid[-1]
                kappa = sim.thermal_conductivity
                flux = kappa[-1]*dT/dx/10000.0
                maxT = max(sim.T)
    
                # write the velocity, temperature, and mole fractions to a CSV file
                sim.save(csv_file, basis="mole")
                
                if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
                    print(ireact, maxT, flux, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['rate-constant']['A'])
                    data = [ireact, maxT, flux, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['rate-constant']['A']]
    
                if rxn_type == "falloff-Troe":    
                    print(ireact, maxT, flux, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['low-P-rate-constant']['A'])
                    data = [ireact, maxT, flux, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['low-P-rate-constant']['A']]          

            except:
                if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
                    data = [ireact, np.nan, np.nan, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['rate-constant']['A']]          
                if rxn_type == "falloff-Troe":    
                    data = [ireact, np.nan, np.nan, custom_reactions[ireact].equation, custom_reactions[ireact].rate.input_data['low-P-rate-constant']['A']]          
                 
            np.savetxt(result_file, data, fmt="%s")
