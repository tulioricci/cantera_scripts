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

#phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30, 2.0])
phi_array = np.asarray([0.70, 1.0, 1.30, 2.0])

color_list = ["red", "green", "blue", "black", "orange", "magenta"]

####################################

surf_temp = 300.0
burner_temp = 300.0

mech = "uiuc_20sp_UQ.yaml"

transp_model = 'mixture-averaged'
#transp_model = 'multicomponent'
#enable_soret = True

os.system("mkdir -p csv")

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

gas0 = ct.Solution(mech)
species = gas0.species()
reactions = gas0.reactions()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

selected_reacts = np.arange(0,10,1)
p_or_m = ["-","+","+","-","+","-","-","+","-","+"]

for i in range(0,10):
    print(reactions[i], p_or_m[i])

from itertools import chain, combinations

def powerset(iterable):
    s = list(iterable)  # allows duplicate elements
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

list_of_combos = []
list_of_combos.append([])
list_of_combos.append((0, 1, 2, 3, 4, 5, 6, 7, 8, 9))

for phi in phi_array:

    maxT = []
    idx = []
    #for i, combo in enumerate(powerset(selected_reacts), 1):
    for i, combo in enumerate(list_of_combos):
        print('combo #{}: {}'.format(i, combo))

        factor = 1.1
        custom_reactions = [r for r in reactions]
        for ii in range(0,len(combo)):

            ireact = combo[ii]

            if p_or_m[ireact] == "+":
                coeff = 1.0/factor
            if p_or_m[ireact] == "-":
                coeff = 1.0*factor

            custom_reactions[ireact] = ct.Reaction(
                reactions[ireact].reactants,
                reactions[ireact].products,
                ct.ArrheniusRate(coeff*reactions[ireact].rate.input_data['rate-constant']['A'],
                                 1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                                 1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
                third_body=custom_reactions[ireact].third_body)

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas', 
                           species=species, reactions=custom_reactions)
        gas.transport_model = transp_model

        air = "O2:0.21,N2:0.79"
        fuel = "C2H4:1"
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

        gas.TP = burner_temp, 101325.0
        temp, rho_int, x = gas.TDX

        gas.write_yaml("uiuc_20sp_UQ_" + str("%04d" % i) + ".yaml",
                       units={"length": "cm", "quantity": "mol", "activation-energy": "cal/mol"})

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

        maxT.append(max(sim.T))
        idx.append(i-1)

        # write the velocity, temperature, and mole fractions to a CSV file
        sim.write_csv('./csv/stagnation_flame_uiuc20sp' +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) +
            '_S' + str('%04d' % i) +
            '.csv', species="X", quiet=False)

    print(np.argmax(maxT), max(maxT))
