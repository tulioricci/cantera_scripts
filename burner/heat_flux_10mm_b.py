import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import *
import cantera as ct
import sys
from statistics import median

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

species_output = "X"

volume_per_min = 25
mdot_rate = 0.1693

width = 0.010  # m

phi_array = np.array([0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 2.0])

####################################

ratio = 2
slope = 0.1
curve = 0.2
loglevel = 0

r_int = 2.38*25.4/2000 #radius, actually
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

os.system("mkdir -p csv_10mm")

mechanism_list = ["uiuc_20sp", "wang99_22sp", "wang99_22sp_b", "wang99_27sp", "wang99_reduced_mod",
                  "wang99_32sp", "wang99_34sp", "wang99_36sp",
                  "wang99_51sp", "wang99_75sp"]

for mechanism in mechanism_list:
    for phi in phi_array:

        filename = (
            './csv_10mm/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        if os.path.exists(filename):
            print(filename)

        else:
            gas = ct.Solution(mechanism + ".yaml")
            gas.transport_model = 'mixture-averaged'

            air = "O2:0.21,N2:0.79"
            fuel = "C2H4:1"
            gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

            gas.TP = 300.0, 101325.0
            temp, rho_int, x = gas.TDX

            mass_unb = volume_per_min/np.sum(x[:])
            r_int = 2.38*25.4/2000 #radius, actually
            A_int = np.pi*r_int**2
            lmin_to_m3s = 1.66667e-5
            u_int = mass_unb*lmin_to_m3s/A_int
            rhoU_int = rho_int*u_int
            mdot_int = rho_int*u_int*A_int

            try:
                sim = ct.ImpingingJet(gas=gas, width=width)
                sim.inlet.mdot = mdot_rate*1.5
                sim.surface.T = 300.0
                sim.set_initial_guess(products='equil')

                sim.inlet.mdot = mdot_rate*1.0
                sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=0.0)
                sim.solve(loglevel, refine_grid=True, auto=True)
                print(sim.T.shape)

                # write the velocity, temperature, and mole fractions to a CSV file
                sim.write_csv('./csv_10mm/stagnation_flame_' + mechanism +
                    '_m' + str('%4.2f' % volume_per_min) + 
                    '_phi' + str('%4.2f' % phi) + 
                    '.csv', species=species_output, quiet=False)

                print('m =', u_int*A_int/lmin_to_m3s,', phi =',phi, ', Tb =', sim.T[0], sim.T[-1], sim.velocity[0])

            except:
                sys.exit()

