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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for phi in phi_array:
    for index in [1024]:
        for p_or_m in ["p", "m"]:

            mech = "yaml_" + p_or_m + "/uiuc_20sp_UQ_" + str("%04d" % index) + ".yaml"

            gas = ct.Solution(mech)
            gas.transport_model = transp_model

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

            # write the velocity, temperature, and mole fractions to a CSV file
            sim.write_csv('./csv/stagnation_flame_uiuc20sp' +
                '_m' + str('%4.2f' % volume_per_min) + 
                '_phi' + str('%4.2f' % phi) +
                '_S' + str('%04d' % index) + p_or_m +
                '.csv', species="X", quiet=False)

