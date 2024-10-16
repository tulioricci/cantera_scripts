"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0

Sensitivity analysis
https://cantera.org/examples/jupyter/flames/flame_speed_with_sensitivity_analysis.ipynb.html
"""

import cantera as ct
import matplotlib
import sys
import numpy as np

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
width = 0.02  # m
loglevel = 0  # amount of diagnostic output (0 to 8)

import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.dpi"] = 120


#def compute_grad(q, x):

#    N = q.shape[0]
#    grad = np.zeros((N))

#    grad[ 0] = (q[ 1]-q[ 0])/(x[ 1]-x[ 0])
#    for i in range(1,N-1):
#        grad[i] = (q[i+1]-q[i-1])/(x[i+1]-x[i-1])
#    grad[-1] = (q[-1]-q[-2])/(x[-1]-x[-2])

#    idx_xloc = grad.argmax(axis=0)
#    xloc = x[idx_xloc]

#    return grad, xloc

H2 = 0.0
for phi in np.arange(0.85,1.45,0.15):
    for pressure in [1.0]:

        jj = 0
        for transport in ["mix"]:

            plt.close('all')
            fig, ax = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1, 2]})

            p = pressure * ct.one_atm  # pressure
            for mechanism in ['aramco']:

                rxnmech = mechanism.replace("../","")
                #print(rxnmech, H2, pressure, phi)

                # Solution object used to compute mixture properties
                # set to the state of the upstream fuel-air mixture
                gas = ct.Solution(mechanism + '.yaml')

                #https://cantera.org/examples/python/thermo/equivalenceRatio.py.html
                air = "O2:0.21,N2:0.79"
                fuel = "C2H4:1"
                gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

                gas.TP = Tin, p

                # Set up flame object
                f = ct.FreeFlame(gas, width=width)

                # Solve with mixture-averaged transport model
                if transport == "mix":
                    f.transport_model = 'mixture-averaged'
                if transport == "Le1":
                    f.transport_model = 'unity-Lewis-number' 

                f.energy_enabled = False
                f.set_refine_criteria(ratio=4, slope=0.5, curve=0.5)
                f.solve(loglevel=loglevel, refine_grid=True, auto=True)
                print("coarse ok")

                f.energy_enabled = True
                f.set_refine_criteria(ratio=2, slope=0.15, curve=0.15)
                f.solve(loglevel=loglevel, refine_grid=True, auto=True)

                print('flamespeed = {0:7f} m/s'.format(f.velocity[0]))

                # write the velocity, temperature, density, and mole fractions
                filename = (rxnmech + '_phi' + str('%4.2f' % phi) + 
                    '_p' + str('%4.2f' % pressure) + 
                    '_E:H' + str('%3.1f' % H2) + '_' + transport + '.csv')
                f.write_csv('./csv/adiabatic_flame_' + filename,
                            species="X", quiet=False)

#                nx = f.grid.shape[0]
#                nspecies = gas.n_species
#                transport_data = np.zeros((nx, 1+1+1+1+nspecies))
#                transport_data[:,0] = f.grid
#                transport_data[:,1] = f.T
#                transport_data[:,2] = f.viscosity
#                transport_data[:,3] = f.thermal_conductivity
#                transport_data[:,4:] = f.mix_diff_coeffs.T

#                header = ["x","T","mu","kappa"]
#                for i in range(nspecies):
#                    header.append("diff_" + gas.species_name(i))
#                
#                np.savetxt('./transport/flame1D_transp_' + filename,
#                           transport_data, delimiter=",",
#                           header=', '.join(str(x) for x in header))


