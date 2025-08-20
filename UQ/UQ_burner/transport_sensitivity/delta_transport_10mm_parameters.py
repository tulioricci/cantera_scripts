import numpy as np
import matplotlib.pyplot as plt
import os
import cantera as ct
import sys


def chebyshev_nodes(n, a, b):
    """Generates n Chebyshev nodes in the interval [a, b]."""
    k = np.arange(1, n + 1)
    nodes = 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k - 1) * np.pi / (2 * n))
    return nodes

import getopt
arg_list = sys.argv[1:]
opts, args = getopt.getopt(arg_list,":",["parameter="])
for opt, arg in opts:
    if opt in ("--parameter"):
        parameter = arg
        print ('Parameter = ', parameter)


plt.close("all")

volume_per_min = 25

coeff = 0.99

mech = "uiuc_20sp"

width = 10 # m

phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])

####################################

surf_temp = 300.0
burner_temp = 300.0

transp_model = 'mixture-averaged'

use_radiation = True

os.system("mkdir -p csv")
os.system("mkdir -p flux")

##############

loglevel = 0

r_int = 2.38*25.4/2000
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

gas0 = ct.Solution(mech + ".yaml")
nspecies = gas0.n_species
del gas0

for spc_idx in range(nspecies):

    gas0 = ct.Solution(mech + ".yaml")
    gas0.transport_model = "mixture-averaged"
    species = gas0.species()
    reactions = gas0.reactions()
    
    newSpecies = []
    for spc in species:
        newSpecies.append(spc)
    print(parameter)
    print(species[spc_idx].input_data["transport"])     

    if parameter == "diameter": 
        old_parameter = species[spc_idx].transport.diameter
        newSpecies[spc_idx].transport.diameter = species[spc_idx].transport.diameter*coeff
        new_parameter = newSpecies[spc_idx].transport.diameter
    if parameter == "well-depth": 
        old_parameter = species[spc_idx].transport.well_depth
        newSpecies[spc_idx].transport.well_depth = species[spc_idx].transport.well_depth*coeff
        new_parameter = newSpecies[spc_idx].transport.well_depth
    if parameter == "rotational-relaxation":
        old_parameter = species[spc_idx].transport.rotational_relaxation 
        newSpecies[spc_idx].transport.rotational_relaxation = species[spc_idx].transport.rotational_relaxation*coeff
        new_parameter = newSpecies[spc_idx].transport.rotational_relaxation
        
    if parameter == "dipole":
        old_parameter = species[spc_idx].transport.dipole 
        newSpecies[spc_idx].transport.dipole = species[spc_idx].transport.dipole*coeff
        new_parameter = newSpecies[spc_idx].transport.dipole
    if parameter == "polarizability":
        old_parameter = species[spc_idx].transport.polarizability 
        newSpecies[spc_idx].transport.polarizability = species[spc_idx].transport.polarizability*coeff
        new_parameter = newSpecies[spc_idx].transport.polarizability

    gas = ct.Solution(thermo='ideal-gas', kinetics='gas', 
                      species=newSpecies, reactions=reactions)
    gas.transport_model = "mixture-averaged"

    x = np.zeros(gas.n_species,)
    x[gas.species_index('N2')] = 1.0
    gas.TPX = 273.15, 101325.0, x
    rho_ref = gas.density

    nGrid = 640
    nodes = np.hstack((1.0, chebyshev_nodes(nGrid, -1, 1), -1.0))[::-1]
    nodes = 0.5*(nodes+1) * width/1000
    
    for phi in phi_array:
        print("phi=", phi)

        result_file = ('./flux/stagnation_flame_' + mech +
                       '_m' + str('%4.2f' % volume_per_min) + 
                       '_phi' + str('%4.2f' % phi) +
                       '_' + parameter + gas.species_name(spc_idx) + 
                       '_coeff' + str('%4.2f' % coeff) +                          
                       '.dat')

        csv_file = ('./csv/stagnation_flame_' + mech +
                       '_m' + str('%4.2f' % volume_per_min) + 
                       '_phi' + str('%4.2f' % phi) +
                       '_' + parameter + gas.species_name(spc_idx) + 
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

            sim = ct.ImpingingJet(gas=gas, grid=nodes)
            sim.inlet.T = burner_temp
            sim.surface.T = surf_temp
            sim.set_initial_guess(products='equil')
            sim.radiation_enabled = use_radiation

            try:

                sim.inlet.mdot = rhoU_ref
                
                sim.solve(loglevel, refine_grid=False, auto=False)
                os.system("date")

                dT = sim.T[-2] - sim.T[-1]
                dx = sim.grid[-2] - sim.grid[-1]
                kappa = sim.thermal_conductivity
                flux = kappa[-1]*dT/dx/10000.0
                maxT = max(sim.T)

                # write the velocity, temperature, and mole fractions to a CSV file
                sim.save(csv_file, basis="mole")
                
                print(spc_idx, maxT, flux, old_parameter, new_parameter, parameter)
                data = [spc_idx, maxT, flux, old_parameter, new_parameter, parameter]        
            
            except:
                print(spc_idx, np.nan, np.nan, old_parameter, new_parameter, parameter)
                data = [spc_idx, np.nan, np.nan, old_parameter, new_parameter, parameter]      
            
            np.savetxt(result_file, data, fmt="%s")
            
    del gas, gas0, species, newSpecies
