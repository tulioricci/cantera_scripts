import numpy as np
import matplotlib.pyplot as plt
import os
import cantera as ct
import sys

plt.close("all")

volume_per_min = 25

coeff = 0.9

mech = "uiuc_20sp"

width = 0.010 # m

phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])

####################################

surf_temp = 300.0
burner_temp = 300.0

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
nspecies = gas0.n_species
del gas0

for spc_idx in range(nspecies):
    for parameter in ["diameter","well-depth","dipole","polarizability","rotational-relaxation"]:
    
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
