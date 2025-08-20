import numpy as np
import matplotlib.pyplot as plt
import os
import cantera as ct
import sys

plt.close("all")

volume_per_min = 25

coeff = 1.01

mech = "uiuc_20sp"

width = 10 # mm

#phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])
phi_array = np.asarray([1.0])

####################################

def chebyshev_nodes(n, a, b):
    """Generates n Chebyshev nodes in the interval [a, b]."""
    k = np.arange(1, n + 1)
    nodes = 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k - 1) * np.pi / (2 * n))
    return nodes


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
slope = 0.03
curve = 0.03
loglevel = 0

if transp_model == 'mixture-averaged':
   enable_soret = False

r_int = 2.38*25.4/2000
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

gas0 = ct.Solution(mech + ".yaml")
nspecies = gas0.n_species
del gas0

nGrid = 640
nodes = np.hstack((1.0, chebyshev_nodes(nGrid, -1, 1), -1.0))[::-1]
nodes = 0.5*(nodes+1) * width/1000

for parameter in ["mu","kappa","diffO2-","diffN2-","diffCO2-","diffH2O-","diffC2H4-"]:
    for spc_idx in range(nspecies):
        
        gas = ct.Solution(mech + ".yaml")
        gas.transport_model = "mixture-averaged"

        x = np.zeros(gas.n_species,)
        x[gas.species_index('N2')] = 1.0
        gas.TPX = 273.15, 101325.0, x
        rho_ref = gas.density
        
        if parameter == "mu":
            mu_arr = gas.get_viscosity_polynomial(spc_idx)
            gas.set_viscosity_polynomial(spc_idx, mu_arr*1.01)        
        
        if parameter == "kappa":
            kappa_arr = gas.get_thermal_conductivity_polynomial(spc_idx)
            gas.set_thermal_conductivity_polynomial(spc_idx, kappa_arr*1.01)
            
        if parameter[:4] == "diff":
            dilutentIndx = gas.species_index(parameter[:-1].replace("diff",""))
            print(dilutentIndx, gas.species_name(dilutentIndx))
            diff_arr = gas.get_binary_diff_coeffs_polynomial(spc_idx, dilutentIndx)
            gas.set_binary_diff_coeffs_polynomial(spc_idx, dilutentIndx, diff_arr*1.01)
            gas.set_binary_diff_coeffs_polynomial(dilutentIndx, spc_idx, diff_arr*1.01)

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
                sim.soret_enabled = enable_soret
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
                    
                    print(spc_idx, gas.species_name(spc_idx), maxT, flux, parameter)
                    data = [spc_idx, maxT, flux, parameter]        
                
                except:
                    #print(spc_idx, gas.species_name(spc_idx), np.nan, np.nan, parameter)
                    #data = [spc_idx, np.nan, np.nan, parameter]
                    import sys
                    sys.exit()
                
                np.savetxt(result_file, data, fmt="%s")
                
        del gas
