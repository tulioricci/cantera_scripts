import numpy as np
import os
import cantera as ct
import sys

new_A_parameters = [10, 1234, 666, 3.141592]
list_of_reactions = [1, 3, 4, 10]
UQ_sample = 1  # identifier of UQ sample space

####################################

nreactions = len(new_A_parameters)

volume_per_min = 25
width = 0.010
phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])
surf_temp = 300.0
burner_temp = 300.0
mech = "uiuc_20sp"

transp_model = 'mixture-averaged'
use_radiation = True
enable_soret = False

os.system("mkdir -p csv")
os.system("mkdir -p flux")

##############

ratio = 2
slope = 0.05
curve = 0.05
loglevel = 0

r_int = 2.38*25.4/2000
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

gas0 = ct.Solution(mech + ".yaml")
species = gas0.species()
reactions = gas0.reactions()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flux = []
maxT = []
idx = []

custom_reactions = [r for r in reactions]

# modify all the desired reactions:
for ii, _ireact in enumerate(list_of_reactions):
    ireact = _ireact - 1 # convert to 0-index
    
    rxn_type = custom_reactions[ireact].reaction_type

    if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
        A_parameter = reactions[ireact].rate.input_data['rate-constant']['A']
    if rxn_type == "falloff-Troe":
        A_parameter = reactions[ireact].rate.input_data['low-P-rate-constant']['A']
    if np.abs(A_parameter) > 0.0:
        print(ireact, reactions[ireact].equation, A_parameter, new_A_parameters[ii])

    if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(new_A_parameters[ii],
                             #coeff*reactions[ireact].rate.input_data['rate-constant']['A'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
            third_body=custom_reactions[ireact].third_body)
            
    if rxn_type == "falloff-Troe":

        low = ct.Arrhenius(new_A_parameters[ii],
                           #coeff*reactions[ireact].rate.input_data['low-P-rate-constant']['A'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['b'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['Ea'])

        high = ct.Arrhenius(new_A_parameters[ii],
                            #coeff*reactions[ireact].rate.input_data['high-P-rate-constant']['A'],
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
                   '_S' + str('%04d' % UQ_sample) +                       
                   '.csv') 

    csv_file = ('./csv/stagnation_flame_' + mech +
               '_m' + str('%4.2f' % volume_per_min) + 
               '_phi' + str('%4.2f' % phi) +
               '_S' + str('%04d' % UQ_sample) +
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
            
            data = [ireact, maxT, flux]       

        except:
            data = [ireact, np.nan, np.nan]         
             
        np.savetxt(result_file, data, fmt="%s")
