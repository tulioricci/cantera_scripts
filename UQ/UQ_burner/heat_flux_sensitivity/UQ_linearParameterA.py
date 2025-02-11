import numpy as np
import os
import cantera as ct
import sys

new_A_parameters = [(6.30768457e+01, 3.01726231e+01), 3.09058991e+01, 3.28621297e+01, 1.07705880e+01, 6.38283073e+00, 2.51576477e+01, 3.83229660e+01, 9.64601125e+01, 3.04355929e+01, 3.14596625e+01]
new_b_parameters = [(-4.76000000e+00, -6.30000000e-01), -6.70700000e-01, -1.00000000e+00, 1.22800000e+00, 2.43300000e+00, 0.00000000e+00, -2.00000000e+00, -9.14700000e+00, -7.60000000e-01, -1.39000000e+00]
new_E_parameters = [(1.02089600e+07, 1.60247200e+06), 7.12995440e+07, 7.11280000e+07, 2.92880000e+05, 2.23852368e+08, 0.00000000e+00, 0.00000000e+00, 1.96229600e+08, 0.00000000e+00, 4.22584000e+06]

list_of_reactions = [69, 1, 25, 18, 12, 37, 5, 50, 10, 47]

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
    print(custom_reactions[ireact].equation)

    #if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
    if type(new_A_parameters[ii]) is float:
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(new_A_parameters[ii], new_b_parameters[ii], new_E_parameters[ii]),
            third_body=custom_reactions[ireact].third_body)
            
    #if rxn_type == "falloff-Troe":
    if type(new_A_parameters[ii]) is tuple:
        low = ct.Arrhenius(new_A_parameters[ii][0], new_b_parameters[ii][0], new_E_parameters[ii][0])
        high = ct.Arrhenius(new_A_parameters[ii][1], new_b_parameters[ii][1], new_E_parameters[ii][1])

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
                   #'_S' + str('%04d' % UQ_sample) +                       
                   '.csv') 

    csv_file = ('./csv/stagnation_flame_' + mech +
               '_m' + str('%4.2f' % volume_per_min) + 
               '_phi' + str('%4.2f' % phi) +
               #'_S' + str('%04d' % UQ_sample) +
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
            
            data = [maxT, flux]

        except:
            data = [np.nan, np.nan]
             
        np.savetxt(result_file, data, fmt="%s")
