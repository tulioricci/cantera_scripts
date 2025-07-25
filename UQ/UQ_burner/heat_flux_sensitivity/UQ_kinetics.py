import numpy as np
import os
import cantera as ct
import sys

#mech = "wang99_51sp"
mech = "uiuc_20sp"

####################################

#if mech == "wang99_51sp":
#    new_A_parameters = [(63.076845661346454, 30.172623109393093),  25.142106444743007,  10.770588040219511,  (56.204000710479825, 19.771347927429105),  23.43131603804862,  8.1886891244442,  32.86212973278313,  (83.0753849045796, 20.80022687808254),  10.833681189579274,  22.33070174670984]
#    new_b_parameters = [(-4.76, -0.63),  0.0,  1.228,  (-3.4, 1.62),  0.0,  2.0,  -1.0,  (-7.62, 0.454),  1.93,  0.1]
#    new_E_parameters = [(10208960.0, 1602472.0000000002),  60303992.00000001,  292880.0, (149781844.48000002, 155009752.48000002), 2510400.0, 10460000.000000002, 71128000.0, (29162480.0, 7614880.000000001), 54182800.0, 44350400.0]
#    list_of_reactions = [84, 1, 30, 151, 63, 224, 49, 219, 220, 102]

if mech == "uiuc_20sp":
    #new_lnA_parameters = [(6.30768457e+01, 3.01726231e+01), 3.09058991e+01, 3.28621297e+01, 1.07705880e+01, 6.38283073e+00, 2.51576477e+01, 3.83229660e+01, 9.64601125e+01, 3.04355929e+01, 3.14596625e+01]
    new_lnA_parameters = [30.905899130912033, 10.770588040219511, (63.076845661346454, 30.172623109393093), 32.86212973278313, 25.15764770195422, 22.751414084238696, 30.43559288402138, 25.03875194781205, 31.45966251241764, 24.124463218608565]    
#    new_b_parameters = [(-4.76000000e+00, -6.30000000e-01), -6.70700000e-01, -1.00000000e+00, 1.22800000e+00, 2.43300000e+00, 0.00000000e+00, -2.00000000e+00, -9.14700000e+00, -7.60000000e-01, -1.39000000e+00]
#    new_E_parameters = [(1.02089600e+07, 1.60247200e+06), 7.12995440e+07, 7.11280000e+07, 2.92880000e+05, 2.23852368e+08, 0.00000000e+00, 0.00000000e+00, 1.96229600e+08, 0.00000000e+00, 4.22584000e+06]

    list_of_reactions = [1, 18, 69, 25, 37, 26, 10, 14, 47, 44]
    

nreactions = len(list_of_reactions)
if len(new_lnA_parameters) != nreactions:
    print("new_lnA_parameters is wrong")
    sys.exit()
#if len(new_b_parameters) != nreactions:
#    print("new_b_parameters is wrong")
#    sys.exit()
#if len(new_E_parameters) != nreactions:
#    print("new_E_parameters is wrong")
#    sys.exit()

volume_per_min = 25
width = 0.010
phi_array = np.asarray([0.55, 0.70, 0.85, 1.0, 1.15, 1.30])
surf_temp = 300.0
burner_temp = 300.0

transp_model = 'mixture-averaged'
use_radiation = True
enable_soret = False

#os.system("mkdir -p csv")
os.system("mkdir -p flux")

##############

ratio = 2
slope = 0.03
curve = 0.03
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
    if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
        print(custom_reactions[ireact].input_data['rate-constant']['A'])
    if rxn_type == "falloff-Troe":
        print(custom_reactions[ireact].input_data['low-P-rate-constant']['A'])
        print(custom_reactions[ireact].input_data['high-P-rate-constant']['A'])
    
    #continue

    #if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
    if type(new_lnA_parameters[ii]) is float:
        new_A = np.exp(new_lnA_parameters[ii])
        new_b = custom_reactions[ireact].input_data['rate-constant']['b']
        new_E = custom_reactions[ireact].input_data['rate-constant']['Ea']        
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(new_A, new_b, new_E),
            third_body=custom_reactions[ireact].third_body)
            
    #if rxn_type == "falloff-Troe":
    if type(new_lnA_parameters[ii]) is tuple:
        new_A_low = np.exp(new_lnA_parameters[ii][0])
        new_b_low = custom_reactions[ireact].input_data['low-P-rate-constant']['b']
        new_E_low = custom_reactions[ireact].input_data['low-P-rate-constant']['Ea']
        low = ct.Arrhenius(new_A_low, new_b_low, new_E_low)

        new_A_high = np.exp(new_lnA_parameters[ii][1])
        new_b_high = custom_reactions[ireact].input_data['high-P-rate-constant']['b']
        new_E_high = custom_reactions[ireact].input_data['high-P-rate-constant']['Ea']        
        high = ct.Arrhenius(new_A_high, new_b_high, new_E_high)    

        falloff_coeffs = np.array([
                reactions[ireact].rate.input_data["Troe"]["A"],
                reactions[ireact].rate.input_data["Troe"]["T3"],
                reactions[ireact].rate.input_data["Troe"]["T1"],
                reactions[ireact].rate.input_data["Troe"]["T2"]])
                
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.TroeRate(high=high, low=low, falloff_coeffs=falloff_coeffs),
            third_body=custom_reactions[ireact].third_body
            )

gas = ct.Solution(thermo='ideal-gas', kinetics='gas', 
                   species=species, reactions=custom_reactions)
gas.transport_model = transp_model

x = np.zeros(gas.n_species,)
x[gas.species_index('N2')] = 1.0
gas.TPX = 273.15, 101325.0, x
rho_ref = gas.density

for phi in phi_array:
    print("phi =", phi)

    result_file = ('./flux/stagnation_flame_' + mech +
                   '_m' + str('%4.2f' % volume_per_min) + 
                   '_phi' + str('%4.2f' % phi) +
                   #'_S' + str('%04d' % UQ_sample) +                       
                   '.csv') 

    #stag_csv_file = ('./csv/stagnation_flame_' + mech +
    #           '_m' + str('%4.2f' % volume_per_min) + 
    #           '_phi' + str('%4.2f' % phi) +
    #           #'_S' + str('%04d' % UQ_sample) +
    #           '.csv')
               
    #flame_csv_file = ('./csv/adiabatic_flame_' + mech +
    #           '_m' + str('%4.2f' % volume_per_min) + 
    #           '_phi' + str('%4.2f' % phi) +
    #           #'_S' + str('%04d' % UQ_sample) +
    #           '.csv')               

    if True:
    #if os.path.exists(result_file) is False:

        air = "O2:0.21,N2:0.79"
        fuel = "C2H4:1"
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)
        
        #~~~ first, solve the 1D flame speed:
        f = ct.FreeFlame(gas, width=0.02)

        try:
            f.energy_enabled = False
            f.set_refine_criteria(ratio=4, slope=0.3, curve=0.3)
            f.solve(loglevel=loglevel, refine_grid=True, auto=True)

            f.energy_enabled = True
            f.set_refine_criteria(ratio=3, slope=0.15, curve=0.15)
            f.solve(loglevel=loglevel, refine_grid=True, auto=True)        

            f.set_refine_criteria(ratio=2, slope=0.05, curve=0.05)
            f.solve(loglevel=loglevel, refine_grid=True, auto=True)

            flame_speed = f.velocity[0]
            #print(flame_speed)
            #f.save(flame_csv_file, basis="mole")
        except:
            flame_speed = np.nan

        #~~~ then, solve the stagnation flow
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)
        gas.TP = burner_temp, 101325.0
        #rho_int = gas.density

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
            
            sim.set_refine_criteria(ratio=4, slope=0.15, curve=0.15,
                                    prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)

            sim.set_refine_criteria(ratio=3, slope=0.10, curve=0.10,
                                    prune=0.05)
            sim.solve(loglevel, refine_grid=True, auto=True)

            sim.set_refine_criteria(ratio=ratio, slope=0.05, curve=0.05,
                                    prune=0.025)
            sim.solve(loglevel, refine_grid=True, auto=True)

            sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve,
                                    prune=0.025)
            sim.solve(loglevel, refine_grid=True, auto=True)

            dT = sim.T[-2] - sim.T[-1]
            dx = sim.grid[-2] - sim.grid[-1]
            kappa = sim.thermal_conductivity
            flux = kappa[-1]*dT/dx/10000.0
            
            idx = np.argmin(np.abs(sim.grid - 0.005))
            temperature = sim.T[idx]
            X_CO = sim.X[gas.species_index("CO"),idx]
            X_CO2 = sim.X[gas.species_index("CO2"),idx]
            
            #print(temperature)
            #print(X_CO)
            #print(X_CO2)
            #print(flux)
            #~~ write the velocity, temperature, and mole fractions to a CSV file
            #sim.save(stag_csv_file, basis="mole")
            
        except:
            temperature = np.nan
            X_CO = np.nan
            X_CO2 = np.nan
            flux = np.nan
        
        data = [flame_speed, temperature, X_CO, X_CO2, flux]
        np.savetxt(result_file, data, fmt="%s")
        
        #sim.save("debug" + str('%4.2f' % phi) + ".csv", basis="mole", overwrite=True)
