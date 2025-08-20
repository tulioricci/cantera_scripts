import cantera as ct
import numpy as np
#from functools import partial
import numpy.linalg as la
import sys

# get reaction coeffs
# https://groups.google.com/g/cantera-users/c/ID-oraSnUrU


np.set_printoptions(precision=8, linewidth=256)

# Cantera solution object
gas0 = ct.Solution('mechanism.yaml')
species = gas0.species()
reactions = gas0.reactions()

coeff = 1.10

perturb_species = True
perturb_reaction = False


if perturb_species:

    for spc_idx in range(0, len(species)):
        perturbed_species = gas0.species_name(spc_idx)
        print(spc_idx, perturbed_species)

        counter = 0
        custom_reactions = [r for r in reactions]

        for ireact in range(0,len(reactions)):
        
            ii = 0
            temp_array = np.linspace(1000,2800,29)
            temp_history = np.zeros((temp_array.shape[0],gas0.n_species+1))
            pressure = 101325.0

            rxn_type = custom_reactions[ireact].reaction_type

            change = False
            if perturbed_species in reactions[ireact].reactants:
                #print(reactions[ireact])
                counter += 1
                change = True

            if perturbed_species in reactions[ireact].products:
                #print(reactions[ireact])
                counter += 1
                change = True
                
            if change:
            
                if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
                    custom_reactions[ireact] = ct.Reaction(
                        reactions[ireact].reactants,
                        reactions[ireact].products,
                        ct.ArrheniusRate(coeff*reactions[ireact].rate.input_data['rate-constant']['A'],
                                         1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                                         1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
                        third_body=custom_reactions[ireact].third_body)

                    #print(reactions[ireact].rate.input_data['rate-constant']['A'])
                    #print(custom_reactions[ireact].rate.input_data['rate-constant']['A'])

                if rxn_type == "falloff-Troe":

                    low = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['low-P-rate-constant']['A'],
                                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['b'],
                                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['Ea'])

                    high = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['high-P-rate-constant']['A'],
                                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['b'],
                                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['Ea'])
                    falloff_coeffs = np.array([
                            reactions[ireact].rate.input_data["Troe"]["A"],
                            reactions[ireact].rate.input_data["Troe"]["T3"],
                            reactions[ireact].rate.input_data["Troe"]["T1"],
                            reactions[ireact].rate.input_data["Troe"]["T2"]])
                    custom_reactions[ireact] = ct.Reaction(
                        reactants=reactions[ireact].reactants,
                        products=reactions[ireact].products,
                        rate=ct.TroeRate(low=low, high=high, falloff_coeffs=falloff_coeffs)
                        )

                    #print(reactions[ireact].rate.input_data['low-P-rate-constant']['A'])
                    #print(custom_reactions[ireact].rate.input_data['low-P-rate-constant']['A'])     
                    
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=custom_reactions)

        print("counter =", counter)
        import matplotlib.pyplot as plt
        plt.close("all")
        
        data = np.genfromtxt('timescale.dat')
        temp = data[:,0]
        timescale = np.abs(data[:,:])
        nspecies = timescale.shape[1]
        for i in range(0,51):
            plt.scatter(temp, timescale[:,i], marker='o',color='black')    
        
        from random import randint
        colors = []
        for i in range(100):
            colors.append('#%06X' % randint(0, 0xFFFFFF))    
        for maxT in temp_array:

            # mass ratios
            air = "O2:1.0,N2:3.76"
            fuel = "C2H4:1"
            gas.set_equivalence_ratio(phi=1.0, fuel=fuel, oxidizer=air)

            gas.TP = 900.0, 101325.0

            # Define a reactor
            #r = ct.IdealGasConstPressureReactor(gas, name="Batch Reactor")
            r = ct.IdealGasReactor(gas, name="Batch Reactor")
            net = ct.ReactorNet([r])
            T = r.T
            while T < maxT:
                net.step()
                T = r.T

            jac_ct = gas.net_production_rates_ddX
            jac_matrix = jac_ct/gas.molecular_weights*gas.mean_molecular_weight

            w = np.real(la.eigvals(jac_matrix))

            timescale = np.nan_to_num(1.0/w, copy=True, nan=1.0, posinf=None, neginf=None)

            #temp_history[ii, 0] = T
            #temp_history[ii,1:] = timescale
            #ii = ii + 1

            for i in range(0,51-1):
                plt.scatter(T, np.abs(timescale[i]), marker='*',color=colors[i])
                
        plt.yscale('log')
        plt.ylim([1e-9,0.01])
        plt.savefig("pert_" + perturbed_species + ".png")
        plt.close()

        #np.savetxt(f'pert_timescale_species' + "CH2star" + '.dat', temp_history[:ii])                
                
                
                
                

if perturb_reaction:

    for ireact in range(0,len(reactions)):
        
        ii = 0
        temp_array = np.linspace(1000,2800,11)
        temp_history = np.zeros((temp_array.shape[0],gas0.n_species+1))
        pressure = 101325.0

        custom_reactions = [r for r in reactions]
        rxn_type = custom_reactions[ireact].reaction_type

        rxn_type = custom_reactions[ireact].reaction_type
        #print(ireact, reactions[ireact].equation, reactions[ireact].reaction_type)

        if rxn_type == "Arrhenius" or rxn_type == "three-body-Arrhenius":
            custom_reactions[ireact] = ct.Reaction(
                reactions[ireact].reactants,
                reactions[ireact].products,
                ct.ArrheniusRate(coeff*reactions[ireact].rate.input_data['rate-constant']['A'],
                                 1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                                 1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
                third_body=custom_reactions[ireact].third_body)

            print(reactions[ireact].rate.input_data['rate-constant']['A'])
            print(custom_reactions[ireact].rate.input_data['rate-constant']['A'])

        if rxn_type == "falloff-Troe":

            low = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['low-P-rate-constant']['A'],
                                   1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['b'],
                                   1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['Ea'])

            high = ct.Arrhenius(coeff*reactions[ireact].rate.input_data['high-P-rate-constant']['A'],
                                    1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['b'],
                                    1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['Ea'])
            falloff_coeffs = np.array([
                    reactions[ireact].rate.input_data["Troe"]["A"],
                    reactions[ireact].rate.input_data["Troe"]["T3"],
                    reactions[ireact].rate.input_data["Troe"]["T1"],
                    reactions[ireact].rate.input_data["Troe"]["T2"]])
            custom_reactions[ireact] = ct.Reaction(
                reactants=reactions[ireact].reactants,
                products=reactions[ireact].products,
                rate=ct.TroeRate(low=low, high=high, falloff_coeffs=falloff_coeffs)
                )

            print(reactions[ireact].rate.input_data['low-P-rate-constant']['A'])
            print(custom_reactions[ireact].rate.input_data['low-P-rate-constant']['A'])

        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                           species=species, reactions=custom_reactions)

        for maxT in temp_array:

            #print(maxT)

            # mass ratios
            air = "O2:1.0,N2:3.76"
            fuel = "C2H4:1"
            gas.set_equivalence_ratio(phi=1.0, fuel=fuel, oxidizer=air)

            gas.TP = 900.0, 101325.0

            # Define a reactor
            #r = ct.IdealGasConstPressureReactor(gas, name="Batch Reactor")
            r = ct.IdealGasReactor(gas, name="Batch Reactor")
            net = ct.ReactorNet([r])
            T = r.T
            while T < maxT:
                net.step()
                T = r.T

            mass_fractions = r.Y
            temperature = T

            y = mass_fractions
            y = np.where(np.less(y,1e-16), 1e-16, y)
            y = np.nan_to_num(y, copy=True, nan=1e-16, posinf=None, neginf=None)
            guess_temp = temperature

            jac_ct = gas.net_production_rates_ddX
            jac_ct = jac_ct/gas.molecular_weights*gas.mean_molecular_weight

            jac_matrix = jac_ct

            w = np.real(la.eigvals(jac_matrix))

            timescale = np.nan_to_num(1.0/w, copy=True, nan=1.0, posinf=None, neginf=None)

            temp_history[ii, 0] = temperature
            temp_history[ii,1:] = timescale
            ii = ii + 1

            np.savetxt(f'pert_timescale_reaction{ireact:03d}.dat', temp_history[:ii])

