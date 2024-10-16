import cantera as ct
import numpy as np
from functools import partial
import numpy.linalg as la

# get reaction coeffs
# https://groups.google.com/g/cantera-users/c/ID-oraSnUrU


np.set_printoptions(precision=8, linewidth=256)

# Cantera solution object
gas0 = ct.Solution('mechanism.yaml')
species = gas0.species()
reactions = gas0.reactions()

for ireact in range(0,len(reactions)):
    
    ii = 0
    temp_array = np.linspace(1000,2800,11)
    temp_history = np.zeros((temp_array.shape[0],gas0.n_species+1))
    pressure = 101325.0

    print(ireact, reactions[ireact].equation, reactions[ireact].reaction_type)

    custom_reactions = [r for r in reactions]
    if custom_reactions[ireact].reaction_type.find("falloff") < 0:
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(1.1*reactions[ireact].rate.input_data['rate-constant']['A'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
            third_body=custom_reactions[ireact].third_body)

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

