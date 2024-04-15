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

eps = 0.9

perturbed_species = "CH2*"


if perturbed_species is None:
    reaction_list = reactions
    idx_reacts = 0 #FIXME
else:
    ireact = -1
    idx_reacts = []
    reaction_list = []
    for R in reactions:
        ireact += 1
        if any(perturbed_species in reactant for reactant in R.reactants) or any(perturbed_species in product for product in R.products):
            reaction_list.append(R)
            idx_reacts.append(ireact)
        else:
            continue
        print(R)


for idx in range(0,len(reaction_list)):
    ireact = idx_reacts[idx]
    print(idx, ireact, reactions[ireact], reactions[ireact].reaction_type)

    ii = 0
    temp_array = np.linspace(1000,2800,11)
    temp_history = np.zeros((temp_array.shape[0],gas0.n_species+1))
    pressure = 101325.0

    custom_reactions = [r for r in reactions]

    #print(ireact, reactions[ireact].equation, reactions[ireact].reaction_type)

    if reactions[ireact].reaction_type == "Arrhenius":
        custom_reactions[ireact] = ct.Reaction(
            reactions[ireact].reactants,
            reactions[ireact].products,
            ct.ArrheniusRate(eps*reactions[ireact].rate.input_data['rate-constant']['A'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                             1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']))

    if reactions[ireact].reaction_type == "three-body-Arrhenius":
        print(reactions[ireact].third_body.efficiencies)
        effs = reactions[ireact].third_body.efficiencies
        if len(effs) == 1:
            eff = ct.ThirdBody(collider=list(effs.keys())[0])            
            custom_reactions[ireact] = ct.Reaction(
                    equation=reactions[ireact].equation,
                    rate=ct.ArrheniusRate(eps*reactions[ireact].rate.input_data['rate-constant']['A'],
                                          1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                                          1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
                    )
        else:
            eff = ct.ThirdBody()
            eff.efficiencies = reactions[ireact].third_body.efficiencies
            custom_reactions[ireact] = ct.Reaction(
                equation=reactions[ireact].equation,
                rate=ct.ArrheniusRate(eps*reactions[ireact].rate.input_data['rate-constant']['A'],
                                      1.0*reactions[ireact].rate.input_data['rate-constant']['b'],
                                      1.0*reactions[ireact].rate.input_data['rate-constant']['Ea']),
                third_body=eff
                )


    if reactions[ireact].reaction_type == "falloff-Troe":
        low = ct.Arrhenius(eps*reactions[ireact].rate.input_data['low-P-rate-constant']['A'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['b'],
                           1.0*reactions[ireact].rate.input_data['low-P-rate-constant']['Ea'])
        high = ct.Arrhenius(eps*reactions[ireact].rate.input_data['high-P-rate-constant']['A'],
                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['b'],
                            1.0*reactions[ireact].rate.input_data['high-P-rate-constant']['Ea'])
        Troe = [reactions[ireact].rate.input_data['Troe']['A'],
                reactions[ireact].rate.input_data['Troe']['T3'],
                reactions[ireact].rate.input_data['Troe']['T1'],
                reactions[ireact].rate.input_data['Troe']['T2']]
        eff = ct.ThirdBody()
        eff.efficiencies = reactions[ireact].third_body.efficiencies

        custom_reactions[ireact] = ct.FalloffReaction(
            equation=reactions[ireact].equation,
            rate=ct.TroeRate(low=low, high=high, falloff_coeffs=Troe),
            third_body=eff)

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

