import cantera as ct
import numpy as np

np.set_printoptions(precision=8, linewidth=256)

# Cantera solution object
sol = ct.Solution('mechanism.yaml')

import numpy.linalg as la

temp_array = np.linspace(1000,2800,101)
temp_history = np.zeros((temp_array.shape[0],sol.n_species+1))
mass_history = np.zeros((temp_array.shape[0],sol.n_species+1))
mole_history = np.zeros((temp_array.shape[0],sol.n_species+1))

ii = 0
for maxT in temp_array:

    print(maxT)

    # mass ratios
    air = "O2:1.0,N2:3.76"
    fuel = "C2H4:1"
    sol.set_equivalence_ratio(phi=1.0, fuel=fuel, oxidizer=air)

    sol.TP = 900.0, 101325.0

    # Define a reactor
    #r = ct.IdealGasConstPressureReactor(sol, name="Batch Reactor")
    r = ct.IdealGasReactor(sol, name="Batch Reactor")
    net = ct.ReactorNet([r])
    T = r.T
    while T < maxT:
        net.step()
        T = r.T

    mass_fractions = r.Y
    temperature = T

#    jac_pm = chemical_jacobian(y)
    jac_ct = sol.net_production_rates_ddX
    jac_ct = jac_ct/sol.molecular_weights*sol.mean_molecular_weight

    jac_matrix = jac_ct

    w = np.real(la.eigvals(jac_matrix))

    timescale = np.nan_to_num(1.0/w, copy=True, nan=1.0, posinf=None, neginf=None)

    temp_history[ii, 0] = temperature
    temp_history[ii,1:] = timescale
    mass_history[ii, 0] = temperature
    mass_history[ii,1:] = sol.Y
    #mole_history[ii, 0] = temperature
    #mole_history[ii,1:] = sol.X
    ii = ii + 1

    np.savetxt('timescale.dat', temp_history[:ii])
    np.savetxt('species_mass.dat', mass_history[:ii])
    np.savetxt('species_mole.dat', mole_history[:ii])

