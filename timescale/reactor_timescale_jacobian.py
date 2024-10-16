import cantera as ct
import pyrometheus as pyro
from pyro_jax import *
from functools import partial
from jax import jit
from jax import jacfwd
import numpy as np

np.set_printoptions(precision=8, linewidth=256)

def newton(fun, jac, y, *args, num_it=10, tol=1e-6):
    """Newton method"""
    for it in range(num_it):    
        dy = jnp.linalg.solve(jac(y, *args), fun(y, *args))
        y -= dy
        if jnp.linalg.norm(dy) < tol:
            return y


# Cantera solution object
sol = ct.Solution('mechanism.yaml')

# Pyrometheus-generated code
pyro_class = pyro.codegen.python.get_thermochem_class(sol)
pyro_gas = make_jax_pyro_class(pyro_class, jnp)

import numpy.linalg as la

temp_array = np.linspace(1000,2800,101)
temp_history = np.zeros((temp_array.shape[0],sol.n_species+1))
mass_history = np.zeros((temp_array.shape[0],sol.n_species+1))
mole_history = np.zeros((temp_array.shape[0],sol.n_species+1))

# Thermodynamic conditions
pressure = pyro_gas.one_atm

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

    # internal energy
    energy = pyro_gas.get_mixture_internal_energy_mass(temperature, mass_fractions)

    # density
    density = pyro_gas.get_density(pressure, temperature, mass_fractions)

    # JIT-compiled functions
    @jit
    def get_temperature(mass_fractions):
        return pyro_gas.get_temperature_energy(mass_fractions, energy, guess_temp)

    @jit
    def get_pressure(mass_fractions, temperature):
        return pyro_gas.get_pressure(density, temperature, mass_fractions)

    @jit
    def chemical_source_term(mass_fractions):
        temperature = pyro_gas.get_temperature_energy(mass_fractions, energy,
                                                   guess_temp)
        return (
            pyro_gas.get_net_production_rates(density, temperature, mass_fractions)
        )
#        return (
#            pyro_gas.wts * pyro_gas.get_net_production_rates(density,
#                temperature, mass_fractions) / density
#        )

#    chemical_jacobian = jit(jacfwd(chemical_source_term))

    y = mass_fractions
    y = np.where(np.less(y,1e-16), 1e-16, y)
    y = np.nan_to_num(y, copy=True, nan=1e-16, posinf=None, neginf=None)
    guess_temp = temperature

#    jac_pm = chemical_jacobian(y)
    jac_ct = sol.net_production_rates_ddX
    jac_ct = jac_ct/sol.molecular_weights*sol.mean_molecular_weight

    jac_matrix = jac_ct

    w = np.real(la.eigvals(jac_matrix))

    timescale = np.nan_to_num(1.0/w, copy=True, nan=1.0, posinf=None, neginf=None)

    temp_history[ii, 0] = temperature
    temp_history[ii,1:] = timescale
    mass_history[ii, 0] = temperature
    mass_history[ii,1:] = y
    mole_history[ii, 0] = temperature
    mole_history[ii,1:] = sol.X
    ii = ii + 1

    np.savetxt('timescale.dat', temp_history[:ii])
    np.savetxt('species_mass.dat', mass_history[:ii])
    np.savetxt('species_mole.dat', mole_history[:ii])

