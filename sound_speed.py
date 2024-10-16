"""
Compute the "equilibrium" and "frozen" sound speeds for a gas

Requires: cantera >= 3.0.0
Keywords: thermodynamics, equilibrium
"""

import cantera as ct
import math
import numpy as np
import matplotlib.pyplot as plt


def equilSoundSpeeds(gas, rtol=1.0e-9, max_iter=20000):
    """
    Returns a tuple containing the equilibrium and frozen sound speeds for a
    gas with an equilibrium composition.  The gas is first set to an
    equilibrium state at the temperature and pressure of the gas, since
    otherwise the equilibrium sound speed is not defined.
    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate('TP', rtol=rtol, max_iter=max_iter)

    # save properties
    s0 = gas.s
    p0 = gas.P
    r0 = gas.density

    # perturb the pressure
    p1 = p0*1.0001

    # set the gas to a state with the same entropy and composition but
    # the perturbed pressure
    gas.SP = s0, p1

    # frozen sound speed
    afrozen = math.sqrt((p1 - p0)/(gas.density - r0))

    # now equilibrate the gas holding S and P constant
    gas.equilibrate('SP', rtol=rtol, max_iter=max_iter)

    # equilibrium sound speed
    aequil = math.sqrt((p1 - p0)/(gas.density - r0))

    # check against the built-in sound speed function
    #afrozen2 = gas.sound_speed
    afrozen2 = math.sqrt(gas.cp/gas.cv*ct.gas_constant/gas.mean_molecular_weight*gas.T)

    return aequil, afrozen, afrozen2


# test program
if __name__ == "__main__":
    gas = ct.Solution('gri30.yaml')
    gas.X = 'CH4:1.00, O2:2.0, N2:7.52'
    T_range = np.linspace(300, 5001, 100)
    aux = np.zeros((4, T_range.shape[0]))
    ii = 0
    for T in T_range:
        gas.TP = T, ct.one_atm
        #print(T, equilSoundSpeeds(gas))
        aequil, afrozen, afrozen2 = equilSoundSpeeds(gas)
        aux[0,ii] = T
        aux[1,ii] = aequil
        aux[2,ii] = afrozen
        aux[3,ii] = afrozen2
        ii = ii + 1

    plt.plot(aux[0,:], aux[1,:])
    plt.plot(aux[0,:], aux[2,:])
    plt.plot(aux[0,:], aux[3,:])
    plt.show()
