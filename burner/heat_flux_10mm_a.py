import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import *
import cantera as ct
import sys
from statistics import median

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"


splm250_exp_055 = np.array([4.03, 4.27, 4.38, 4.61, 4.38, 4.39, 4.48, 4.42, 4.35, 4.86])
splm250_exp_070 = np.array([5.01, 4.77, 4.60, 4.91, 4.54, 4.84, 4.80, 4.77, 4.56, 5.08])
splm250_exp_085 = np.array([5.38, 5.36, 5.39, 5.52, 5.41, 5.35, 5.26, 5.50, 5.78, 5.34])
splm250_exp_100 = np.array([6.97, 5.90, 5.96, 6.63, 6.21, 6.62, 6.40, 6.60, 6.46, 5.95])
splm250_exp_115 = np.array([5.77, 6.15, 5.60, 5.93, 5.84, 5.45, 6.47, 6.42, 5.70, 5.93])
splm250_exp_130 = np.array([6.16, 5.95, 6.11, 5.80, 5.75, 6.07, 6.27, 6.12, 5.46, 6.02])

#splm250_7sp = np.array([0.55, 2.74, 0.70, 3.04, 0.85, 3.35, 1.0, 3.81, 1.15, 3.52, 1.30, 3.46]).reshape(6,2)
splm250_20sp = np.array([0.55, 4.80, 0.70, 5.23, 0.85, 5.75, 1.0, 6.39, 1.15, 5.90, 1.30, 6.0]).reshape(6,2)

species_output = "Y"

volume_per_min = 25
mdot_rate = 0.1693

width = 0.010  # m

use_radiation = False

#transp_model = 'unity-Lewis-number'
#enable_soret = False

transp_model = 'mixture-averaged'
enable_soret = False

#transp_model = 'multicomponent'
##enable_soret = False
#enable_soret = True

####################################

path = "csv_10mm_" + species_output + "_" + transp_model
if enable_soret:
    path = "csv_10mm_" + species_output + "_" + transp_model + "_soret"

if use_radiation:
    path = path + "_radiation"

os.system("mkdir -p " + path)

#phi_array = np.hstack((np.arange(0.55,1.351,0.05),
#                       np.asarray([1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1])))
phi_array = np.array([0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 2.0])

ratio = 2
slope = 0.1
curve = 0.2
loglevel = 0

r_int = 2.38*25.4/2000 #radius, actually
A_int = np.pi*r_int**2
lmin_to_m3s = 1.66667e-5

print("")
cantera_uiuc_7sp = []

gas = ct.Solution("uiuc_7sp.yaml")
gas.transport_model = transp_model

burner_temp = 300.0
surf_temp = 300.0

for phi in phi_array:

    filename = (
        path + '/stagnation_flame_uiuc_7sp' +
        '_m' + str('%4.2f' % volume_per_min) + 
        '_phi' + str('%4.2f' % phi) + '.csv')

    if os.path.exists(filename):
        continue
    else:
        air = "O2:0.21,N2:0.79"
        fuel = "C2H4:1"
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

        gas.TP = 300.0, 101325.0
        temp, rho_int, x = gas.TDX

        mass_unb = volume_per_min/np.sum(x[:])
        r_int = 2.38*25.4/2000 #radius, actually
        A_int = np.pi*r_int**2
        lmin_to_m3s = 1.66667e-5
        u_int = mass_unb*lmin_to_m3s/A_int
        rhoU_int = rho_int*u_int
        mdot_int = rho_int*u_int*A_int

        try:
            sim = ct.ImpingingJet(gas=gas, width=width)
            sim.inlet.mdot = mdot_rate*1.5
            sim.inlet.T = burner_temp
            sim.surface.T = surf_temp
            sim.soret_enabled = enable_soret
            sim.set_initial_guess(products='equil')
            sim.radiation_enabled = use_radiation

            sim.set_refine_criteria(ratio=12, slope=1.0, curve=1.0, prune=0.2)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.set_refine_criteria(ratio=6, slope=1.0, curve=1.0, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.set_refine_criteria(ratio=4, slope=0.9, curve=0.9, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.set_refine_criteria(ratio=2, slope=0.5, curve=0.5, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.set_refine_criteria(ratio=2, slope=0.4, curve=0.4, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.refine()

            sim.inlet.mdot = mdot_rate*1.4
            sim.set_refine_criteria(ratio=2, slope=0.40, curve=0.38, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.inlet.mdot = mdot_rate*1.3
            sim.set_refine_criteria(ratio=2, slope=0.38, curve=0.38, prune=0.1)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.refine()

            sim.inlet.mdot = mdot_rate*1.2
            sim.set_refine_criteria(ratio=ratio, slope=0.3, curve=0.3, prune=0.075)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.inlet.mdot = mdot_rate*1.15
            sim.set_refine_criteria(ratio=ratio, slope=0.25, curve=0.25, prune=0.075)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            sim.inlet.mdot = mdot_rate*1.1
            sim.set_refine_criteria(ratio=ratio, slope=0.2, curve=curve, prune=0.05)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

#            sim.refine()

#            sim.inlet.mdot = mdot_rate*1.05
#            sim.set_refine_criteria(ratio=ratio, slope=0.15, curve=curve, prune=0.025)
#            sim.solve(loglevel, refine_grid=True, auto=True)
#            print(sim.T.shape)

            sim.inlet.mdot = mdot_rate*1.0
            sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=0.0)
            sim.solve(loglevel, refine_grid=True, auto=True)
            print(sim.T.shape)

            nspecies = gas.n_species
            idx = 4+nspecies
            ngrid = sim.T.shape[0]
            data = np.zeros((ngrid, 4+nspecies+2))
            data[:,0] = sim.grid
            data[:,1] = sim.velocity
            data[:,2] = sim.T
            data[:,3] = sim.density
            data[:,4:idx] = (sim.Y).T
            data[:,idx] = sim.viscosity
            data[:,idx+1] = sim.thermal_conductivity

            header = "grid, u, T, D, Y_" + ', Y_'.join(str(x) for x in gas.species_names) + ", mu, kappa"
            np.savetxt(filename, data, fmt="%18.12e", delimiter=", ", header=header)
            print("Solution saved to " + filename)

            print('m =', u_int*A_int/lmin_to_m3s,', phi =',phi, ', Tb =', sim.T[0], sim.T[-1], sim.velocity[0])

        except:
            print("Meh")
            continue
            #sys.exit()

print("")
cantera_uiuc_20sp = []

gas = ct.Solution("uiuc_20sp.yaml")
gas.transport_model = transp_model
for phi in phi_array:

    filename = (
        path + '/stagnation_flame_uiuc_20sp' +
        '_m' + str('%4.2f' % volume_per_min) + 
        '_phi' + str('%4.2f' % phi) + '.csv')

    if os.path.exists(filename):
        continue

    else:
        air = "O2:0.21,N2:0.79"
        fuel = "C2H4:1"
        gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

        gas.TP = 300.0, 101325.0
        temp, rho_int, x = gas.TDX

        mass_unb = volume_per_min/np.sum(x[:])
        r_int = 2.38*25.4/2000 #radius, actually
        A_int = np.pi*r_int**2
        lmin_to_m3s = 1.66667e-5
        u_int = mass_unb*lmin_to_m3s/A_int
        rhoU_int = rho_int*u_int
        mdot_int = rho_int*u_int*A_int

        sim = ct.ImpingingJet(gas=gas, width=width)
        sim.inlet.mdot = mdot_rate
        sim.inlet.T = burner_temp
        sim.surface.T = surf_temp
        sim.soret_enabled = enable_soret
        sim.set_initial_guess(products='equil')
        sim.radiation_enabled = use_radiation

        sim.set_refine_criteria(ratio=4, slope=0.9, curve=0.9, prune=0.1)
        sim.solve(loglevel, refine_grid=True, auto=True)
        print(sim.T.shape)

        sim.set_refine_criteria(ratio=2, slope=0.5, curve=0.5, prune=0.1)
        sim.solve(loglevel, refine_grid=True, auto=True)
        print(sim.T.shape)

        sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=0.)
        sim.solve(loglevel, refine_grid=True, auto=False)
        print(sim.T.shape)

        nspecies = gas.n_species
        idx = 4+nspecies
        ngrid = sim.T.shape[0]
        data = np.zeros((ngrid, 4+nspecies+2))
        data[:,0] = sim.grid
        data[:,1] = sim.velocity
        data[:,2] = sim.T
        data[:,3] = sim.density
        data[:,4:idx] = (sim.Y).T
        data[:,idx] = sim.viscosity
        data[:,idx+1] = sim.thermal_conductivity

        header = "grid, u, T, D, Y_" + ', Y_'.join(str(x) for x in gas.species_names) + ", mu, kappa"
        np.savetxt(filename, data, fmt="%18.12e", delimiter=", ", header=header)
        print("Solution saved to " + filename)

        print('m =', u_int*A_int/lmin_to_m3s,', phi =',phi, ', Tb =', sim.T[0], sim.T[-1], sim.velocity[0])


my_dict = {}

for mechanism in ["uiuc_7sp", "uiuc_20sp", "wang99_51sp", "wang99_75sp", "gri30", "Blanquart2018", "AramcoMech3.0"]:
#for mechanism in ["uiuc_7sp", "uiuc_20sp", "wang99_51sp", "wang99_75sp", "gri30", "Blanquart2018"]:

    mech_data = []
    for phi in phi_array:

        filename = (
            path + '/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        print("")
        gas = ct.Solution(mechanism + ".yaml")
        gas.transport_model = transp_model

        print(os.path.exists(filename))
        if os.path.exists(filename):
            print("Solution read from " + filename)
            data = np.loadtxt(filename, skiprows=1, delimiter=",")

            nspecies = data.shape[1]-6
            idx = 4+nspecies
            grid = data[:,0]
            velocity = data[:,1]
            T = data[:,2]
            density = data[:,3]
            Y = data[:,4:idx]
            mu = data[:,idx]
            kappa = data[:,idx+1]

            gas.TPY = T[0], 101325.0, Y[0]
            _, rho_int, X = gas.TDX

            mass_unb = volume_per_min/np.sum(X[:])
            u_int = mass_unb*lmin_to_m3s/A_int
            rhoU_int = rho_int*u_int
            mdot_int = rho_int*u_int*A_int

            dT = T[-2] - T[-1]
            dx = grid[-2] - grid[-1]
            flux_N = kappa[-1]*dT/dx/10000.0
            flux_0 = -kappa[0]*(T[1] - T[0])/(grid[1] - grid[0])/10000.0
            
            mech_data.append([phi, mu[-1], kappa[-1], dT/dx, flux_N, flux_0])

            print('m =', u_int*A_int/lmin_to_m3s,', phi =', phi, ', Tb =', T[0], T[-1], velocity[0])

        else:
            air = "O2:0.21,N2:0.79"
            fuel = "C2H4:1"
            gas.set_equivalence_ratio(phi=phi, fuel=fuel, oxidizer=air)

            gas.TP = 300.0, 101325.0
            temp, rho_int, x = gas.TDX

            mass_unb = volume_per_min/np.sum(x[:])
            r_int = 2.38*25.4/2000 #radius, actually
            A_int = np.pi*r_int**2
            lmin_to_m3s = 1.66667e-5
            u_int = mass_unb*lmin_to_m3s/A_int
            rhoU_int = rho_int*u_int
            mdot_int = rho_int*u_int*A_int

            sim = ct.ImpingingJet(gas=gas, width=width)
            sim.inlet.mdot = mdot_rate
            sim.inlet.T = burner_temp
            sim.surface.T = surf_temp
            sim.soret_enabled = enable_soret
            sim.set_initial_guess(products='equil')
            sim.radiation_enabled = use_radiation

            sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=0.0)
            sim.solve(loglevel, refine_grid=True, auto=True)

            nspecies = gas.n_species
            idx = 4+nspecies
            ngrid = sim.T.shape[0]
            data = np.zeros((ngrid, 4+nspecies+2))
            data[:,0] = sim.grid
            data[:,1] = sim.velocity
            data[:,2] = sim.T
            data[:,3] = sim.density
            data[:,4:idx] = (sim.Y).T
            data[:,idx] = sim.viscosity
            data[:,idx+1] = sim.thermal_conductivity

            header = "grid, u, T, D, Y_" + ', Y_'.join(str(x) for x in gas.species_names) + ", mu, kappa"
            np.savetxt(filename, data, fmt="%18.12e", delimiter=", ", header=header)
            print("Solution saved to " + filename)

            print('m =', u_int*A_int/lmin_to_m3s,', phi =',phi, ', Tb =', sim.T[0], sim.T[-1], sim.velocity[0])
    
    my_dict.update({mechanism: mech_data})
    
"""########################################################################"""


cantera_uiuc_7sp = np.array(my_dict["uiuc_7sp"])
cantera_uiuc_20sp = np.array(my_dict["uiuc_20sp"])
cantera_wang_51sp = np.array(my_dict["wang99_51sp"])
cantera_wang_75sp = np.array(my_dict["wang99_75sp"])
if "gri30" in my_dict.keys():
    cantera_gri30 = np.array(my_dict["gri30"])
if "Blanquart2018" in my_dict.keys():
    cantera_blanq = np.array(my_dict["Blanquart2018"])
if "AramcoMech3.0" in my_dict.keys():
    cantera_aramco = np.array(my_dict["AramcoMech3.0"])

plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["figure.dpi"] = 200

plt.close('all')
fig = plt.figure(1,figsize=[6.4,5.2])
ax1 = fig.add_subplot(111)
ax2 = plt.twiny()
ax3 = plt.twinx()
ax1.set_position([0.15,0.12,0.8,0.84])
ax2.set_position([0.15,0.12,0.8,0.84])
ax3.set_position([0.15,0.12,0.8,0.84])
ax1.tick_params(axis='both', labelsize=16, bottom=True, top=False, left=True, right=False)
ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
ax3.yaxis.set_minor_locator(AutoMinorLocator())

ax1.plot(         splm250_20sp[:,0],           splm250_20sp[:,1],  '-', lw=3.0, color='darkgreen', marker='^', markersize=8, label='MIRGE-Com, 20sp (2D)')
#ax1.plot(          splm250_7sp[:,0],            splm250_7sp[:,1],  '-', lw=3.0, color='blue', marker='v', markersize=8, label='MIRGE-Com, 7sp (2D)')
#ax1.plot(        grad_25slpm_7sp[:,0],         -grad_25slpm_7sp[:,4], '-.', lw=3.0, color='dodgerblue', marker='X', markersize=8, label='Cantera, 7sp (1D)')
ax1.plot(  cantera_uiuc_7sp[:,0],  -cantera_uiuc_7sp[:,4], '-.', lw=3.0, color='dodgerblue', marker='X', markersize=8, label='Cantera, 7sp (1D)')
ax1.plot( cantera_uiuc_20sp[:,0], -cantera_uiuc_20sp[:,4], '-.', lw=3.0, color='forestgreen', marker='P', markersize=8, label='Cantera, 20sp (1D)')
ax1.plot( cantera_wang_51sp[:,0], -cantera_wang_51sp[:,4], '-.', lw=3.0, color='tomato', marker='d', markersize=6, label='Cantera, 51sp (1D)')
ax1.plot( cantera_wang_75sp[:,0], -cantera_wang_75sp[:,4], '-.', lw=3.0, color='red', marker='s', markersize=6, label='Cantera, 75sp (1D)')
if "gri30" in my_dict.keys():
    ax1.plot( cantera_gri30[:,0], -cantera_gri30[:,4], '-.', lw=3.0, color='purple', marker='s', markersize=6, label='Cantera, GRI 3.0 (1D)')
if "Blanquart2018" in my_dict.keys():
    ax1.plot( cantera_blanq[:,0], -cantera_blanq[:,4], '-.', lw=3.0, color='gold', marker='s', markersize=6, label='Cantera, Blanquart (1D)')
if "AramcoMech3.0" in my_dict.keys():
    ax1.plot( cantera_aramco[:,0], -cantera_aramco[:,4], '-.', lw=3.0, color='silver', marker='s', markersize=6, label='Cantera, Aramco (1D)')

ax1.plot(splm250_exp_055[:]*0 + 0.55, splm250_exp_055[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black', label="Experiment")
ax1.plot(splm250_exp_070[:]*0 + 0.70, splm250_exp_070[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black')
ax1.plot(splm250_exp_085[:]*0 + 0.85, splm250_exp_085[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black')
ax1.plot(splm250_exp_100[:]*0 + 1.00, splm250_exp_100[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black')
ax1.plot(splm250_exp_115[:]*0 + 1.15, splm250_exp_115[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black')
ax1.plot(splm250_exp_130[:]*0 + 1.30, splm250_exp_130[:], ' ', lw=0, marker='o', markerfacecolor='black', color='black')

ax1.plot(0.55, median(splm250_exp_055[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black', label="Exp. median")
ax1.plot(0.70, median(splm250_exp_070[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black')
ax1.plot(0.85, median(splm250_exp_085[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black')
ax1.plot(1.00, median(splm250_exp_100[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black')
ax1.plot(1.15, median(splm250_exp_115[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black')
ax1.plot(1.30, median(splm250_exp_130[:]), ' ', lw=0, marker='o', markerfacecolor='red', color='black')

ax1.set_xlim(0.48, 1.32)
ax2.set_xlim(0.48, 1.32)
ax1.set_ylim(4.0, 8.0)
ax3.set_ylim(4.0, 8.0)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
ax1.set_ylabel(r'$\bf{Heat \;  flux \; (W/cm^2)}$')
ax1.legend(loc='upper left', bbox_to_anchor=(-0.01,1.01), fontsize=13)
##ax3.legend(loc='upper left', Bbox_to_anchor=(-0.01,0.79), fontsize=13)
legend = ax1.get_legend()
kk = 0
for text in legend.get_texts():
    text.set_color(ax1.get_lines()[kk].get_color())
    kk = kk + 1
    if kk == 8:
        break
fig_name = 'burner_heat_flux' + '_m' + str('%4.2f' % volume_per_min) + "_" + transp_model
if enable_soret:
    fig_name = fig_name + "_soret"
if use_radiation:
    fig_name = fig_name + "_radiation"
plt.savefig(fig_name + '.png')
plt.show()
#plt.close()

plt.close('all')
fig = plt.figure(1,figsize=[6.4,5.2])
ax1 = fig.add_subplot(111)
ax2 = plt.twiny()
ax3 = plt.twinx()
ax1.set_position([0.15,0.12,0.8,0.84])
ax2.set_position([0.15,0.12,0.8,0.84])
ax3.set_position([0.15,0.12,0.8,0.84])
ax1.tick_params(axis='both', labelsize=16, bottom=True, top=False, left=True, right=False)
ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)
ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)  
ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
ax3.yaxis.set_minor_locator(AutoMinorLocator())
#ax1.plot([0.55, 0.70, 0.85, 1.0, 1.15, 1.30],
#         [1.025266617454133583, 1.025348659710407707, 1.025283900253967023, 1.02525578816328228, 1.02648720446326346, 1.02722594938282046], lw=3.0,
#         label="MIRGE-Com, 20sp (2D)", color="darkgreen", marker='^', markersize=8)
#ax1.plot([0.55, 0.70, 0.85, 1.0, 1.15, 1.30],
#         [1.025653007979627426, 1.025440074666339553, 1.02523589857371942, 1.02504272339411241, 1.025172615483559122, 1.025119101582028323], lw=3.0,
#         label="MIRGE-Com, 7sp (2D)", color="blue", marker='v', markersize=8)

ax1.plot( cantera_uiuc_7sp[:,0],  cantera_uiuc_7sp[:,2], '-.', lw=3.0, color='dodgerblue', marker='X', markersize=8, label='Cantera, 7sp (1D)')
ax1.plot(cantera_uiuc_20sp[:,0], cantera_uiuc_20sp[:,2], '-.', lw=3.0, color='forestgreen', marker='P', markersize=8, label='Cantera, 20sp (1D)')
ax1.plot(cantera_wang_51sp[:,0], cantera_wang_51sp[:,2], '-.', lw=3.0, color='tomato', marker='d', markersize=6, label='Cantera, 51sp (1D)')
ax1.plot(cantera_wang_75sp[:,0], cantera_wang_75sp[:,2], '-.', lw=3.0, color='red', marker='s', markersize=6, label='Cantera, 75sp (1D)')
if "gri30" in my_dict.keys():
    ax1.plot( cantera_gri30[:,0], cantera_gri30[:,2], '-.', lw=3.0, color='purple', marker='s', markersize=6, label='Cantera, GRI 3.0 (1D)')
if "Blanquart2018" in my_dict.keys():
    ax1.plot( cantera_blanq[:,0], cantera_blanq[:,2], '-.', lw=3.0, color='gold', marker='s', markersize=6, label='Cantera, Blanquart (1D)')
if "AramcoMech3.0" in my_dict.keys():
    ax1.plot( cantera_aramco[:,0], cantera_aramco[:,2], '-.', lw=3.0, color='silver', marker='s', markersize=6, label='Cantera, Aramco (1D)')

#ax1.plot( cantera_aramco[:,0], cantera_aramco[:,2], '-.', lw=3.0, color='silver', marker='s', markersize=6, label='Cantera, Aramco (1D)')
#ax1.plot([0.55, 0.70, 0.85, 1.0, 1.15, 1.30],
#         [0.025654016692188017, 0.025466471385977768, 0.025300670314592, 0.025112590915185175, 0.026825361474584836, 0.02784529872452974], lw=3.0,
#         color="darkgreen", marker='^')
#ax1.plot([0.55, 0.70, 0.85, 1.0, 1.15, 1.30],
#         [0.025653007979627426, 0.025440074666339553, 0.02523589857371942, 0.02504272339411241, 0.025172615483559122, 0.025119101582028323], lw=3.0,
#         color="blue", marker='v')
ax1.set_xlim(0.48, 2.02)
ax2.set_xlim(0.48, 2.02)
if volume_per_min == 17:
    ax1.set_ylim(0.024, 0.034)
    ax3.set_ylim(0.024, 0.034)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
#ax3.set_ylabel(r'$\bf{Viscosity \; (xxx)}$')
ax1.set_ylabel(r'$\bf{Thermal \; conductivity \; (W/m-K)}$')
ax1.legend(loc='upper left', fontsize=13)
legend = ax1.get_legend()
kk = 0
for text in legend.get_texts():
    text.set_color(ax1.get_lines()[kk].get_color())
    kk = kk + 1
plt.savefig('burner_heat_flux_kappa' + '_m' + str('%4.2f' % volume_per_min) + '.png')
plt.show()
#plt.close()

