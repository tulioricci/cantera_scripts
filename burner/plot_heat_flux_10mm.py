import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import *

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

import sys

#mechanism = "uiuc20sp"
mechanism = "wang99_51sp"
volume_per_min = 25

phi_array = np.array([0.55, 0.7, 0.85, 1.0, 1.15, 1.3, 2.0])
path_list = ["csv_10mm_Y_unity-Lewis-number",
             "csv_10mm_Y_mixture-averaged",
             "csv_10mm_Y_multicomponent",
             "csv_10mm_Y_multicomponent_soret"
            ]

color_list = {
              "csv_10mm_Y_unity-Lewis-number": "red",
              "csv_10mm_Y_mixture-averaged": "green",
              "csv_10mm_Y_multicomponent": "blue",
              "csv_10mm_Y_multicomponent_soret": "black"
             }

label_list = {
              "csv_10mm_Y_unity-Lewis-number": "Le=1",
              "csv_10mm_Y_mixture-averaged": "MixAvg",
              "csv_10mm_Y_multicomponent": "Multicomp.",
              "csv_10mm_Y_multicomponent_soret": "Mc. + Soret"
             }

##########################################################

plt.rcParams["axes.labelsize"] = 18
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["figure.dpi"] = 200

plt.close('all')
fig = plt.figure(1,figsize=[6.4,5.2])
ax1 = fig.add_subplot(111)
ax2 = plt.twiny()
ax3 = plt.twinx()
ax1.set_position([0.15,0.12,0.8,0.84])
ax2.set_position([0.15,0.12,0.8,0.84])
ax3.set_position([0.15,0.12,0.8,0.84])
ax1.tick_params(axis='both', bottom=True, top=False, left=True, right=False)
ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
ax3.yaxis.set_minor_locator(AutoMinorLocator())

for path in path_list:

    plot = []

    for phi in phi_array:

        filename = (
            path + '/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        data = np.loadtxt(filename, skiprows=1, delimiter=",")
        print("Solution read from " + filename)

        nspecies = data.shape[1]-6
        idx = 4+nspecies
        grid = data[:,0]
        velocity = data[:,1]
        T = data[:,2]
        density = data[:,3]
        Y = data[:,4:idx]
        mu = data[:,idx]
        kappa = data[:,idx+1]

        dT = T[-2] - T[-1]
        dx = grid[-2] - grid[-1]
        flux_N = kappa[-1]*dT/dx/10000.0
        flux_0 = -kappa[0]*(T[1] - T[0])/(grid[1] - grid[0])/10000.0
            
        plot.append([phi, mu[-1], kappa[-1], dT/dx, flux_N, flux_0])

    plot = np.array(plot)

    ax1.plot(plot[:,0], -plot[:,4], '-', lw=3.0, color=color_list[path],
             marker='s', markersize=6, label=label_list[path])

ax1.set_xlim(0.48, 1.32)
ax2.set_xlim(0.48, 1.32)
ax1.set_ylim(4.2, 6.0)
ax3.set_ylim(4.2, 6.0)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
ax1.set_ylabel(r'$\bf{Heat \;  flux \; (W/cm^2)}$')
ax1.legend(loc='upper left', bbox_to_anchor=(-0.01,1.01), fontsize=13)
legend = ax1.get_legend()
#kk = 0
#for text in legend.get_texts():
#    text.set_color(ax1.get_lines()[kk].get_color())
#    kk = kk + 1
#    if kk == 7:
#        break
fig_name = 'heat_flux_' + mechanism + '_m' + str('%4.2f' % volume_per_min)
plt.savefig(fig_name + '.png')
plt.show()


#if plot_all:
#    plt.close('all')
#    fig = plt.figure(1,figsize=[6.4,5.2])
#    ax1 = fig.add_subplot(111)
#    ax2 = plt.twiny()
#    ax3 = plt.twinx()
#    ax1.set_position([0.15,0.12,0.8,0.84])
#    ax2.set_position([0.15,0.12,0.8,0.84])
#    ax3.set_position([0.15,0.12,0.8,0.84])
#    ax1.tick_params(axis='both', labelsize=16, bottom=True, top=False, left=True, right=False)
#    ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#    ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#    ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
#    ax2.xaxis.set_minor_locator(AutoMinorLocator())
#    ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
#    ax3.yaxis.set_minor_locator(AutoMinorLocator())

#    #ax1.plot(     cantera_llnl[:,0],      -cantera_llnl[:,5], '-.', lw=3.0, color='magenta', marker='s', markersize=6, label='Cantera, LLNL (1D)')
#    #ax1.plot(cantera_gri30_mod[:,0], -cantera_gri30_mod[:,5], '-.', lw=3.0, color='cyan', marker='s', markersize=6, label='Cantera, GRI mod (1D)')
#    ax1.plot(cantera_wang_51sp[:,0], -cantera_wang_51sp[:,5], '-.', lw=3.0, color='red', marker='s', markersize=6, label='Cantera, 52sp (1D)')
#    ax1.plot(cantera_wang_51sp_N[:,0], -cantera_wang_51sp_N[:,5], ':', lw=3.0, color='gray', marker='s', markersize=6, label='Cantera, 52sp (1D)')
#    ax1.plot(cantera_wang_22sp[:,0], -cantera_wang_22sp[:,5], '-.', lw=3.0, color='tomato', marker='s', markersize=6, label='Cantera, 22sp (1D)')
#    ax1.plot(cantera_uiuc_20sp[:,0], -cantera_uiuc_20sp[:,5], '-.', lw=3.0, color='forestgreen', marker='P', markersize=8, label='Cantera, 20sp (1D)')
#    #ax1.plot( cantera_uiuc_7sp[:,0],  -cantera_uiuc_7sp[:,5], '-.', lw=3.0, color='dodgerblue', marker='X', markersize=8, label='Cantera, 7sp (1D)')
#    ax1.set_xlim(0.48, 1.32)
#    ax2.set_xlim(0.48, 1.32)
#    ax1.set_ylim(1.0, 13.0)
#    ax3.set_ylim(1.0, 13.0)
#    ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
#    ax1.set_ylabel(r'$\bf{Heat \;  flux \; (W/cm^2)}$')
#    ax1.legend(loc='upper left', bbox_to_anchor=(-0.01,1.01), fontsize=13)
#    #ax3.legend(loc='upper left', Bbox_to_anchor=(-0.01,0.79), fontsize=13)
#    legend = ax1.get_legend()
#    kk = 0
#    for text in legend.get_texts():
#        text.set_color(ax1.get_lines()[kk].get_color())
#        kk = kk + 1
#        if kk == 7:
#            break
#    plt.savefig('burner_heat_flux_m' + str('%4.2f' % volume_per_min) + '_0.png')
#    plt.show()


#if plot_all:
#    plt.close('all')
#    fig = plt.figure(1,figsize=[6.4,5.2])
#    ax1 = fig.add_subplot(111)
#    ax2 = plt.twiny()
#    ax3 = plt.twinx()
#    ax1.set_position([0.15,0.12,0.8,0.84])
#    ax2.set_position([0.15,0.12,0.8,0.84])
#    ax3.set_position([0.15,0.12,0.8,0.84])
#    ax1.tick_params(axis='both', labelsize=16, bottom=True, top=False, left=True, right=False)
#    ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#    ax3.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False) 
#    ax1.xaxis.set_minor_locator(AutoMinorLocator()) 
#    ax2.xaxis.set_minor_locator(AutoMinorLocator())
#    ax1.yaxis.set_minor_locator(AutoMinorLocator()) 
#    ax3.yaxis.set_minor_locator(AutoMinorLocator())

#    #ax1.plot(     cantera_llnl[:,0],      -cantera_llnl[:,3], '-.', lw=3.0, color='magenta', marker='s', markersize=6, label='Cantera, LLNL (1D)')
#    #ax1.plot(cantera_gri30_mod[:,0], -cantera_gri30_mod[:,3], '-.', lw=3.0, color='cyan', marker='s', markersize=6, label='Cantera, GRI mod (1D)')
#    ax1.plot(cantera_wang_51sp[:,0], -cantera_wang_51sp[:,3], '-.', lw=3.0, color='red', marker='s', markersize=6, label='Cantera, 52sp (1D)')
#    ax1.plot(cantera_wang_51sp_N[:,0], -cantera_wang_51sp_N[:,3], ':', lw=3.0, color='gray', marker='s', markersize=6, label='Cantera, 52sp (1D)')
#    ax1.plot(cantera_wang_22sp[:,0], -cantera_wang_22sp[:,3], '-.', lw=3.0, color='tomato', marker='s', markersize=6, label='Cantera, 22sp (1D)')
#    ax1.plot(cantera_uiuc_20sp[:,0], -cantera_uiuc_20sp[:,3], '-.', lw=3.0, color='forestgreen', marker='P', markersize=8, label='Cantera, 20sp (1D)')
#    #ax1.plot( cantera_uiuc_7sp[:,0],  -cantera_uiuc_7sp[:,3], '-.', lw=3.0, color='dodgerblue', marker='X', markersize=8, label='Cantera, 7sp (1D)')
#    ax1.set_xlim(0.48, 1.32)
#    ax2.set_xlim(0.48, 1.32)
#    ax1.set_ylim(0.8e6, 1.42e6)
#    ax3.set_ylim(0.8e6, 1.42e6)
#    ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
#    ax1.set_ylabel(r'$\bf{Temperature \; gradient \; (K/m^2)}$')
#    ax1.legend(loc='upper left', bbox_to_anchor=(-0.01,1.01), fontsize=13)
#    #ax3.legend(loc='upper left', Bbox_to_anchor=(-0.01,0.79), fontsize=13)
#    legend = ax1.get_legend()
#    kk = 0
#    for text in legend.get_texts():
#        text.set_color(ax1.get_lines()[kk].get_color())
#        kk = kk + 1
#        if kk == 7:
#            break
#    plt.savefig('burner_grad_T_m' + str('%4.2f' % volume_per_min) + '.png')
#    plt.show()



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

for path in path_list:

    plot = []

    for phi in phi_array:

        filename = (
            path + '/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        data = np.loadtxt(filename, skiprows=1, delimiter=",")
        print("Solution read from " + filename)

        nspecies = data.shape[1]-6
        idx = 4+nspecies
        grid = data[:,0]
        velocity = data[:,1]
        T = data[:,2]
        density = data[:,3]
        Y = data[:,4:idx]
        mu = data[:,idx]
        kappa = data[:,idx+1]

        dT = T[-2] - T[-1]
        dx = grid[-2] - grid[-1]
        flux_N = kappa[-1]*dT/dx/10000.0
        flux_0 = -kappa[0]*(T[1] - T[0])/(grid[1] - grid[0])/10000.0
            
        plot.append([phi, mu[-1], kappa[-1], dT/dx, flux_N, flux_0])

    plot = np.array(plot)

    ax1.plot(plot[:,0], plot[:,2], '-', lw=3.0, color=color_list[path],
             marker='s', markersize=6, label=label_list[path])

ax1.set_xlim(0.48, 1.32)
ax2.set_xlim(0.48, 1.32)
ax1.set_ylim(0.024, 0.03)
ax3.set_ylim(0.024, 0.03)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
ax1.set_ylabel(r'$\bf{Thermal \; conductivity \; (W/m-K)}$')
ax1.legend(loc='upper left', fontsize=13)
legend = ax1.get_legend()
#kk = 0
#for text in legend.get_texts():
#    text.set_color(ax1.get_lines()[kk].get_color())
#    kk = kk + 1
fig_name = 'heat_flux_transp_kappa_' + mechanism + '_m' + str('%4.2f' % volume_per_min)
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

for path in path_list:

    plot = []

    for phi in phi_array:

        filename = (
            path + '/stagnation_flame_' + mechanism +
            '_m' + str('%4.2f' % volume_per_min) + 
            '_phi' + str('%4.2f' % phi) + '.csv')

        data = np.loadtxt(filename, skiprows=1, delimiter=",")
        print("Solution read from " + filename)

        nspecies = data.shape[1]-6
        idx = 4+nspecies
        grid = data[:,0]
        velocity = data[:,1]
        T = data[:,2]
        density = data[:,3]
        Y = data[:,4:idx]
        mu = data[:,idx]
        kappa = data[:,idx+1]

        dT = T[-2] - T[-1]
        dx = grid[-2] - grid[-1]
        flux_N = kappa[-1]*dT/dx/10000.0
        flux_0 = -kappa[0]*(T[1] - T[0])/(grid[1] - grid[0])/10000.0
            
        plot.append([phi, mu[-1], kappa[-1], dT/dx, flux_N, flux_0])

    plot = np.array(plot)

    ax1.plot(plot[:,0], plot[:,1], '-', lw=3.0, color=color_list[path],
             marker='s', markersize=6, label=label_list[path])

ax1.set_xlim(0.48, 1.32)
ax2.set_xlim(0.48, 1.32)
#if volume_per_min == 17:
    #ax1.set_ylim(0.024, 0.034)
    #ax3.set_ylim(0.024, 0.034)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \mathbf{\Phi}$')
ax1.set_ylabel(r'$\bf{Viscosity \; (Pa \, s)}$')
ax1.legend(loc='upper right', fontsize=13)
legend = ax1.get_legend()
kk = 0
for text in legend.get_texts():
    text.set_color(ax1.get_lines()[kk].get_color())
    kk = kk + 1
plt.savefig('heat_flux_transp_mu_' + mechanism + '_m' + str('%4.2f' % volume_per_min) + '.png')
#plt.show()
plt.close()


"""
x = np.linspace(0.55,2.05,101)
a = np.interp(x, grad_25slpm_wang_51sp[:,0], grad_25slpm_wang_51sp[:,2])
b = np.interp(x,       grad_25slpm_7sp[:,0],       grad_25slpm_7sp[:,2])
c = np.interp(x,       grad_25slpm_20sp[:,0],     grad_25slpm_20sp[:,2])
"""

"""
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
ax1.plot(x, (a-b)/a, '--', lw=3.0, color='limegreen', label='Cantera, 7sp')
ax1.plot(x, (a-c)/a, '--', lw=3.0, color='lightcoral', label='Cantera, 20sp')
ax1.set_xlim(0.48, 2.02)
ax2.set_xlim(0.48, 2.02)
if volume_per_min == 17:
    ax1.set_ylim(-0.005, 0.255)
    ax3.set_ylim(-0.005, 0.255)
ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \mathbf{\Phi}$')
#ax3.set_ylabel(r'$\bf{Viscosity \; (xxx)}$')
ax1.set_ylabel(r'$\bf{Thermal \; conductivity \; error}$')
ax1.legend(loc='upper left', fontsize=13)
legend = ax1.get_legend()
kk = 0
for text in legend.get_texts():
    text.set_color(ax1.get_lines()[kk].get_color())
    kk = kk + 1
plt.savefig('burner_heat_flux_transport_error' + '_m' + str('%4.2f' % volume_per_min) + '.png')
#plt.show()
plt.close()
"""
