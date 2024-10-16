"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0

Sensitivity analysis
https://cantera.org/examples/jupyter/flames/flame_speed_with_sensitivity_analysis.ipynb.html
"""

import cantera as ct
import matplotlib
import sys

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
width = 0.02  # m
loglevel = 0  # amount of diagnostic output (0 to 8)

import matplotlib.pyplot as plt
from matplotlib.ticker import *

#plt.rc('text', usetex=True)
#plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#plt.rcParams['text.latex.preamble']=[r'\boldmath']
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["figure.dpi"] = 200
plt.rcParams["legend.fontsize"] = 15

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

ax1.tick_params(axis='both', bottom=True, top=True, left=True, right=True,
                labelbottom=True, labeltop=False, labelleft=True,
                labelright=False)

import os
import numpy as np

linestyle = ['-', '--', '-', '-', '-']
color = {"uiuc_7sp": 'black',
         "uiuc_7sp_b": 'gray',
         #"uiuc_7sp_cp": 'gray',
         #"uiuc_7sp_powerLaw": 'silver',
         #"uiuc_10sp_a": 'blue',
         #"uiuc_11sp_a": 'teal',
         #"uiuc_12sp": 'mediumturquoise',
         #"uiuc_13sp_b": 'aquamarine',
         #"uiuc_13sp_C": 'gold',
         "uiuc_18sp": 'orange',
         #"uiuc_20sp_v0": 'tomato',
         "uiuc_20sp": 'green',
         "wang99_75sp": 'magenta',
         #"wang99_reduced": 'red',
         #"wang99_noFallOff": 'magenta',
         "gri30": 'tomato',
         "gri30_33sp": "firebrick",
         "ucsd": 'cyan',
         "blanquart": 'black',
         "Davis2005": 'blue',
         "aramco": 'gray'
         }
mech_label = {"uiuc_7sp": 'UIUC (7 sp)',
              "uiuc_7sp_b": 'UIUC (7 sp, b)',
              "uiuc_7sp_cp": r'UIUC$_{cp}$ (7 sp)',
              "uiuc_7sp_powerLaw": r'UIUC$_{PL}$ (7 sp)',
              "uiuc_20sp_v0": 'UIUC (20 sp, v0)',
              "uiuc_20sp": 'UIUC (20 sp)',
              "uiuc_18sp": 'UIUC (18 sp)',
              "uiuc_13sp_C": 'UIUC (13 sp, C)',
              "uiuc_13sp_b": 'Davis (13 sp, b)',
              "uiuc_12sp": 'Davis (12 sp)',              
              "uiuc_11sp_a": 'Davis (11 sp, a)',
              "uiuc_10sp_a": 'Davis (10 sp, a)',
              "wang99_75sp": 'Wang 99',
              "wang99_reduced": 'Wang 99 (50 sp)',
              "gri30": 'GRI 3.0',
              "gri30_33sp": "GRI 3.0 (33sp)",
              "ucsd": 'UCSD',
              "blanquart": 'Blanquart 2018',
              "Davis2005": 'Davis (H2)',
              "aramco": 'Aramco 3.0'
              }


mech_list = [#'gri30',
             #'gri30_33sp',
             'wang99_75sp',
             'blanquart',
             'aramco',
             '../uiuc_20sp',
            ]

jj = 0
pressure = 1.0
H2 = 0.0
transport = "mix"
fuel = "ethylene"
    
for mechanism in mech_list:

    rxnmech = mechanism.replace("../","")
    
    ii = 0
    phi_array = np.arange(0.55,1.80,0.05)
    speed = np.zeros((phi_array.shape[0]))

    for phi in phi_array:

        filename = (fuel + '/csv/adiabatic_flame_' +
            rxnmech +
            '_phi' + str('%4.2f' % phi) + 
            '_p' + str('%4.2f' % pressure) + 
            '_E:H' + str('%3.1f' % H2) + '_' + transport + '.csv')
        if os.path.exists(filename):
            #print(filename)
            with open(filename) as data:
                first_line = data.readline()
                header = first_line.split(',')
            data = np.genfromtxt(filename, skip_header=1, delimiter=',')
            speed[ii] = data[0,1]

        ii = ii + 1  #~~ phi

    speed = np.where(np.less(speed, 0.001), np.nan, speed)
    ax1.plot(phi_array, speed, color=color[rxnmech], linestyle="-", marker="", label=mech_label[rxnmech])

ax1.set_xlabel(r'$\bf{Equivalence \; ratio} \; \Phi$')
ax1.set_ylabel(r'$\bf{Flame \; speed} \; S_L \; (m/s)$')

## egolfopoulos
#exp = np.asarray([
#0.5017899761336515, 10.413223140495868,
#0.6002386634844868, 21.5702479338843,
#0.7022673031026253, 35.20661157024794,
#0.8025059665871122, 46.11570247933884,
#0.899164677804296, 58.01652892561984,
#1.0011933174224343, 64.46280991735537,
#1.0495226730310263, 66.19834710743802,
#1.1032219570405728, 68.67768595041322,
#1.2016706443914082, 70.41322314049587,
#1.301909307875895, 66.19834710743802,
#1.402147971360382, 60.0,
#1.5041766109785204, 53.80165289256198,
#1.5525059665871122, 47.35537190082645,
#1.701073985680191, 30.24793388429752,
#1.8997613365155133, 22.31404958677686,
#0.5304295942720764, 14.62809917355372,
#]).reshape(16,2)
#ax1.scatter(exp[:,0],exp[:,1]/100, s=48, color="red", label="Egolfopoulos et al 90")

exp = np.asarray([
0.4635451505016723, 8.645320197044335,
0.49698996655518396, 12.192118226600986,
0.5321070234113713, 16.625615763546797,
0.5521739130434783, 19.285714285714285,
0.5939799331103679, 24.16256157635468,
0.6508361204013378, 31.699507389162562,
0.6959866220735786, 36.354679802955665,
0.7963210702341137, 47.88177339901478,
0.8531772575250836, 52.98029556650246,
0.8933110367892977, 59.85221674876848,
0.9953177257525083, 66.2807881773399,
1.045484949832776, 68.27586206896552,
1.0989966555183948, 70.04926108374384,
1.1993311036789298, 72.04433497536947,
1.2963210702341137, 71.60098522167488,
1.3983277591973244, 65.83743842364532,
1.4986622073578597, 55.862068965517246,
1.5505016722408027, 49.21182266009853,
1.6943143812709032, 34.35960591133005,
1.8933110367892976, 24.605911330049263,
2.10066889632107, 18.177339901477833,
2.19933110367893, 15.517241379310345,
]).reshape(22,2)
ax1.scatter(exp[:,0],exp[:,1]/100, s=48, marker="o", color="black", facecolor="none", label="Egolfopoulos et al 90")


# hassan
exp = np.asarray([
0.8, 0.49,
0.9, 0.59,
1.0, 0.65,
1.1, 0.68,
1.2, 0.65,
1.4, 0.58,
1.6, 0.33,
]).reshape(7,2)
ax1.scatter(exp[:,0],exp[:,1], s=48, marker="d", color="red", facecolor="none", label="Hassan et al 98")


# jomaas
exp = np.asarray([
0.6, 23.180212014134277,
0.7019230769230769, 35.47703180212014,
0.798076923076923, 49.57597173144876,
0.8999999999999999, 56.572438162544174,
1.001923076923077, 63.03886925795053,
1.0999999999999999, 64.62897526501767,
1.2, 63.886925795053,
1.2999999999999998, 60.3886925795053,
1.4, 52.544169611307424,
]).reshape(9,2)
ax1.scatter(exp[:,0],exp[:,1]/100, s=48, marker="^", color="green", facecolor="none", label="Jomaas et al 05")

# kumar
exp = np.asarray([
0.5014207650273225, 17.38680465717982,
0.6010928961748634, 31.46183699870634,
0.7007650273224044, 43.26002587322122,
0.8021857923497269, 53.91979301423027,
0.9036065573770492, 62.82018111254852,
1.0032786885245901, 68.82276843467012,
1.1029508196721312, 71.72056921086676,
1.2008743169398908, 70.78913324708927,
1.3005464480874318, 63.95860284605434,
1.4037158469945354, 59.301423027166884,
]).reshape(10,2)
ax1.scatter(exp[:,0],exp[:,1]/100, s=48, marker="v", color="blue", facecolor="none", label="Kumar et al 08")



#ax1.legend(loc='center', bbox_to_anchor=(1.0,1.0))
ax1.legend(loc='lower center')

#leg = ax1.get_legend()
#kk = 0
#for mechanism in mech_list:
#    rxnmech = mechanism[3:]
#    leg.legendHandles[kk].set_color(color[rxnmech])
#    kk += 1

#legend = ax1.get_legend()
#kk = 0
#for text in legend.get_texts():
#    mechanism = mech_list[kk]
#    rxnmech = mechanism.replace("../","")
#    text.set_color(color[rxnmech])
#    kk = kk + 1


ax1.set_xlim(0.5,1.7)
ax1.set_ylim(0.0,0.8)
ax2.set_xlim(0.5,1.7)
ax2.set_ylim(0.0,0.8)
ax3.set_xlim(0.5,1.7)
ax3.set_ylim(0.0,0.8)
plt.savefig('ethylene_flame_speed_p' + str('%4.2f' % pressure) + '.png', dpi=200)
plt.show()


