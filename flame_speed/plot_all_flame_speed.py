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

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams['text.latex.preamble']=[r'\boldmath']

plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["legend.fontsize"] = 15
#plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.dpi"] = 120

fig = plt.figure(1,figsize=[9.6,7.2])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
ax1.set_position([0.10,0.52,0.4,0.4])
ax2.set_position([0.52,0.52,0.4,0.4])
ax3.set_position([0.10,0.10,0.4,0.4])
ax4.set_position([0.52,0.10,0.4,0.4])

ax1.tick_params(axis='both', bottom=True, top=True, left=True, right=True,
                labelbottom=False, labeltop=False, labelleft=True, labelright=False)
ax2.tick_params(axis='both', bottom=True, top=True, left=True, right=True,
                labelbottom=False, labeltop=False, labelleft=False, labelright=True) 
ax3.tick_params(axis='both', bottom=True, top=True, left=True, right=True,
                labelbottom=True, labeltop=False, labelleft=True, labelright=False) 
ax4.tick_params(axis='both', bottom=True, top=True, left=True, right=True,
                labelbottom=True, labeltop=False, labelleft=False, labelright=True) 

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
         "Davis2005": 'blue'
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
              "blanquart": 'Blanquart',
              "Davis2005": 'Davis (H2)'
              }


mech_list = ['gri30',
             'gri30_33sp',
             #'../ucsd',
             'wang99_75sp',
             'blanquart',
             #'Davis2005'
             #'../wang99_reduced',
             #'../wang99_noFallOff',
             '../uiuc_20sp',
             #'../uiuc_20sp_v0',
             #'../uiuc_18sp',
#             '../uiuc_13sp_C',
#             '../uiuc_13sp_b',
#             '../uiuc_12sp',
#             '../uiuc_11sp_a',
#             '../uiuc_10sp_a',
             #'../uiuc_7sp',
#             '../uiuc_7sp_b'
            ]


ax1.text(0.6, 2.00, "H2")
ax3.text(0.6, 0.20, "CH4")
ax2.text(0.6, 1.45, "C2H2")
ax4.text(0.6, 0.80, "C2H4")

jj = 0
pressure = 1.0
H2 = 0.0
transport = "mix"
axis_list = [ax1, ax2, ax3, ax4]
for fuel in ["hydrogen", "acetylene", "methane", "ethylene"]:
    kk = 0

    axis = axis_list[jj]
    print(axis)

    for mechanism in mech_list:

        rxnmech = mechanism.replace("../","")
        
        ii = 0
        phi_array = np.hstack(( np.linspace(0.55,1.35,17), [1.5, 1.75, 1.9]))
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
        if jj == 0:
            axis.plot(phi_array, speed, color=color[rxnmech], linestyle=linestyle[kk], label=mech_label[rxnmech])
        else:
            axis.plot(phi_array, speed, color=color[rxnmech], linestyle=linestyle[kk])

        kk = kk + 1

    jj = jj + 1

#ax3.set_xlabel(r'Eq. ratio $\phi$')
#ax3.set_ylabel(r'Flame speed $S_L$ (m/s)')

#ax1.legend(loc='center', bbox_to_anchor=(1.0,1.0))
ax1.legend(loc='lower right')

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

# egolfopoulos
exp = np.asarray([
0.5017899761336515, 10.413223140495868,
0.6002386634844868, 21.5702479338843,
0.7022673031026253, 35.20661157024794,
0.8025059665871122, 46.11570247933884,
0.899164677804296, 58.01652892561984,
1.0011933174224343, 64.46280991735537,
1.0495226730310263, 66.19834710743802,
1.1032219570405728, 68.67768595041322,
1.2016706443914082, 70.41322314049587,
1.301909307875895, 66.19834710743802,
1.402147971360382, 60.0,
1.5041766109785204, 53.80165289256198,
1.5525059665871122, 47.35537190082645,
1.701073985680191, 30.24793388429752,
1.8997613365155133, 22.31404958677686,
0.5304295942720764, 14.62809917355372,
]).reshape(16,2)
ax4.scatter(exp[:,0],exp[:,1]/100, s=48, color="black")

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
ax4.scatter(exp[:,0],exp[:,1], s=48, marker="d", color="red")


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
ax4.scatter(exp[:,0],exp[:,1]/100, s=48, marker="x", color="blue")


exp = np.asarray([
0.5497847919655667, 6.7039800995024805,
0.5979913916786227, 12.101990049751237,
0.6995695839311336, 20.199004975124375,
0.7494978479196557, 24.2910447761194,
0.8011477761836443, 28.121890547263682,
0.9010043041606888, 35.69651741293532,
0.999139167862267, 40.13681592039801,
1.1041606886657103, 41.61691542288557,
1.2022955523672885, 36.39303482587064,
1.3004304160688667, 26.99004975124378,
1.4002869440459111, 17.41293532338308,
1.4984218077474893, 9.838308457711442
]).reshape(12,2)

ax3.scatter(exp[:,0],exp[:,1]/100, s=48, color="black")


exp = np.asarray([
0.5914893617021277, 0.6997319034852547,
0.795744680851064, 1.059249329758713,
0.9945288753799394, 1.3584450402144772,
1.2024316109422495, 1.5056300268096516,
1.399392097264438, 1.5587131367292226,
1.5945288753799396, 1.4332439678284183,
1.7951367781155017, 1.1557640750670242,
]).reshape(7,2)

ax2.scatter(exp[:,0],exp[:,1], s=48, color="black")

exp = np.asarray([
0.3511111111111111, 0.1360381861575179,
0.4, 0.25059665871121717,
0.4488888888888889, 0.3794749403341289,
0.49777777777777776, 0.5513126491646778,
0.5466666666666666, 0.7589498806682577,
0.6311111111111111, 1.0310262529832936,
0.6977777777777778, 1.224343675417661,
0.8088888888888889, 1.4677804295942722,
0.8266666666666667, 1.632458233890215,
0.9955555555555555, 2.0763723150357998,
1.0977777777777777, 2.2768496420047732,
0.25577464788732396, 0.07734806629834254,
0.2743661971830986, 0.11602209944751381,
0.29295774647887324, 0.13052486187845302,
0.3233802816901409, 0.19820441988950274,
0.3656338028169014, 0.2803867403314917,
0.4028169014084507, 0.343232044198895,
0.4416901408450704, 0.4060773480662983,
0.47887323943661975, 0.5269337016574586,
0.5194366197183099, 0.6477900552486188,
0.5549295774647888, 0.7783149171270718,
0.5971830985915494, 0.9136740331491713,
0.6343661971830986, 1.0393646408839778,
0.6969014084507043, 1.242403314917127,
0.7611267605633802, 1.4116022099447514,
0.8118309859154931, 1.4937845303867403,
0.8304225352112677, 1.6484806629834254,
0.9064788732394367, 1.8611878453038673,
1.0180281690140847, 2.0980662983425415,
1.1025352112676057, 2.286602209944751,
1.2022535211267606, 2.4316298342541436,
1.4980281690140844, 2.852209944751381,
]).reshape(32,2)

ax1.scatter(exp[:,0],exp[:,1], s=48, color="black")

#plt.legend(loc='lower right', fontsize=14)
ax1.set_xlim(0.5,1.5)
ax3.set_xlim(0.5,1.5)
ax2.set_xlim(0.5,1.5)
ax4.set_xlim(0.5,1.5)
ax1.set_ylim(0.5,3.0)
ax3.set_ylim(0.0,0.5)
ax2.set_ylim(0.5,1.75)
ax4.set_ylim(0.0,1.0)
plt.savefig('flame_speed_p' + str('%4.2f' % pressure) + '.png', dpi=200)
plt.show()


