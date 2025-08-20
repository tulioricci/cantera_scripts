"""
Minimal reaction path diagram example
Author: Tilman Bremer
Date: June 21st 2017
 
Written for Python 3.6 with cantera 2.3.0 but might as well work with other versions.
"""
import os
import cantera as ct
import graphviz
import sys
 
# Define a gas mixture at a high temperature that will undergo a reaction:
#gas = ct.Solution('uiuc_based_on_usc_v5a.yaml')
gas = ct.Solution('../wang99_reduced.yaml')

air = "O2:0.21,N2:0.79"
fuel = "C2H4:1"
gas.set_equivalence_ratio(phi=1, fuel=fuel, oxidizer=air)

gas.TP = 800, 101325

for maxT in [800, 1200, 1600, 2000, 2400]:
 
    # Define a reactor, let it react until the temperature reaches 1800 K:
    r = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([r])
    T = r.T
    while T < maxT:
        net.step()
        T = r.T


    


    # Define the element to follow in the reaction path diagram:
    element = 'H'

    print( r.Y[gas.species_index('C2H4')] )
     
    # Initiate the reaction path diagram:
    diagram = ct.ReactionPathDiagram(gas, element)
     
    # Options for cantera:
    diagram.show_details = False
    diagram.font='CMU Serif Roman'
    diagram.threshold=0.05
    diagram.dot_options='node[fontsize=20,shape="box"]'
    diagram.title = 'Reaction path diagram following {0}'.format(element)
     
    # Define the filenames:
    dot_file = 'ReactionPathDiagram_' + element + '_' + str('%4i' % T) + '.dot'
    img_file = 'ReactionPathDiagram_' + element + '_' + str('%4i' % T) + '.png'
     
    # Write the dot-file first, then create the image from the dot-file with customizable
    # parameters:
    diagram.write_dot(dot_file)
    # The last command requires dot to be in your system path variables, or your system
    # will not undersand the command "dot".
    # The command -Tpng defines the filetype and needs to fit your filename defined above,
    # or else you will get errors opening the file later.
    # The command -Gdpi sets the resolution of the generated image in dpi.
    os.system('dot {0} -Gdpi=300 -Tpng -o{1}'.format(dot_file, img_file))

sys.exit()
