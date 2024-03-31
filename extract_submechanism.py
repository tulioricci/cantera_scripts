import cantera as ct

input_file = 'uiuc_with_O'
all_species = ct.Species.list_from_file(input_file + ".yaml")
species = []

# Filter species
for S in all_species:
#    comp = S.composition
#    if 'C' in comp and 'H' in comp:
#        # Exclude all hydrocarbon species
#        continue
#    if 'N' in comp and comp != {'N': 2}:
#        # Exclude all nitrogen compounds except for N2
#        continue
#    if 'Ar' in comp:
#        # Exclude Argon
#        continue

    name = S.name

    if name == 'O2':
        species.append(S)

    if name == 'CO2':
        species.append(S)

species_names = {S.name for S in species}
print('Species: {0}'.format(', '.join(S.name for S in species)))

# Filter reactions, keeping only those that only involve the selected species
ref_phase = ct.Solution(thermo='ideal-gas', kinetics='gas', species=all_species)
all_reactions = ct.Reaction.list_from_file(input_file + ".yaml", ref_phase)
reactions = []

print('\nReactions:')
for R in all_reactions:
    if not all(reactant in species_names for reactant in R.reactants):
        continue

    if not all(product in species_names for product in R.products):
        continue

    reactions.append(R)
    print(R.equation)
print('\n')

gas1 = ct.Solution(input_file + ".yaml")
gas2 = ct.Solution(name=input_file + "-submech",
                   thermo="ideal-gas", kinetics="gas",
                   transport_model="mixture-averaged",
                   species=species, reactions=reactions)

# Save the resulting mechanism for later use
gas2.update_user_header({"description": "Submechanism extracted from " + input_file + ".yaml"})
gas2.write_yaml(input_file + "-submech.yaml", header=True)
