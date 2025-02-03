import yaml
with open('uiuc_20sp.yaml') as mechYaml:
     my_dict = yaml.safe_load(mechYaml)
mechYaml.close()

newYaml = open("new_mech.yaml","w") 

newYaml.write("description:\n")
string = my_dict["description"]
string = string.split("\n")
for line in string:
    newYaml.write("  " + line + "\n")

newYaml.write("\n")
newYaml.write("generator: " + my_dict["generator"] + "\n")

newYaml.write("\n")
newYaml.write("units: " + str(my_dict["units"]) + "\n")

newYaml.write("\n")
string = my_dict["phases"]
newYaml.write("phases:\n")
nphases = len(string)
if nphases != 1:
    sys.exit()
phasesKeys = string[0].keys()
for keys in phasesKeys:
    if keys == "name":
        newYaml.write("- name: " + str(string[0]["name"] + "\n"))
    elif keys == "state":
        newYaml.write("  state:\n")
        newYaml.write("    T: " + str(string[0]["state"]["T"]) + "\n")
        newYaml.write("    P: " + str(string[0]["state"]["P"]) + "\n")        
    elif keys == "elements":
        newYaml.write("  elements: " + str(string[0]["elements"]).replace("'","") + "\n")
    elif keys == "species":
        newYaml.write("  species: " + str(string[0]["species"]).replace("'","") + "\n")
    else:
        newYaml.write("  " + keys + ": " + str(string[0][keys]) + "\n" )
        
newYaml.write("\n")
species_dict = my_dict["species"]
newYaml.write("species:\n")
nspecies = len(species_dict)
for string in species_dict:
    speciesKeys = string.keys()
    for keys in speciesKeys:
        if keys == "name":
            newYaml.write("- name: " + str(string["name"] + "\n"))
        if keys == "composition":
            newYaml.write("  composition: " + str(string["composition"]).replace("'","") + "\n")
        if keys == "thermo":
            newYaml.write("    model: " + str(string["thermo"]["model"]) + "\n")
            newYaml.write("    temperature-ranges: " + str(string["thermo"]["temperature-ranges"]) + "\n")
            newYaml.write("    data:\n")
            newYaml.write("    - " + str(string["thermo"]["data"][0]) + "\n")
            newYaml.write("    - " + str(string["thermo"]["data"][1]) + "\n")
        if keys == "transport":
            newYaml.write("  transport:\n")
            for transKeys in string["transport"]:
                newYaml.write("    " + transKeys + ": " + str(string["transport"][transKeys]) + "\n")
        if keys == "note":
            newYaml.write("  note: " + str(string["note"] + "\n"))        

newYaml.close()
