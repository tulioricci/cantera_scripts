import glob
import numpy as np
import re
from operator import itemgetter
import os

os.system("mkdir -p centerline_data")

dummy = glob.glob("csv/*phi0.55*")[0]
f = open(dummy)
header = f.readline().strip().split(",")
f.close()

ii = 0
for var in header:
    if var == "X_CO":
       idx_CO = ii
    if var == "X_CO2":
       idx_CO2 = ii
    if var == "T":
       idx_T = ii
    ii += 1

for phi in ["0.55", "0.70", "0.85", "1.00", "1.15", "1.30"]:
    print(phi)
    outputFile = "centerline_data/output_phi" + phi + ".dat"
    f = open(outputFile, "w")

    f.write("scale_coeff reaction array_idx x T X_CO X_CO2\n")
    filelist = glob.glob("csv/*phi" + phi + "*")

    myList = []
    for file in filelist:
        # Extract the integer after 'S'
        s_match = re.search(r'S(\d{4})', file)
        s_value = int(s_match.group(1)) if s_match else None

        # Extract the float after 'coeff'
        coeff_match = re.search(r'coeff([0-9]+\.[0-9]+)', file)
        coeff_value = float(coeff_match.group(1)) if coeff_match else None

        data = np.loadtxt(file, skiprows=1,delimiter=",")
        x = data[:,0]
        idx = np.argmin(np.abs(x-0.005))
        myList.append([coeff_value, s_value, idx, x[idx], data[idx,idx_T],
                       data[idx,idx_CO], data[idx,idx_CO2]])

    mySortedList = sorted(myList, key=itemgetter(4))
    mySortedList.reverse()
    #mySortedList = sorted(myList, key=itemgetter(1))
    #mySortedList = sorted(mySortedList, key=itemgetter(0))
    for item in mySortedList:
        list_as_a_string = ' '.join(str(x) for x in item)
        f.write(list_as_a_string + "\n")
    f.close()

