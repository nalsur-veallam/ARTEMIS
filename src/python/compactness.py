import sys
import numpy as np
import json

print("\nSCRIPT FOR CALCULATING COMPACTNESS IS LAUNCHED\n")

p_charge = ["ARG", "LYS"]
n_charge = ["ASP", "GLU"]

if not ("-n" in sys.argv and len(sys.argv) == 3):
    print("USAGE:\n"+sys.argv[0]+" -n name\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
        
        
map_path =  "output/" + name + "/map/" + name + "_map"

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
     
try:
    names = np.array(data['names'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key\n")
    exit()
   
NResidues = int(np.shape(map_)[0])

mi = 0
repulsive_mi = 0
   
for i in range(NResidues):
    
    for j in range(i+1, NResidues):
        mi += map_[i, j]
        if (names[i] in p_charge and names[j] in p_charge) or (names[i] in n_charge and names[j] in n_charge):
            repulsive_mi += map_[i, j]

print("Full MI is equal to ", round(mi, 3))
print("Repulsive MI is equal to ", round(repulsive_mi, 3))
print("Compactness is equal to ", round((mi - repulsive_mi)/mi, 3))
print("Normed compactness ( (mi - repulsive_mi)/NRes^2) is equal to ", round((mi - repulsive_mi)/NResidues**2, 3), '\n')
