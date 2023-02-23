import sys
import numpy as np
import pandas as pd
import json

width = .6
noseq=0

print("\nSCRIPT FOR CONSTRUCTION OF THE ALLOSTERIC COMMUNICATION INTENSITY TABLE IS LAUNCHED\n")

if not ("-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f_act":
        act_path = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])    
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
    
try:
    NResidues = int(data['NResidues'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'NResidues' key. Continues without this data.")
    NResidues = len(map_)
    
try:
    names = np.array(data['names'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key. Continues without this data.")
    names = np.arange(1, NResidues+1)
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'real_numbers' key. Continues without this data.")
    real_numbers = np.empty(NResidues)


try:
    with open(act_path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", act_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", act_path,"by '" + str(as_name) + "' key\n")
    exit()

act_s_name = []
for i in active_site:
    act_s_name.append(names[i-1] + " (" + str(real_numbers[i-1]) +")")

int_act_site = np.zeros((NResidues, len(active_site)))

intensity = []
for i in range(NResidues):
    
    if i + 1 in active_site:
        intensity.append(0)
    else:
        inten = 0
        k = 0
        for resid in active_site:
            if np.abs(resid - 1 - i) >= noseq:
                inten += map_[resid - 1][i]
                int_act_site[i, k] = map_[resid - 1][i]
            k += 1
        intensity.append(inten)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + " (" + str(real_numbers[i]) +")")

INTENSITY = {}
INTENSITY["Residue"] = new_names
INTENSITY["Intensity"] = np.array(intensity)
INTENSITY["Mean Intensity"] = np.mean(int_act_site, axis=1)
INTENSITY["STD"] = np.std(int_act_site, axis=1)
INTENSITY = pd.DataFrame(data=INTENSITY)

int_act_site = pd.DataFrame(data=int_act_site, columns=act_s_name)
INTENSITY = pd.concat([INTENSITY, int_act_site], axis=1)

INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

INTENSITY = np.round(INTENSITY, 2)        

try:
    INTENSITY.to_csv(out_path + "_intensity_table.csv")
    print("File",out_path + "_intensity_table.csv created\n")
except:
    print("Error writing file",out_path + '_intensity_table.csv\n')
