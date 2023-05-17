import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

diag = True
norm = False

rnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
charges = [0, +1, 0, -1, 0, -1, 0, 0, 0, 0, 0, +1, 0, 0, 0, 0, 0, 0, 0, 0]

print("\nSCRIPT FOR DRAWING CHARGE MATRIX IS LAUNCHED\n")

if not ("-n" in sys.argv and len(sys.argv) >= 3):
    print("USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nodiag":
        diag = False
    if sys.argv[i] == "-norm":
        norm = True
        
        
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/map/" + name + "_charge_map"

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
    exit()

try:
    names = np.array(data['names'])
    NResidues = len(names)
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key. Continues without this data.")
    exit()
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'real_numbers' key. Continues without this data.")
    real_numbers = np.empty(NResidues)

map_ = np.zeros((NResidues, NResidues))

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
    
Charges = np.zeros(NResidues)
for i in range(20):
    resname = rnames[i]
    Charges[np.argwhere(names == resname).reshape(-1)] = charges[i]
    
i, j = np.indices(map_.shape)
map_[i, j] = np.abs(Charges[i] * Charges[j])

if not diag:
        map_[i,i] = np.nan

if norm:
    map_ = map_ / np.nanmax(map_)
    
MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

plt.title('Charges product on residues for ' + name, fontsize=20)
try:
    fig.savefig(out_path + '.pdf')
    print("File",out_path + ".pdf created\n")
except:
    print("Error writing file",out_path + '.pdf\n')
    
new_data = {}
new_data['names'] = data['names']
new_data['NResidues'] = data["NResidues"]
new_data['name'] = name
new_data['map'] = map_.tolist()
new_data['real_numbers'] = data['real_numbers']
try:
    with open(out_path + '.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + ".json created\n")
except:
    print("Error writing file",out_path + '.json\n')
