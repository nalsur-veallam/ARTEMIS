import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore

width = .6
noseq = 0

print("\nSCRIPT FOR Z-SCORE OF ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

if not ("-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)\n")
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
    with open(path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", path, ". USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", path,"by '" + str(as_name) + "' key\n")
    exit()

intensity = []
for i in range(NResidues):
        if i + 1 in active_site:
            intensity.append(0)
        else:
            inten = 0
            for resid in active_site:
                if np.abs(resid - 1 - i) >= noseq:
                    inten += map_[resid - 1][i]
            intensity.append(inten)
            
fig, axs = plt.subplots(figsize=(20, 15), constrained_layout=True)
axs.hist(intensity, bins=50)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.title('Histogram for intensity of connectivity of residues\nwith the active site for ' + name, fontsize=40)
try:
    fig.savefig(out_path + '_intensity_hist.pdf')
    print("File",out_path + "_intensity_hist.pdf created")
except:
    print("Error writing file",out_path + '_intensity_hist.pdf')
            
intensity = zscore(intensity)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(real_numbers[i]) +")")

INTENSITY = {}
INTENSITY["Z-score intensity"] = np.array(intensity)
INTENSITY["Residue"] = new_names

colarr = INTENSITY["Z-score intensity"] - np.min(INTENSITY["Z-score intensity"])
colors = plt.cm.viridis(colarr/np.max(colarr))
colors = sns.color_palette(colors, as_cmap=True)

#INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Z-score intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Intensity of connectivity of residues with the active site zscore for ' + name, fontsize=40)
try:
    fig.savefig(out_path + '_intensity_zscore.pdf')
    print("File",out_path + "_intensity_zscore.pdf created")
except:
    print("Error writing file",out_path + '_intensity_zscore.pdf')

new_data = {}
new_data['name'] = name
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['active_site'] = your_data[as_name]
new_data['Z-score intensity'] = intensity.tolist()
try:
    with open(out_path + '_intensity_zscore.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + "_intensity_zscore.json created\n")
except:
    print("Error writing file",out_path + '_intensity_zscore.json\n')
