import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore

width = .6

if not ("-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name")
    exit()
    
for i in range(1, 7) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
        
if not (name and as_name and path):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name")
    exit()
    
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
NResidues = data['NResidues']

with open(path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

intensity = []
for i in range(NResidues):
        if i + 1 in active_site:
            intensity.append(0)
        else:
            inten = 0
            for resid in active_site:
                inten += map_[resid - 1][i]
            intensity.append(inten)
            
fig, axs = plt.subplots(figsize=(20, 15), constrained_layout=True)
axs.hist(intensity, bins=50)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.title('Histogram for intensity of connectivity of residues\nwith the active site for ' + name, fontsize=40)
fig.savefig(out_path + '_intensity_hist.pdf')
            
intensity = zscore(intensity)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(i+1) +")")

INTENSITY = {}
INTENSITY["Z-score intensity"] = np.array(intensity)
INTENSITY["Residue"] = new_names


#INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Z-score intensity", data=INTENSITY, palette="viridis", dodge=False, hue="Z-score intensity")
axs.legend_.remove()
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Intensity of connectivity of residues with the active site zscore for ' + name, fontsize=40)
fig.savefig(out_path + '_intensity_zscore.pdf')

new_data = {}
new_data['name'] = name
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['active_site'] = your_data[as_name]
new_data['Z-score intensity'] = intensity.tolist()
with open(out_path + '_intensity_zscore.json', 'w') as outfile:
    json.dump(new_data, outfile)
