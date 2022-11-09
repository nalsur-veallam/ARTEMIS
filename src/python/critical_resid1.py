import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

width = .6

if not ("-allsn" in sys.argv and "-f_all" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) == 11):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name")
    exit()
    
for i in range(1, 11) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f_act":
        act_path = sys.argv[i+1]
    if sys.argv[i] == "-allsn":
        alls_name = sys.argv[i+1]
    if sys.argv[i] == "-f_all":
        all_path = sys.argv[i+1]
        
if not (name and as_name and act_path and all_path and alls_name):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name")
    exit()
    
    
map_path = "output/map/" + name + "_map"
out_path = "output/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
NResidues = data['NResidues']

with open(act_path) as json_file:
    active_data = json.load(json_file)

active_site = np.array(active_data[as_name])

with open(all_path) as json_file:
    all_data = json.load(json_file)

allosteric_site = np.array(all_data[alls_name])

connectivity = []

for i in range(NResidues):
    if i + 1 in active_site or i + 1 in allosteric_site:
            connectivity.append(0)
    else:
        connectiv = 0
        for resid in active_site:
            connectiv += map_[resid - 1][i]
        for resid in allosteric_site:
            connectiv += map_[resid - 1][i]
        connectivity.append(connectiv)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(i+1) +")")

CONNECTIVITY = {}
CONNECTIVITY["Connectivity"] = np.array(connectivity) / np.sqrt(sum(abs(np.array(connectivity).flatten())**2))
CONNECTIVITY["Residue"] = new_names
            
            
fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Connectivity", data=CONNECTIVITY, palette="viridis", dodge=False, hue="Connectivity")
axs.legend_.remove()
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Connectivity1 of residues with the active and allosteric sites for ' + name, fontsize=40)
fig.savefig(out_path + '_connectivity1.pdf')

new_data = {}
new_data['name'] = name
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['active_site'] = active_data[as_name]
new_data['active_site'] = all_data[alls_name]
new_data['connectivity'] = connectivity
with open(out_path + '_connectivity1.json', 'w') as outfile:
    json.dump(new_data, outfile)
