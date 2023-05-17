import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

width = .6
filt = False
sasa_filt = False
noseq = 0

print("\nSCRIPT FOR ALLOSTERY ANALYSIS IS LAUNCHED\n")

if not ("-f_all" in sys.argv and "-allsn" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 11):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -noseq noseq -filt -sasa_filt\n")
    exit()
    
for i in range(1, len(sys.argv)) :
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
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
    if sys.argv[i] == "":
        noseq = int(sys.argv[i+1])
    if sys.argv[i] == "-filt":
        filt = True
    if sys.argv[i] == "-sasa_filt":
        sasa_filt = True
    
if sasa_filt:
    map_path =  "output/" + name + "/map/" + name + "_sasa_filt_map"
    out_path =  "output/" + name + "/analysis/" + name + "_sasa_filt"
elif filt:
    map_path =  "output/" + name + "/map/" + name + "_filt_map"
    out_path =  "output/" + name + "/analysis/" + name + "_filt"
else:
    map_path =  "output/" + name + "/map/" + name + "_map"
    out_path =  "output/" + name + "/analysis/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
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
    print("Error reading file", act_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", act_path,"by '" + str(as_name) + "' key\n")
    exit()

intensity = []
for i in range(NResidues):
    if i + 1 in active_site:
        intensity.append(-1)
    else:
        inten = 0
        for resid in active_site:
            if np.abs(resid - 1 - i) >= noseq:
                inten += map_[resid - 1][i]
        intensity.append(inten)

try:
    with open(all_path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", all_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
    exit()
    
try:
    all_site = np.array(your_data[alls_name])
except:
    print("Error: Can't get data from file", all_path,"by '" + str(alls_name) + "' key\n")
    exit()
    
NAct = len(active_site)
NAll = len(all_site)

act_intensity = np.zeros(NAct)
all_intensity = np.zeros(NAll)
act_entropy = map_[active_site-1, active_site-1]
all_entropy = map_[all_site-1, all_site-1]
act_Entropy = 0
all_Entropy = 0
connectivity = 0
act_names = []
all_names = []
for i, res in enumerate(active_site):
    act_names.append(names[res-1] + "\n(" + str(real_numbers[res-1]) +")")
    for j in range(NResidues):
        if not (j+1 in active_site):
            act_intensity[i] += map_[res-1, j]
        elif (j+1 == res):
            act_Entropy += map_[res-1, j]
        else:
            act_Entropy += 1/2*map_[res-1, j]
            
        if (j+1 in all_site):
            connectivity += map_[res-1, j]
            
for i, res in enumerate(all_site):
    all_names.append(names[res-1] + "\n(" + str(real_numbers[res-1]) +")")
    for j in range(NResidues):
        if not (j+1 in all_site):
            all_intensity[i] += map_[res-1, j]
        elif (j+1 == res):
            all_Entropy += map_[res-1, j]
        else:
            all_Entropy += 1/2*map_[res-1, j]
            
act_mi = np.sum(act_intensity)
all_mi = np.sum(all_intensity)

size = np.max([NAct, NAll])
        
print("The mutual information of the active site with the remaining protein is", round(act_mi, 2))
print("The mutual information of the allosteric site with the remaining protein is", round(all_mi, 2))
print("The internal two-dimensional entropy of the active site is", round(act_Entropy, 2))
print("The internal two-dimensional entropy of the allosteric site is", round(all_Entropy, 2))
print("Allosteric connectivity between active and allosteric sites is", round(connectivity, 2), '\n')

fig = plt.figure(figsize=(size*width + 5, 24))
axs = [None for _ in range(4)]

axs[0] = plt.subplot2grid((2,2), (0,0),colspan=1, rowspan=1)
INTENSITY = {"Intensity": np.array(act_intensity), "Residue": act_names}
colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)
axs[0] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[0])
axs[0].set_title('Intensity of connectivity of residues from the active site for ' + name, fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=16)

axs[1] = plt.subplot2grid((2,2), (0,1),colspan=1, rowspan=1)
INTENSITY = {"Intensity": np.array(all_intensity), "Residue": all_names}
colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)
axs[1] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[1])
axs[1].set_title('Intensity of connectivity of residues from the allosteric site for ' + name, fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=16)

axs[2] = plt.subplot2grid((2,2), (1,0),colspan=1, rowspan=1)
INTENSITY = {"Intensity": np.array(act_entropy), "Residue": act_names}
colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)
axs[2] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[2])
axs[2].set_title('Two-dimensional entropy of residues from the active site for ' + name, fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=16)

axs[3] = plt.subplot2grid((2,2), (1,1),colspan=1, rowspan=1)
INTENSITY = {"Intensity": np.array(all_entropy), "Residue": all_names}
colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)
axs[3] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[3])
axs[3].set_title('Two-dimensional entropy of residues from the allosteric site for ' + name, fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=16)

try:
    fig.savefig(out_path + '_allostery.pdf')
    print("File",out_path + "_allostery.pdf created")
except:
    print("Error writing file",out_path + '_allostery.pdf')
