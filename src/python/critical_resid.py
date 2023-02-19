import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

width = .6
noseq = 0

if not (sys.argv[1] == "-n" and len(sys.argv) >= 3):
    print("USAGE:\n"+sys.argv[0]+" -n name -noseq num_of_res(default 0)")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
        
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])
NResidues = data['NResidues']

labels = []
for i in range(len(names)):
    labels.append(names[i] + "\n(" + str(real_numbers[i]) + ")")


#####################################################################################

Max_mie = []
Max_idx = []
Mie = []
intensity = []

for i in range(NResidues):
    mie = 0
    max_mie = 0
    idx  = labels[i] + "\nwith\n" + labels[i]
    for j in range(NResidues):
        if i != j:
            mie += map_[i, j]
            
        if map_[i, j] > max_mie and np.abs(i-j) > noseq:
            max_mie = map_[i, j]
            idx = labels[i] + "\nwith\n" + labels[j]
            
    Max_mie.append(max_mie)
    Mie.append(mie)
    Max_idx.append(idx)
    intensity.append(map_[i,i])

INTENSITY = {}
INTENSITY["Intensity"] = np.array(intensity) #/ np.sqrt(sum(abs(np.array(intensity).flatten())**2))
INTENSITY["Residue"] = labels

colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('The importance of residues for ' + name, fontsize=40)
fig.savefig(out_path + '_entropy.pdf')


INTENSITY = {}
INTENSITY["Intensity"] = np.array(Max_mie)/np.array(intensity)
INTENSITY["Residue"] = np.array(Max_idx)

colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Max MIE of residues for ' + name, fontsize=40)
fig.savefig(out_path + '_maxMie.pdf')

INTENSITY = {}
INTENSITY["Intensity"] = np.array(Mie)
INTENSITY["Residue"] = labels

colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('MIE of residues for ' + name, fontsize=40)
fig.savefig(out_path + '_mie.pdf')
