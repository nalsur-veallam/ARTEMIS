import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

width = .6

if not (sys.argv[1] == "-n" and len(sys.argv) == 3):
    print("USAGE:\n"+sys.argv[0]+" -n name")
    exit()
    
name = sys.argv[2]
map_path = "output/map/" + name + "_map"
out_path = "output/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])
NResidues = data['NResidues']

labels = []
for i in range(len(names)):
    labels.append(names[i] + "\n(" + str(real_numbers[i]) + ")")

intensity = np.sum(map_[::-1, :], axis=0) - 2*np.diagonal(map_[::-1, :])

INTENSITY = {}
INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
INTENSITY["Residue"] = labels


fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette="viridis", dodge=False, hue="Intensity")
axs.legend_.remove()
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('The importance of residues for ' + name, fontsize=40)
fig.savefig(out_path + '_crit.pdf')
