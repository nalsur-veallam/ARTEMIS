import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json
from pymol import cmd, stored

max_sasa = [257.7, 348.8, 287, 264.8, 268, 272.8, 316.9, 227.3, 340.4, 283.1, 308.6, 365.5, 309.1, 367.7, 285.3, 244.1, 267.1, 397.4, 349, 264.8]
rnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

if not ("-sasa" in sys.argv and "-cutoff" in sys.argv and "-strc" in sys.argv and sys.argv[1] == "-n" and len(sys.argv) == 9):
    print("USAGE:\n"+sys.argv[0]+" -n name -strc sctructure.pdb(.gro ...) -cutoff cutoff -sasa sasa_file.xvg")
    exit()
    
for i in range(1, 9) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-sasa":
        sasa_path = sys.argv[i+1]
    if sys.argv[i] == "-cutoff":
        cutoff = float(sys.argv[i+1])
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/map/" + name + "_sasa_filt_map"

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])
NResidues = int(data["NResidues"])

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
    map_[i,i] = np.nan

cmd.load(str_path)
stored.resnames = []
cmd.iterate('name ca', 'stored.resnames.append(resn)')
resnames = np.array(stored.resnames)

f = open(sasa_path, 'r')
skip = 0
for line in f:
    if line[0] == '#' or line[0] == '@':
        skip += 1
    else:
        break
    
f.close()
    
sasa_file = np.loadtxt(sasa_path, skiprows=skip)

sasa_per_residue = 100*sasa_file[:, 1]
max_sasa_per_residue = np.zeros(len(sasa_per_residue))

for i in range(20):
    resname = rnames[i]
    max_sasa_per_residue[np.argwhere(resnames == resname).reshape(-1)] = max_sasa[i]

exp = np.greater(np.array(sasa_per_residue), cutoff*np.array(max_sasa_per_residue)).astype(np.int8)

i, j = np.indices(map_.shape)
new_map_ = map_.copy()
new_map_[i,j] = map_[i, j] * exp[i] * exp[j]
map_[i, j] = map_[i, j] * exp[j]

MIE = pd.DataFrame(data=new_map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

plt.title('Mutual information on residues for ' + name + '\n with filtration', fontsize=20)
fig.savefig(out_path + '.pdf')

new_data = {}
new_data['names'] = data['names']
new_data['NResidues'] = data["NResidues"]
new_data['name'] = name
new_data['map'] = map_.tolist()
new_data['real_numbers'] = data['real_numbers']
new_data['cutoff'] = cutoff
with open(out_path + '.json', 'w') as outfile:
    json.dump(new_data, outfile)
