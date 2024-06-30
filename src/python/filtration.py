import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json
from pymol import cmd, stored
from tqdm import tqdm

print("\nSCRIPT FOR MI MATRIX FILTRATION WITH PYMOL IS LAUNCHED\n")

if not ("-cutoff" in sys.argv and "-strc" in sys.argv and sys.argv[1] == "-n" and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -n name -strc sctructure.pdb(.gro ...) -cutoff cutoff\n")
    exit()
    
for i in range(1, 7) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-cutoff":
        cutoff = float(sys.argv[i+1])
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/map/" + name + "_filt_map"

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -strc sctructure.pdb(.gro ...) -cutoff cutoff\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
    
try:
    NResidues = int(data['NResidues'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'NResidues' key. Continues without this data.\n")
    NResidues = len(map_)
    
try:
    names = np.array(data['names'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key. Continues without this data.\n")
    names = np.arange(1, NResidues+1)
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'real_numbers' key. Continues without this data.\n")
    real_numbers = np.empty(NResidues)

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
    map_[i,i] = np.nan

cmd.load(str_path)
stored.residues = []
stored.reschs = []
stored.resnames = []
cmd.iterate('name ca', 'stored.residues.append(resi)')
cmd.iterate('name ca', 'stored.reschs.append(chain)')
cmd.iterate('name ca', 'stored.resnames.append(resn)')

rnames, idxs = np.unique(stored.resnames, return_index=True)
rnumbers = np.array(stored.residues)[idxs]
rchs = np.array(stored.reschs)[idxs]
sasa = []

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 4)

sasa_per_residue = []
for i in tqdm(range(NResidues)):
    if stored.reschs[i] == '':
        sasa_per_residue.append(float(cmd.get_area('resi '+ str(stored.residues[i]))))
    else:
        sasa_per_residue.append(float(cmd.get_area('resi '+ str(stored.residues[i]) + ' and chain ' + str(stored.reschs[i]))))

for i in range(len(rnames)):
    if rchs[i] == '':
        cmd.remove('all and not (resi '+ str(rnumbers[i]) + ')')
    else:
        cmd.remove('all and not (resi '+ str(rnumbers[i]) + ' and chain ' + str(rchs[i]) + ')')
    cmd.set('dot_solvent', 1)
    cmd.set('dot_density', 4)
    sasa.append(float(cmd.get_area('all')))
    cmd.delete('all')
    cmd.load(str_path)
    
max_sasa_per_residue = []
for i in range(NResidues):
    resn = stored.resnames[i]
    max_sasa_per_residue.append(sasa[np.argwhere(rnames == resn)[0,0]])

exp = np.greater(np.array(sasa_per_residue), cutoff*np.array(max_sasa_per_residue)).astype(np.int8)

i, j = np.indices(map_.shape)
new_map_ = map_.copy()
new_map_[i,j] = map_[i, j] * exp[i] * exp[j]
map_[i, j] = map_[i, j] * exp[j]

MIE = pd.DataFrame(data=new_map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

plt.title('Mutual information on residues for ' + name + '\n with filtration', fontsize=20)
try:
    fig.savefig(out_path + '.pdf')
    print("\nFile",out_path + ".pdf created")
except:
    print("\nError writing file",out_path + '.pdf')

new_data = {}
new_data['names'] = data['names']
new_data['NResidues'] = data["NResidues"]
new_data['name'] = name
new_data['map'] = map_.tolist()
new_data['real_numbers'] = data['real_numbers']
new_data['cutoff'] = cutoff
try:
    with open(out_path + '.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + ".json created\n")
except:
    print("Error writing file",out_path + '.json\n')
