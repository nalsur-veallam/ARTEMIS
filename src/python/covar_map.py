import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
from pymol import cmd, stored
import json
import sys

All = True

print("\nCOVAR (DAT) GROMACS MATRIX CONVERTING SCRIPT IS LAUNCHED\n")

if not ("-f" in sys.argv and "-strc" in sys.argv and "-o" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f dat_file -o out_path\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f":
        i_file = sys.argv[i+1]
    if sys.argv[i] == "-o":
        o_file = sys.argv[i+1]
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]

symbols = []
counts = []

try:
    data = np.loadtxt(i_file)
except:
    print("Error reading file", i_file, ". USAGE:\n"+sys.argv[0]+" -f dat_file -o out_path\n")
    exit()

cmd.load(str_path)
cmd.h_add('all')
stored.residues = []
stored.reschs = []
stored.resnames = []
cmd.iterate('all', 'stored.residues.append(resi)')
cmd.iterate('all', 'stored.reschs.append(chain)')
cmd.iterate('all', 'stored.resnames.append(resn)')

NAtoms = len(stored.resnames)
NStr = len(data)

if 3*NAtoms**2 == NStr:
    All = True
elif 3*(NAtoms + 1)**2 == NStr:
    All = False
    NAtoms += 1
    stored.residues.append(stored.residues[-1])
    stored.reschs.append(stored.reschs[-1])
    stored.resnames.append(stored.resnames[-1])
elif 3*(NAtoms - 1)**2 == NStr:
    All = False
    NAtoms -= 1
    stored.residues.pop(-1)
    stored.reschs.pop(-1)
    stored.resnames.pop(-1)
else:
    print("Error! The structure does not fit the data: different number of degrees of freedom\n")
    exit()

NRes = 0
resi = ''
chain = ''
resn = ''
residue = []
labels = []
for res, ch, rn in zip(stored.residues, stored.reschs, stored.resnames):
    if res != resi or rn != resn or ch != chain:
        resi = res
        resn = rn
        chain = ch
        NRes += 1
        labels.append(resn + " ("+str(NRes) + ")")
        
    residue.append(NRes)
    
map_ = np.zeros((NRes, NRes))
residue = np.array(residue)

for i in range(NRes):
    idxs1 = np.argwhere(residue==i+1)
    for j in range(i+1, NRes):
        idxs2 = np.argwhere(residue==j+1)
        idx1 = np.ones(np.shape(idxs2)).reshape(-1, 1)@idxs1.reshape(1, -1)
        idx2 = idxs2.reshape(-1, 1)@np.ones(np.shape(idxs1)).reshape(1, -1)
        idx = (3*NAtoms*idx1+idx2).flatten().astype(int)
        map_[i, j] = np.sum(data[idx, :].flatten()) + np.sum(data[idx+int(NAtoms), :].flatten()) + np.sum(data[idx+int(2*NAtoms), :].flatten())
        map_[j, i] = map_[i, j]
    

DM = pd.DataFrame(data=map_[::-1, :], index=labels[::-1], columns=labels)


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(DM, annot=False, cmap="icefire")

plt.title('Covariance matrix', fontsize=20)
try:
    fig.savefig(o_file + '.pdf')
    print("File",o_file+".pdf created")
except:
    print("Error writing file", o_file+'.pdf')

new_data = {}
new_data['map'] = map_.tolist()
try:
    with open(o_file+'.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",o_file+".json created\n")    
except:
    print("Error writing file", o_file+'.json\n')
    exit()

