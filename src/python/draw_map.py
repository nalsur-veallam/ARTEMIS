import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

diag = True

if not ("-n" in sys.argv and len(sys.argv) >= 3):
    print("USAGE:\n"+sys.argv[0]+" -n name -nodiag")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nodiag":
        diag = False
        
        
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/map/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
    if not diag:
        map_[i,i] = np.nan

MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

plt.title('Mutual information on residues for ' + name, fontsize=20)
fig.savefig(out_path + '.pdf')
