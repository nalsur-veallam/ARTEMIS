import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

if not (sys.argv[1] == "-n" and len(sys.argv) == 3):
    print("USAGE:\n"+sys.argv[0]+" -n name")
    exit()
    
name = sys.argv[2]
map_path = "output/map/" + name + "_map"
out_path = "output/map/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")

MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

plt.title('Mutual information on residues for ' + name, fontsize=20)
fig.savefig(out_path + '.pdf')
