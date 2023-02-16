import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

width = 0.6
noseq = 0

def max10(array):
    size = len(array)
    top10 = int(0.1*size)
    
    supp = np.ones(size)
    for i in range(top10):
        supp[i] = 1
    
    d = {'data':array, 'index':np.arange(0,NResidues)}
    
    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['data'], ascending=False)
    df['data'] = df['data']*supp
    df = df.sort_values(by=['index'])
    return np.array(df['data'])
    

if not ("-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -noseq num_of_res(default 0)")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
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
real_numbers = np.array(data['real_numbers'])

with open(path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

intensity = np.zeros(NResidues)

for resid in active_site:
    inten = []
    for i in range(NResidues):
        if real_numbers[i] in active_site:
            inten.append(0)
        elif np.abs(resid - 1 - i) >= noseq:
            inten.append(map_[resid - 1][i])
    intensity += max10(inten)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(real_numbers[i]) +")")

INTENSITY = {}
INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
INTENSITY["Residue"] = new_names


#INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette="viridis", dodge=False, hue="Intensity")
axs.legend_.remove()
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Intensity of connectivity of residues with the active site for ' + name, fontsize=40)
fig.savefig(out_path + '_intensity_top10.pdf')

new_data = {}
new_data['name'] = name
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['active_site'] = your_data[as_name]
new_data['intensity'] = intensity.tolist()
with open(out_path + '_intensity_top10.json', 'w') as outfile:
    json.dump(new_data, outfile)
