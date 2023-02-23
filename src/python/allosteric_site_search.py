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

print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

if not ("-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -filt -sasa_filt -noseq num_of_res(default 0)\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-filt":
        filt = True
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
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
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -filt -sasa_filt -noseq num_of_res(default 0)\n")
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
    with open(path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", path, ". USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name -filt -sasa_filt -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", path,"by '" + str(as_name) + "' key\n")
    exit()

intensity = []
for i in range(NResidues):
        if i + 1 in active_site:
            intensity.append(0)
        else:
            inten = 0
            for resid in active_site:
                if np.abs(resid - 1 - i) >= noseq:
                    inten += map_[resid - 1][i]
            intensity.append(inten)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(real_numbers[i]) +")")

INTENSITY = {}
INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
INTENSITY["Residue"] = new_names

colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
colors = sns.color_palette(colors, as_cmap=True)

#INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Intensity of connectivity of residues with the active site for ' + name, fontsize=40)
try:
    fig.savefig(out_path + '_intensity.pdf')
    print("File",out_path + "_intensity.pdf created")
except:
    print("Error writing file",out_path + '_intensity.pdf')

new_data = {}
new_data['name'] = name
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['active_site'] = your_data[as_name]
new_data['intensity'] = intensity
new_data['filtration'] = filt
try:
    with open(out_path + '_intensity.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + "_intensity.json created\n")
except:
    print("Error writing file",out_path + '_intensity.json\n')
    
    
    
    
