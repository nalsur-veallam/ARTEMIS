import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore

width = .6
filt = False
sasa_filt = False
top = None
Zscore = False
table  = False
noseq = 0

def max_top(array, top):
    size = NResidues
    top10 = int(top/100*size)

    supp = np.zeros(size)
    for i in range(top10):
        supp[i] = 1

    d = {'data':array, 'index':np.arange(0,NResidues)}

    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['data'], ascending=False)

    for i in range(NResidues):
        if np.array(df['data'])[i] < 0:
            supp[i] = 1

    df['data'] = df['data']*supp
    df = df.sort_values(by=['index'])
    return np.array(df['data'])

print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

if not ("-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -sasa_filt -filt -asn active_site_name -noseq num_of_res(default 0) -top top -zscore -table\n")
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
    if sys.argv[i] == "-zscore":
        Zscore = True
    if sys.argv[i] == "-table":
        table = True
    if sys.argv[i] == "-top":
        top = float(sys.argv[i+1])
        if top <=0 or top >= 100:
            print("Error: the top parameter must be in the range (0, 100)\n")
            exit()
    
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

if top is None:
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

    if Zscore:
        fig, axs = plt.subplots(figsize=(20, 15), constrained_layout=True)
        axs.hist(intensity, bins=50)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.title('Histogram for intensity of connectivity of residues\nwith the active site for ' + name, fontsize=40)
        try:
            fig.savefig(out_path + '_intensity_hist.pdf')
            print("File",out_path + "_intensity_hist.pdf created")
        except:
            print("Error writing file",out_path + '_intensity_hist.pdf')

        intensity = zscore(intensity)
else:
    intensity = np.zeros(NResidues)

    for resid in active_site:
        inten = []
        for i in range(NResidues):
            if i+1 in active_site:
                inten.append(0)
            elif np.abs(resid - 1 - i) >= noseq:
                inten.append(map_[resid - 1][i])
            else:
                inten.append(0)
        intensity += max_top(inten, top)

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(real_numbers[i]) +")")

unique, counts = np.unique(new_names, return_counts=True)

for i, item in enumerate(unique):
    if counts[i] > 1:
        indeces = np.argwhere(np.array(new_names)==item)
        for j, idx in enumerate(indeces):
            new_names[int(idx[0])] += "[" +str(j+1)+"]"

if top is None and not Zscore:

    INTENSITY = {}
    INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
    INTENSITY["Residue"] = new_names

    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()

    fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
    axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
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
    new_data['real_numbers'] = data['real_numbers']
    try:
        with open(out_path + '_intensity.json', 'w') as outfile:
            json.dump(new_data, outfile)
        print("File",out_path + "_intensity.json created\n")
    except:
        print("Error writing file",out_path + '_intensity.json\n')

elif top is None:

    INTENSITY = {}
    INTENSITY["Z-score intensity"] = np.array(intensity)
    INTENSITY["Residue"] = new_names

    colarr = INTENSITY["Z-score intensity"] - np.min(INTENSITY["Z-score intensity"])
    colors = plt.cm.viridis(colarr/np.max(colarr))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()

    #INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
    fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
    axs = sns.barplot(x="Residue", y="Z-score intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
    plt.tick_params(axis='both', which='major', labelsize=16)

    plt.title('Intensity of connectivity of residues with the active site zscore for ' + name, fontsize=40)
    try:
        fig.savefig(out_path + '_intensity_zscore.pdf')
        print("File",out_path + "_intensity_zscore.pdf created")
    except:
        print("Error writing file",out_path + '_intensity_zscore.pdf')

    new_data = {}
    new_data['name'] = name
    new_data['names'] = data['names']
    new_data['NResidues'] = NResidues
    new_data['active_site'] = your_data[as_name]
    new_data['Z-score intensity'] = intensity.tolist()
    try:
        with open(out_path + '_intensity_zscore.json', 'w') as outfile:
            json.dump(new_data, outfile)
        print("File",out_path + "_intensity_zscore.json created\n")
    except:
        print("Error writing file",out_path + '_intensity_zscore.json\n')

else:

    INTENSITY = {}
    INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
    INTENSITY["Residue"] = new_names

    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()

    #INTENSITY = pd.DataFrame(data=np.array([intensity, np.linspace(1, NResidues, NResidues)]).T, columns=["Intensity", "Residue"])
    fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
    axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
    plt.tick_params(axis='both', which='major', labelsize=16)

    plt.title('Intensity of connectivity of residues with the active site for ' + name, fontsize=40)
    try:
        fig.savefig(out_path + '_intensity_top'+str(top)+'.pdf')
        print("File",out_path + "_intensity_top"+str(top)+".pdf created")
    except:
        print("Error writing file",out_path + '_intensity_top'+str(top)+'.pdf')

    new_data = {}
    new_data['name'] = name
    new_data['names'] = data['names']
    new_data['NResidues'] = NResidues
    new_data['active_site'] = your_data[as_name]
    new_data['intensity'] = intensity.tolist()
    new_data['top'] = top
    try:
        with open(out_path + '_intensity_top'+str(top)+'.json', 'w') as outfile:
            json.dump(new_data, outfile)
        print("File",out_path + "_intensity_top"+str(top)+".json created\n")
    except:
        print("Error writing file",out_path + '_intensity_top'+str(top)+'.json\n')


    
if table:
    act_s_name = []
    for i in active_site:
        act_s_name.append(names[i-1] + " (" + str(real_numbers[i-1]) +")")

    unique, counts = np.unique(act_s_name, return_counts=True)

    for i, item in enumerate(unique):
        if counts[i] > 1:
            indeces = np.argwhere(np.array(act_s_name)==item)
            for j, idx in enumerate(indeces):
                act_s_name[int(idx[0])] += "[" +str(j+1)+"]"

    int_act_site = np.zeros((NResidues, len(active_site)))

    intensity = []
    for i in range(NResidues):

        if i + 1 in active_site:
            intensity.append(0)
        else:
            inten = 0
            k = 0
            for resid in active_site:
                if np.abs(resid - 1 - i) >= noseq:
                    inten += map_[resid - 1][i]
                    int_act_site[i, k] = map_[resid - 1][i]
                k += 1
            intensity.append(inten)

    new_names = []
    for i in range(NResidues):
        new_names.append(names[i] + " (" + str(real_numbers[i]) +")")

    unique, counts = np.unique(new_names, return_counts=True)

    for i, item in enumerate(unique):
        if counts[i] > 1:
            indeces = np.argwhere(np.array(new_names)==item)
            for j, idx in enumerate(indeces):
                new_names[int(idx[0])] += "[" +str(j+1)+"]"

    INTENSITY = {}
    INTENSITY["Residue"] = new_names
    INTENSITY["Intensity"] = np.array(intensity)
    INTENSITY["Mean Intensity"] = np.mean(int_act_site, axis=1)
    INTENSITY["STD"] = np.std(int_act_site, axis=1)
    INTENSITY = pd.DataFrame(data=INTENSITY)

    int_act_site = pd.DataFrame(data=int_act_site, columns=act_s_name)
    INTENSITY = pd.concat([INTENSITY, int_act_site], axis=1)

    INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

    INTENSITY = np.round(INTENSITY, 2)

    try:
        INTENSITY.to_csv(out_path + "_intensity_table.csv")
        print("File",out_path + "_intensity_table.csv created\n")
    except:
        print("Error writing file",out_path + '_intensity_table.csv\n')
    
