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

def search_allostery(args):

    print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

    map_path =  args.files[0]
    out_path =  "./allosteric"

    path = args.files[1]

    as_name = 'act_s'

    try:
        with open(map_path) as json_file:
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
            plt.title('Histogram for intensity of connectivity of residues\nwith the active site', fontsize=40)
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
                new_names[int(idx)] += "[" +str(j+1)+"]"

    if top is None and not Zscore:

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
        INTENSITY["Residue"] = new_names

        colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
        colors = (sns.color_palette(colors, as_cmap=True)).tolist()

        fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
        axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
        plt.tick_params(axis='both', which='major', labelsize=16)

        plt.title('Intensity of connectivity of residues with the active site', fontsize=40)
        try:
            fig.savefig(out_path + '_intensity.pdf')
            print("File",out_path + "_intensity.pdf created")
        except:
            print("Error writing file",out_path + '_intensity.pdf')


