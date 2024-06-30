import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore

width = .6

def max_top(array, top):
    size = len(array)
    top10 = int(top/100*size)

    supp = np.zeros(size)
    for i in range(top10):
        supp[i] = 1

    d = {'data':array, 'index':np.arange(0,size)}

    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['data'], ascending=False)

    for i in range(size):
        if np.array(df['data'])[i] < 0:
            supp[i] = 1

    df['data'] = df['data']*supp
    df = df.sort_values(by=['index'])
    return np.array(df['data'])

def search_allostery(Allostery, out_path, top, noseq, Zscore, Clustering=None, ClusterIdx=None):

    print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

    if Clustering is None or ClusterIdx is None:

        mask = np.ones(Allostery.NResidues)
        for resid in Allostery.active_site:
            for i in range(int(max(0, resid-1-noseq)), int(min(Allostery.NResidues, resid+noseq))):
                mask[i] = 0
    else:

        mask = np.zeros(Allostery.NResidues)
        mask[np.argwhere(Clustering==ClusterIdx)] = 1

        for resid in Allostery.active_site:
            for i in range(max(0, resid-1-noseq), min(Allostery.NResidues, resid+noseq)):
                mask[i] = 0

    if top is None:
        intensity = []
        for i in range(Allostery.NResidues):
            if i + 1 in Allostery.active_site or mask[i] == 0:
                intensity.append(0)
            else:
                inten = 0
                for resid in Allostery.active_site:
                    inten += Allostery.map_[resid - 1][i]
                intensity.append(inten)
    else:
        intensity = np.zeros(Allostery.NResidues)

        for resid in Allostery.active_site:
            inten = []
            for i in range(Allostery.NResidues):
                if i+1 in Allostery.active_site or mask[i] == 0:
                    inten.append(0)
                else:
                    inten.append(Allostery.map_[resid - 1][i])
            intensity += max_top(inten, top)

    new_names = []
    for i in range(Allostery.NResidues):
        new_names.append(Allostery.names[i] + "\n(" + str(Allostery.real_numbers[i]) +")")

    unique, counts = np.unique(new_names, return_counts=True)

    for i, item in enumerate(unique):
        if counts[i] > 1:
            indeces = np.argwhere(np.array(new_names)==item)
            for j, idx in enumerate(indeces):
                new_names[int(idx[0])] += "[" +str(j+1)+"]"

    if top is None and not Zscore:

        Allostery.intensity = np.array(intensity)

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
        INTENSITY["Residue"] = new_names

        colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
        colors = (sns.color_palette(colors, as_cmap=True)).tolist()

        fig, axs = plt.subplots(figsize=(Allostery.NResidues*width, 20), constrained_layout=True)
        axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
        plt.tick_params(axis='both', which='major', labelsize=16)

        plt.title('Intensity of connectivity of residues with the active site', fontsize=40)
        try:
            fig.savefig(out_path)
            print("File",out_path + " created")
        except:
            print("Error writing file",out_path)


    elif top is None:

        intensity = zscore(intensity)

        Allostery.intensity = np.array(intensity)

        INTENSITY = {}
        INTENSITY["Z-score intensity"] = np.array(intensity)
        INTENSITY["Residue"] = new_names

        colarr = INTENSITY["Z-score intensity"] - np.min(INTENSITY["Z-score intensity"])
        colors = plt.cm.viridis(colarr/np.max(colarr))
        colors = (sns.color_palette(colors, as_cmap=True)).tolist()

        fig, axs = plt.subplots(figsize=(Allostery.NResidues*width, 20), constrained_layout=True)
        axs = sns.barplot(x="Residue", y="Z-score intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
        plt.tick_params(axis='both', which='major', labelsize=16)

        plt.title('Intensity of connectivity of residues with the active site zscore', fontsize=40)
        try:
            fig.savefig(out_path)
            print("File",out_path + " created")
        except:
            print("Error writing file",out_path)

    else:

        Allostery.intensity = np.array(intensity)

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity) / np.sqrt(sum(abs(np.array(intensity).flatten())**2))
        INTENSITY["Residue"] = new_names

        colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
        colors = (sns.color_palette(colors, as_cmap=True)).tolist()

        fig, axs = plt.subplots(figsize=(Allostery.NResidues*width, 20), constrained_layout=True)
        axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, hue="Residue", legend=False)
        plt.tick_params(axis='both', which='major', labelsize=16)

        plt.title('Intensity of connectivity of residues with the active site top ' +str(top) + '%', fontsize=40)
        try:
            fig.savefig(out_path)
            print("File",out_path + " created")
        except:
            print("Error writing file",out_path)

    if out_path == "allosteric_intensity.pdf":
        fname = "allostery.json"
    else:
        fname = ''.join(out_path.split('.')[:-1]) + '.json'
    Allostery.write(fname)
