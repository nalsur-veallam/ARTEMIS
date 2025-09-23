import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import zscore

width = .6

def search_critical(Allostery, out_path, noseq):

    print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY CALCULATION IS LAUNCHED\n")

    mask = np.ones(Allostery.NResidues)
    for resid in Allostery.active_site:
        for i in range(int(max(0, resid-1-noseq)), int(min(Allostery.NResidues, resid+noseq))):
            mask[i] = 0

    for resid in Allostery.allosteric_site:
        for i in range(int(max(0, resid-1-noseq)), int(min(Allostery.NResidues, resid+noseq))):
            mask[i] = 0

    intensity = []
    for i in range(Allostery.NResidues):
        if i + 1 in Allostery.active_site or mask[i] == 0 or i + 1 in Allostery.allosteric_site:
            intensity.append(0)
        else:
            mi_act = 0
            mi_all = 0
            entr = Allostery.map_[i][i]
            for resid in Allostery.active_site:
                mi_act += Allostery.map_[resid - 1][i]
            for resid in Allostery.allosteric_site:
                mi_all += Allostery.map_[resid - 1][i]
            intensity.append(mi_all*mi_act/entr)

    new_names = []
    for i in range(Allostery.NResidues):
        new_names.append(Allostery.names[i] + "\n(" + str(Allostery.real_numbers[i]) +")")

    unique, counts = np.unique(new_names, return_counts=True)

    for i, item in enumerate(unique):
        if counts[i] > 1:
            indeces = np.argwhere(np.array(new_names)==item)
            for j, idx in enumerate(indeces):
                new_names[int(idx[0])] += "[" +str(j+1)+"]"

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

    if out_path == "allosteric_intensity.pdf":
        fname = "allostery.json"
    else:
        fname = ''.join(out_path.split('.')[:-1]) + '.json'
    Allostery.write(fname)
