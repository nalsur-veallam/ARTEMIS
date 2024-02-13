import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from artemis.datatypes import MAP
import json

diag = True
norm = False

def draw_map(Map, out_path, diag, norm):

    print("\nSCRIPT FOR DRAWING MI MATRIX IS LAUNCHED\n")

    labels = []
    for i in range(Map.NResidues):
        labels.append(Map.names[i] + " (" + str(Map.real_numbers[i]) + ")")
        if not diag:
            Map.map_[i,i] = np.nan

    if norm:
        Map.map_ = Map.map_ / np.nanmax(Map.map_)

    MIE = pd.DataFrame(data=Map.map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    sns.heatmap(MIE, annot=False, cmap="icefire")

    plt.title('Mutual information on residues', fontsize=20)
    try:
        fig.savefig(out_path)
        print("File",out_path," created\n")
    except:
        print("Error writing file",out_path + '\n')
