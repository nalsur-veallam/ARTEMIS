import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from artemis.datatypes import MAP
import json

Rnames = np.array(['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP'])
N = len(Rnames)

diag = True
norm = False

def contacts_map(Map, out_path, diag, norm, vmax=None):

    print("\nSCRIPT FOR DRAWING CONTACS MATRIX IS LAUNCHED\n")

    i, j = np.indices(Map.map_.shape)
    Map.map_[i==j] = 0

    totalMI = np.sum(Map.map_)/2
    totalPairs = Map.NResidues*(Map.NResidues-1)/2

    pairmap = np.zeros((N,N))
    paircounts = np.zeros((N,N))

    for i in range(N):
        for j in range(i, N):
            idxs1 = np.argwhere(np.array(Map.names) == Rnames[i]).reshape(-1)
            idxs2 = np.argwhere(np.array(Map.names) == Rnames[j]).reshape(-1)

            if i == j:
                paircounts[i, j] = len(idxs1)*(len(idxs1)+1)/2

                for idx1 in idxs1:
                    for idx2 in idxs2:
                        pairmap[i,j] += Map.map_[idx1, idx2]/2

            else:
                for idx1 in idxs1:
                    for idx2 in idxs2:
                        pairmap[i,j] += Map.map_[idx1, idx2]
                        paircounts[i, j] += 1

                        pairmap[j,i] += Map.map_[idx1, idx2]
                        paircounts[j, i] += 1

    inten_theor = paircounts/totalPairs
    inten_real = pairmap/totalMI

    np.seterr(divide='ignore', invalid='ignore')
    enrichment = np.where(inten_theor>0, inten_real/inten_theor, np.nan)

    if norm:
        enrichment = np.where(enrichment>=1, enrichment, np.nan)

    i,j = np.indices(enrichment.shape)

    if not diag:
        enrichment[i==j] = np.nan

    CONTACS = pd.DataFrame(data=enrichment[::-1, :], index=Rnames[::-1], columns=Rnames)

    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    if norm:
        if vmax is None:
            sns.heatmap(CONTACS, annot=False, cmap="Reds", vmin=1)
        else:
            sns.heatmap(CONTACS, annot=False, cmap="Reds", vmin=1, vmax=vmax)
    elif vmax is None:
        sns.heatmap(CONTACS, annot=False, cmap="coolwarm", vmin=0, center=1)
    else:
        sns.heatmap(CONTACS, annot=False, cmap="coolwarm", vmin=0, center=1, vmax=vmax)

    plt.title('Contacs enrichment matrix', fontsize=20)
    try:
        fig.savefig(out_path)
        print("File",out_path,"created\n")
    except:
        print("Error writing file",out_path + '\n')
