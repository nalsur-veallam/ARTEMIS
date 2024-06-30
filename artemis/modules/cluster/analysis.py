import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import json

heigh = .1
colors = ['blue', 'green', 'red', 'cyan', 'yellow', 'magenta', 'salmon', 'lime', 'hotpink', 'orange', 'chartreuse', 'limegreen', 'olive', 'purple', 'teal', 'gray', 'pink', 'firebrick', 'chocolate', 'brown', 'wheat', 'violet', 'aquamarine', 'palegreen', 'lightblue', 'lightpink', 'skyblue', 'silver', 'gold']

def analyze_clustering(Clusters, noseq, out_path):

    print("\nSCRIPT FOR CLUSTER ANALYSIS IS LAUNCHED\n")

    MIi = np.zeros(Clusters.NClusters)
    MIi_ = np.zeros(Clusters.NClusters)
    MIe = np.zeros(Clusters.NClusters)
    MIe_ = np.zeros(Clusters.NClusters)
    H = np.zeros(Clusters.NClusters)
    MIE = np.zeros((Clusters.NClusters, Clusters.NClusters))
    intensity_act = np.zeros(Clusters.NClusters)
    intensity_all = np.zeros(Clusters.NClusters)
    intensity_act_ = np.zeros(Clusters.NClusters)
    intensity_all_ = np.zeros(Clusters.NClusters)
    clAct = np.zeros(Clusters.NClusters)
    clAll = np.zeros(Clusters.NClusters)
    clNames = []
    clColors = []
    quantity = np.zeros(Clusters.NClusters)

    for cl in range(Clusters.NClusters):
        quantity[cl] = np.sum(Clusters.clustering_labels == cl+1)

        if len(colors[cl]) > 4:
            clNames.append(colors[cl][:4] + " cl(" + str(cl+1)+")")
        else:
            clNames.append(colors[cl] + " cl(" + str(cl+1)+")")

        clColors.append(colors[cl])

        for res in Clusters.active_site:
            if Clusters.clustering_labels[res-1] == cl+1:
                clAct[cl] += 1
        for res in Clusters.allosteric_site:
            if Clusters.clustering_labels[res-1] == cl+1:
                clAll[cl] += 1

    for i in range(Clusters.NResidues):
        H[Clusters.clustering_labels[i]-1]+= Clusters.map_[i,i]
        for j in range(i, Clusters.NResidues):
            if Clusters.clustering_labels[i] == Clusters.clustering_labels[j]:
                MIi[Clusters.clustering_labels[i]-1]+= Clusters.map_[i,j]
                MIi_[Clusters.clustering_labels[i]-1]+= Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]
            else:
                MIe[Clusters.clustering_labels[i]-1]+= Clusters.map_[i,j]
                MIe[Clusters.clustering_labels[j]-1]+= Clusters.map_[i,j]
                MIe_[Clusters.clustering_labels[i]-1]+= Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]
                MIe_[Clusters.clustering_labels[j]-1]+= Clusters.map_[i,j]/quantity[Clusters.clustering_labels[j]-1]
                MIE[Clusters.clustering_labels[i]-1, Clusters.clustering_labels[j]-1] += Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]/quantity[Clusters.clustering_labels[j]-1]
                MIE[Clusters.clustering_labels[j]-1, Clusters.clustering_labels[i]-1] += Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]/quantity[Clusters.clustering_labels[j]-1]

            if j in Clusters.active_site and np.abs(j - i) >= noseq:
                intensity_act[Clusters.clustering_labels[i]-1] += Clusters.map_[i,j]
                intensity_act_[Clusters.clustering_labels[i]-1] += Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]

            if j in Clusters.allosteric_site and np.abs(j - i) >= noseq:
                intensity_all[Clusters.clustering_labels[i]-1] += Clusters.map_[i,j]
                intensity_all_[Clusters.clustering_labels[i]-1] += Clusters.map_[i,j]/quantity[Clusters.clustering_labels[i]-1]

    clColors = sns.color_palette(clColors, as_cmap=True)

    if Clusters.NClusters > 7:
        fig = plt.figure(figsize=(24, 36))
    else:
        fig = plt.figure(figsize=(16, 24))

    axs = [None for _ in range(13)]

    axs[0] = plt.subplot2grid((7,4), (0,0),colspan=2, rowspan=2)
    CLUSTERING = pd.DataFrame(data=MIE[::-1, :], index=np.array(clNames)[::-1], columns=np.array(clNames))
    axs[0] = sns.heatmap(CLUSTERING, cmap="OrRd", annot=True, ax=axs[0], cbar=False)
    axs[0].set_title('Normalized MI matrix for clusters')

    if Clusters.NClusters > 7:
        font = {'size'   : 20}
        plt.rc('font', **font)

    axs[1] = plt.subplot2grid((7,4), (0,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(quantity), "Names": clNames}
    axs[1] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[1], hue="Names", legend=False)
    axs[1].set_title('Quantity of residues in clusters' )
    axs[1].set_xlabel('')
    axs[1].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[2] = plt.subplot2grid((7,4), (1,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(H), "Names": clNames}
    axs[2] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[2], hue="Names", legend=False)
    axs[2].set_title('Informational entropy of clusters' )
    axs[2].set_xlabel('')
    axs[2].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[3] = plt.subplot2grid((7,4), (2,0),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(intensity_act) / np.sqrt(sum(abs(np.array(intensity_act).flatten())**2)), "Names": clNames}
    axs[3] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[3], hue="Names", legend=False)
    axs[3].set_title('Intensity of clusters with the active site' )
    axs[3].set_xlabel('')
    axs[3].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[4] = plt.subplot2grid((7,4), (2,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(intensity_act_) / np.sqrt(sum(abs(np.array(intensity_act_).flatten())**2)), "Names": clNames}
    axs[4] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[4], hue="Names", legend=False)
    axs[4].set_title('Intensity (normalized) of clusters with the active site' )
    axs[4].set_xlabel('')
    axs[4].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[5] = plt.subplot2grid((7,4), (3,0),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(intensity_all) / np.sqrt(sum(abs(np.array(intensity_all).flatten())**2)), "Names": clNames}
    axs[5] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[5], hue="Names", legend=False)
    axs[5].set_title('Intensity of clusters with the allosteric site' )
    axs[5].set_xlabel('')
    axs[5].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[6] = plt.subplot2grid((7,4), (3,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(intensity_all_) / np.sqrt(sum(abs(np.array(intensity_all_).flatten())**2)), "Names": clNames}
    axs[6] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[6], hue="Names", legend=False)
    axs[6].set_title('Intensity (normalized) of clusters with the allosteric site' )
    axs[6].set_xlabel('')
    axs[6].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[7] = plt.subplot2grid((7,4), (4,0),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(MIi), "Names": clNames}
    axs[7] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[7], hue="Names", legend=False)
    axs[7].set_title('MI inside of clusters' )
    axs[7].set_xlabel('')
    axs[7].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[8] = plt.subplot2grid((7,4), (4,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(MIi_), "Names": clNames}
    axs[8] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[8], hue="Names", legend=False)
    axs[8].set_title('MI inside (normalized) of clusters' )
    axs[8].set_xlabel('')
    axs[8].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[9] = plt.subplot2grid((7,4), (5,0),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(MIe), "Names": clNames}
    axs[9] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[9], hue="Names", legend=False)
    axs[9].set_title('MI outside of clusters' )
    axs[9].set_xlabel('')
    axs[9].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[10] = plt.subplot2grid((7,4), (5,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(MIe_), "Names": clNames}
    axs[10] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[10], hue="Names", legend=False)
    axs[10].set_title('MI outside (normalized) of clusters' )
    axs[10].set_xlabel('')
    axs[10].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[11] = plt.subplot2grid((7,4), (6,0),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(clAct), "Names": clNames}
    axs[11] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[11], hue="Names", legend=False)
    axs[11].set_title('Number of active site residues in a cluster' )
    axs[11].set_xlabel('')
    axs[11].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    axs[12] = plt.subplot2grid((7,4), (6,2),colspan=2, rowspan=1)
    INTENSITY = {"Intensity": np.array(clAll), "Names": clNames}
    axs[12] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[12], hue="Names", legend=False)
    axs[12].set_title('Number of allosteric site residues in a cluster' )
    axs[12].set_xlabel('')
    axs[12].set_ylabel('')
    plt.xticks(rotation=45, ha='right');

    plt.subplots_adjust(wspace=1, hspace=1)

    try:
        fig.savefig(out_path)
        print("File",out_path+ " created\n")
    except:
        print("Error writing file",out_path,'\n')
