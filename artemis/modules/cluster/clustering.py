import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
import json

colors = ['blue', 'green', 'red', 'cyan', 'yellow', 'magenta', 'salmon', 'lime', 'hotpink', 'orange', 'chartreuse', 'limegreen', 'olive', 'purple', 'teal', 'gray', 'pink', 'firebrick', 'chocolate', 'brown', 'wheat', 'violet', 'aquamarine', 'palegreen', 'lightblue', 'lightpink', 'skyblue', 'silver', 'gold']

def do_clustering(Clusters, out_path):

    distance_threshold = None

    print("\nSCRIPT OF AGGLOMERATIVE CLUSTERING IS LAUNCHED\n")

    if Clusters.reference_group is not None:

        map_ = Clusters.map_[:, Clusters.reference_group-1] #Residual Distance Matrix

    else:

        map_ = Clusters.map_ #Residual Distance Matrix=

    if Clusters.restriction_group is not None:

        dists = []
        for i in Clusters.restriction_group:
            for j in Clusters.restriction_group:
                dists.append(np.linalg.norm(map_[i-1] - map_[j-1]))

        distance_threshold = max(dists)+1e-15

    if Clusters.restriction_group is None:
        clustering = AgglomerativeClustering(n_clusters = int(Clusters.NClusters)).fit(map_)
        Clusters.NClusters = clustering.n_clusters_
    else:
        clustering = AgglomerativeClustering(distance_threshold=distance_threshold, n_clusters=None, linkage='complete', compute_full_tree=True).fit(map_)
        Clusters.NClusters = clustering.n_clusters_

    unique, counts = np.unique(clustering.labels_, return_counts=True)

    idxs = np.argsort(-1*counts)
    unique = unique[idxs]
    Clusters.clustering_labels = []

    for l in clustering.labels_:
        real_l = int(np.argwhere(unique==l))+1
        Clusters.clustering_labels.append(int(real_l))

    Clusters.clustering_labels = np.array(Clusters.clustering_labels)

    MAP = np.zeros((Clusters.NResidues, Clusters.NResidues))
    for j in range(Clusters.NResidues):
        for k in range(Clusters.NResidues):
            MAP[k][j] = np.nan
            if Clusters.clustering_labels[k] == Clusters.clustering_labels[j]:
                MAP[k][j] = Clusters.clustering_labels[k] + 1

    labels = []
    for i in range(Clusters.NResidues):
        labels.append(Clusters.names[i] + " (" + str(Clusters.real_numbers[i]) + ")")

    clColors = []
    for i in range(int(Clusters.NClusters)):
        clColors.append(colors[i])

    CLUSTERING = pd.DataFrame(data=MAP[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    sns.heatmap(CLUSTERING, annot=False, cmap=clColors, linecolor='black', cbar=False)#, linewidths = 0.1)

    plt.title('Clustering matrix with ' + str(Clusters.NClusters) + " clusters", fontsize=20)
    try:
        fig.savefig(out_path)
        print("File",out_path, "created")
    except:
        print("Error writing file", out_path, "\n")


    if out_path == "clustering.pdf":
        fname = "clustering.json"
    else:
        fname = ''.join(out_path.split('.')[:-1]) + '.json'
    Clusters.write(fname)


def do_spectral_clustering(Clusters, out_path):

    Clusters.restriction_group = None
    Clusters.reference_group = None

    print("\nSCRIPT OF SPECTRAL CLUSTERING IS LAUNCHED\n")

    clustering = SpectralClustering(n_clusters = int(Clusters.NClusters), affinity='precomputed').fit(Clusters.map_)
    Clusters.NClusters = int(Clusters.NClusters)

    unique, counts = np.unique(clustering.labels_, return_counts=True)

    idxs = np.argsort(-1*counts)
    unique = unique[idxs]
    Clusters.clustering_labels = []

    for l in clustering.labels_:
        real_l = int(np.argwhere(unique==l))+1
        Clusters.clustering_labels.append(int(real_l))

    Clusters.clustering_labels = np.array(Clusters.clustering_labels)

    MAP = np.zeros((Clusters.NResidues, Clusters.NResidues))
    for j in range(Clusters.NResidues):
        for k in range(Clusters.NResidues):
            MAP[k][j] = np.nan
            if Clusters.clustering_labels[k] == Clusters.clustering_labels[j]:
                MAP[k][j] = Clusters.clustering_labels[k] + 1

    labels = []
    for i in range(Clusters.NResidues):
        labels.append(Clusters.names[i] + " (" + str(Clusters.real_numbers[i]) + ")")

    clColors = []
    for i in range(int(Clusters.NClusters)):
        clColors.append(colors[i])

    CLUSTERING = pd.DataFrame(data=MAP[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    sns.heatmap(CLUSTERING, annot=False, cmap=clColors, linecolor='black', cbar=False)#, linewidths = 0.1)

    plt.title('Clustering matrix with ' + str(Clusters.NClusters) + " clusters", fontsize=20)
    try:
        fig.savefig(out_path)
        print("File",out_path, "created")
    except:
        print("Error writing file", out_path, "\n")


    if out_path == "clustering.pdf":
        fname = "clustering.json"
    else:
        fname = ''.join(out_path.split('.')[:-1]) + '.json'
    Clusters.write(fname)
