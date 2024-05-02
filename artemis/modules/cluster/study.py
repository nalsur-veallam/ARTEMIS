import numpy as np
import pylab as plt
from sklearn.cluster import AgglomerativeClustering, SpectralClustering
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from sklearn.metrics import silhouette_score
import json
from tqdm import tqdm
import pandas as pd

def centeroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    return sum_x/length, sum_y/length

def study_clustering(Clusters, min_of_clust, max_of_clust, out_path, spectral):
    print("\nSCRIPT FOR SEARCHING THE OPTIMAL NUMBER OF CLUSTERS IS LAUNCHED\n")
    heigh = .1

    out_path = ''.join(out_path.split('.')[:-1])

    if Clusters.reference_group is not None and not spectral:

        map_ = Clusters.map_[:, Clusters.reference_group-1] #Residual Distance Matrix

    else:

        map_ = Clusters.map_ #Residual Distance Matrix=

    metric_quotient = []
    metric_with_log = []
    metric_inertia = []

    for i in tqdm(range(min_of_clust, max_of_clust + 1)):
        internal_mi = 0
        external_mi = 0
        if spectral:
            clustering = SpectralClustering(n_clusters = i, affinity='precomputed').fit(map_)
        else:
            clustering = AgglomerativeClustering(n_clusters = i).fit(map_)

        for j in range(Clusters.NResidues):
            for k in range(Clusters.NResidues):
                if clustering.labels_[k] != clustering.labels_[j]:
                    external_mi += Clusters.map_[k][j]
                if clustering.labels_[k] == clustering.labels_[j]:
                    internal_mi += Clusters.map_[k][j]
        metric_quotient.append(internal_mi/external_mi)
        metric_with_log.append(external_mi/np.log(i))

        inertia = 0

        for j in range(i):
            centroid = np.zeros(map_.shape[1])
            L = len(np.argwhere(clustering.labels_ == j))
            for k, l in enumerate(clustering.labels_):
                if l == j:
                    centroid += map_[k]/L

            for k, l in enumerate(clustering.labels_):
                if l == j:
                    inertia += np.linalg.norm(map_[k] - centroid)**2

        metric_inertia.append(inertia)


    fig, ax = plt.subplots(figsize=(15,12), constrained_layout=True)

    x = np.linspace(min_of_clust, max_of_clust, max_of_clust - min_of_clust + 1)

    plt.plot(x, metric_quotient)
    plt.grid()
    plt.xlabel('number of clusters', fontsize=20)
    plt.ylabel('$MI_{ex}$/$MI_{in}$', fontsize=20)
    plt.title('$MI_{ex}$/$MI_{in}$ from number of clusters', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    try:
        fig.savefig(out_path + '_metric_quotient.pdf')
        print("\nFile",out_path + "_metric_quotient.pdf created")
    except:
        print("\nError writing file",out_path + '_metric_quotient.pdf')

    fig, ax = plt.subplots(figsize=(15,12), constrained_layout=True)

    x = np.linspace(min_of_clust, max_of_clust, max_of_clust - min_of_clust + 1)

    plt.plot(x, metric_with_log)
    plt.grid()
    plt.xlabel('number of clusters (N)', fontsize=20)
    plt.ylabel('$MI_{ex}$/log(N)', fontsize=20)
    plt.title('$MI_{ex}$/log(N) from number of clusters', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    try:
        fig.savefig(out_path + '_metric_with_log.pdf')
        print("File",out_path + "_metric_with_log.pdf created")
    except:
        print("Error writing file",out_path + '_metric_with_log.pdf')

    fig, ax = plt.subplots(figsize=(15,12), constrained_layout=True)

    plt.plot(x, metric_inertia)
    plt.grid()
    plt.xlabel('number of clusters', fontsize=20)
    plt.ylabel('Inertia', fontsize=20)
    plt.title('K-means inertia from number of clusters', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    try:
        fig.savefig(out_path + '_metric_inertia.pdf')
        print("\nFile",out_path + "_metric_inertia.pdf created")
    except:
        print("\nError writing file",out_path + '_metric_inertia.pdf')

    # Plotting dendrogram

    def plot_dendrogram(model, **kwargs):
        # Create linkage matrix and then plot the dendrogram

        # create the counts of samples under each node
        counts = np.zeros(model.children_.shape[0])
        n_samples = len(model.labels_)
        for i, merge in enumerate(model.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  # leaf node
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count

        linkage_matrix = np.column_stack(
            [model.children_, model.distances_, counts]
        ).astype(float)

        # Plot the corresponding dendrogram
        dendrogram(linkage_matrix, **kwargs)

    if not spectral:
        new_clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(pd.DataFrame(data=map_, index=Clusters.names))

        fig, axs = plt.subplots(figsize=(12, heigh*Clusters.NResidues), constrained_layout=True)
        plt.tick_params(axis='both', which='major', labelsize=20)

        labels = []
        for i in range(Clusters.NResidues):
            labels.append(Clusters.names[i] + " (" + str(Clusters.real_numbers[i]) + ")")

        plot_dendrogram(new_clustering, truncate_mode="level", orientation="right", labels=labels)

        plt.title('Dendrogram', fontsize=20)
        plt.xlabel('Distance', fontsize=20)
        plt.ylabel('Residue', fontsize=20)
        try:
            fig.savefig(out_path + '_dendrogram.pdf')
            print("File",out_path + "_dendrogram.pdf created\n")
        except:
            print("Error writing file",out_path + '_dendrogram.pdf\n')
