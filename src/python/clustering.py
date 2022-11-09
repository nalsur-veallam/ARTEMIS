import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import json

heigh = .1

if not ("-n" in sys.argv and "-nclust" in  sys.argv  and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -n name -nclust num_of_clust")
    exit()
for i in range(1, 5) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nclust":
        NClusters = sys.argv[i+1]

if not (name and NClusters):
    print("USAGE:\n"+sys.argv[0]+" -n name -nclust num_of_clust")
    exit()

map_path = "output/map/" + name + "_map"
out_path = "output/clustering/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
NResidues = data['NResidues']
real_numbers = np.array(data['real_numbers'])

frob = np.sqrt(sum(abs(map_.flatten())**2)) #Frobenius norm of a matrix

distances = np.ones((NResidues, NResidues)) - map_/frob #Residual Distance Matrix

clustering = AgglomerativeClustering(n_clusters = int(NClusters)).fit(distances)

MAP = np.zeros((NResidues, NResidues))
for j in range(NResidues):
    for k in range(NResidues):
        MAP[k][j] = np.nan
        if clustering.labels_[k] == clustering.labels_[j]:
            MAP[k][j] = clustering.labels_[k] + 1
            
labels = []
for i in range(NResidues):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
            
CLUSTERING = pd.DataFrame(data=MAP[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(CLUSTERING, annot=False, cmap="tab20b", linecolor='black')#, linewidths = 0.1)

plt.title('Clustering matrix for ' + name + " with " + NClusters + " clusters", fontsize=20)
fig.savefig(out_path + "_"+ NClusters + '_clustering.pdf')


labels = []
for l in clustering.labels_:
    labels.append(int(l)+1)

new_data = {}
new_data['clustering_labels'] = labels
new_data['names'] = data['names']
new_data['NResidues'] = NResidues
new_data['name'] = name
new_data['NClusters'] = NClusters
new_data['real_numbers'] = data['real_numbers']
with open(out_path + "_"+ NClusters + '_clustering.json', 'w') as outfile:
    json.dump(new_data, outfile)
    
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

new_clustering = AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(pd.DataFrame(data=distances, index=names, columns=names))

fig, axs = plt.subplots(figsize=(12, heigh*NResidues), constrained_layout=True)
plt.tick_params(axis='both', which='major', labelsize=20)
labels = []
for i in range(NResidues):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")

plot_dendrogram(new_clustering, truncate_mode="level", orientation="right", labels=labels)

plt.title('Dendrogram for ' + name, fontsize=20)
plt.xlabel('Distance', fontsize=20)
plt.ylabel('Residue', fontsize=20)
fig.savefig(out_path + '_dendrogram.pdf')
