import sys
import numpy as np
import pylab as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
import json
from tqdm import tqdm
import pandas as pd

print("\nSCRIPT FOR SEARCHING THE OPTIMAL NUMBER OF CLUSTERS IS LAUNCHED\n")
heigh = .1

if not (sys.argv[1] == "-n" and len(sys.argv) >= 3):
    print("USAGE:\n"+sys.argv[0]+" -n name -min min_num_of_clust -max max_num_of_clust")
    exit()

min_of_clust = 2
max_of_clust = 0

for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-min":
        min_of_clust = int(sys.argv[i+1])
    if sys.argv[i] == "-max":
        max_of_clust = int(sys.argv[i+1])
    
name = sys.argv[2]
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/clustering/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -min min_num_of_clust -max max_num_of_clust\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
    
try:
    NResidues = int(data['NResidues'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'NResidues' key. Continues without this data.\n")
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

if not max_of_clust:
    max_of_clust = NResidues

frob = np.sqrt(sum(abs(map_.flatten())**2)) #Frobenius norm of a matrix

distances = np.ones((NResidues, NResidues)) - map_/frob #Residual Distance Matrix

metric_quotient = []
metric_with_log = []

for i in tqdm(range(min_of_clust, max_of_clust + 1)):
    internal_mi = 0
    external_mi = 0
    clustering = AgglomerativeClustering(n_clusters = i).fit(distances)
    for j in range(NResidues):
        for k in range(NResidues):
            if clustering.labels_[k] != clustering.labels_[j]:
                external_mi += map_[k][j]
            if clustering.labels_[k] == clustering.labels_[j]:
                internal_mi += map_[k][j]
    metric_quotient.append(internal_mi/external_mi)
    metric_with_log.append(external_mi/np.log(i))

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
try:
    fig.savefig(out_path + '_dendrogram.pdf')
    print("File",out_path + "_dendrogram.pdf created\n")
except:
    print("Error writing file",out_path + '_dendrogram.pdf\n')
