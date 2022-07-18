import sys
import numpy as np
import pylab as plt
from sklearn.cluster import AgglomerativeClustering
import json


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
map_path = "output/map/" + name + "_map"
out_path = "output/opt_num_of_clust/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
NResidues = data['NResidues']
if not max_of_clust:
    max_of_clust = NResidues

frob = np.sqrt(sum(abs(map_.flatten())**2)) #Frobenius norm of a matrix

distances = np.ones((NResidues, NResidues)) - map_/frob #Residual Distance Matrix

metric_quotient = []
metric_with_log = []

for i in range(min_of_clust, max_of_clust + 1):
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
fig.savefig(out_path + '_metric_quotient.pdf')

fig, ax = plt.subplots(figsize=(15,12), constrained_layout=True)

x = np.linspace(min_of_clust, max_of_clust, max_of_clust - min_of_clust + 1)

plt.plot(x, metric_with_log)
plt.grid()
plt.xlabel('number of clusters (N)', fontsize=20)
plt.ylabel('$MI_{ex}$/log(N)', fontsize=20)
plt.title('$MI_{ex}$/log(N) from number of clusters', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=16)
fig.savefig(out_path + '_metric_with_log.pdf')
