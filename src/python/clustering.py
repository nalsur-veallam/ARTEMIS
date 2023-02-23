import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import json

colors = ['blue', 'red', 'green', 'cyan', 'hotpink', 'orange', 'olive', 'aquamarine', 'yellow', 'gray', 'purple']
print("\nSCRIPT OF AGGLOMERATIVE CLUSTERING IS LAUNCHED\n")

if not ("-n" in sys.argv and "-nclust" in  sys.argv  and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -n name -nclust num_of_clust")
    exit()
for i in range(1, 5) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nclust":
        NClusters = sys.argv[i+1]

map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/clustering/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nclust num_of_clust\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
    
try:
    NResidues = int(data['NResidues'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'NResidues' key. Continues without this data.")
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
    
clColors = []
for i in range(int(NClusters)):
    clColors.append(colors[i])
            
CLUSTERING = pd.DataFrame(data=MAP[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(CLUSTERING, annot=False, cmap=clColors, linecolor='black', cbar=False)#, linewidths = 0.1)

plt.title('Clustering matrix for ' + name + " with " + str(NClusters) + " clusters", fontsize=20)
try:
    fig.savefig(out_path + "_"+ str(NClusters) + '_clustering.pdf')
    print("File",out_path + "_"+ str(NClusters) + "_clustering.pdf")
except:
    print("Error writing file",out_path + "_"+ str(NClusters) + '_clustering.pdf')


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
try:
    with open(out_path + "_"+ str(NClusters) + '_clustering.json', 'w') as outfile:
        json.dump(new_data, outfile)
        print("File",out_path + "_"+ str(NClusters) + "_clustering.json\n")
except:
    print("Error writing file",out_path + "_"+ str(NClusters) + '_clustering.json\n')
