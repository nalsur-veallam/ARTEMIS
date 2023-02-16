import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
import json

heigh = .1
noseq = 0
colors = ['blue', 'red', 'green', 'cyan', 'hotpink', 'orange', 'olive', 'aquamarine', 'yellow', 'gray', 'purple']

if not ("-allsn" in sys.argv and "-f_all" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and "-nclust" in  sys.argv  and len(sys.argv) >= 13):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -asn active_site_name 
          -f_all allosteric_site.json -allsn allosteric_site_name -n name -nclust num_of_clust -noseq num_of_res(default 0)")
    exit()
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nclust":
        NClusters = int(sys.argv[i+1])
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f_act":
        act_path = sys.argv[i+1]
    if sys.argv[i] == "-allsn":
        alls_name = sys.argv[i+1]
    if sys.argv[i] == "-f_all":
        all_path = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
        
if not (name and as_name and act_path and all_path and alls_name and NClusters):
    print("USAGE:\n"+sys.argv[0]+"USAGE:\n"+sys.argv[0]+" -f_act active_site.json -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -n name -nclust num_of_clust")
    exit()

map_path =  "output/" + name + "/map/" + name + "_map"
labels_path =  "output/" + name + "/clustering/" + name + "_" + str(NClusters) + "_clustering"
out_path =  "output/" + name + "/clustering/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
NResidues = data['NResidues']

with open(labels_path + '.json') as json_file:
    data = json.load(json_file)
    
labels = np.array(data['clustering_labels'])

with open(act_path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

with open(all_path) as json_file:
    your_data = json.load(json_file)

all_site = np.array(your_data[alls_name])

MIi = np.zeros(NClusters)
MIi_ = np.zeros(NClusters)
MIe = np.zeros(NClusters)
MIe_ = np.zeros(NClusters)
H = np.zeros(NClusters)
MIE = np.zeros((NClusters, NClusters))
intensity_act = np.zeros(NClusters)
intensity_all = np.zeros(NClusters)
intensity_act_ = np.zeros(NClusters)
intensity_all_ = np.zeros(NClusters)
clAct = np.zeros(NClusters)
clAll = np.zeros(NClusters)
clNames = []
clColors = []
quantity = np.zeros(NClusters)

for cl in range(NClusters):
    quantity[cl] = np.sum(labels == cl+1)
            
    clNames.append(colors[cl] + " cl(" + str(cl+1)+")")
    if colors[cl] == 'aquamarine':
        clNames[cl] = "aquam cl(8)"
    clColors.append(colors[cl])
    
    for res in active_site:
        if labels[res-1] == cl+1:
            clAct[cl] += 1
    for res in all_site:
         if labels[res-1] == cl+1:
            clAll[cl] += 1

for i in range(NResidues):
    H[labels[i]-1]+= map_[i,i]
    for j in range(i, NResidues):
        if labels[i] == labels[j]:
            MIi[labels[i]-1]+= map_[i,j]
            MIi_[labels[i]-1]+= map_[i,j]/quantity[labels[i]-1]
        else:
            MIe[labels[i]-1]+= map_[i,j]
            MIe[labels[j]-1]+= map_[i,j]
            MIe_[labels[i]-1]+= map_[i,j]/quantity[labels[i]-1]
            MIe_[labels[j]-1]+= map_[i,j]/quantity[labels[j]-1]
            MIE[labels[i]-1, labels[j]-1] += map_[i,j]
            MIE[labels[j]-1, labels[i]-1] += map_[i,j]
            
        if j in active_site and np.abs(j - i) >= noseq:
            intensity_act[labels[i]-1] += map_[i,j]
            intensity_act_[labels[i]-1] += map_[i,j]/quantity[labels[i]-1]
            
        if j in all_site and np.abs(j - i) >= noseq:
            intensity_all[labels[i]-1] += map_[i,j]
            intensity_all_[labels[i]-1] += map_[i,j]/quantity[labels[i]-1]

clColors = sns.color_palette(clColors, as_cmap=True)

if NClusters > 7:
    fig = plt.figure(figsize=(24, 36))
else:
    fig = plt.figure(figsize=(16, 24))

axs = [None for _ in range(13)]

axs[0] = plt.subplot2grid((7,4), (0,0),colspan=2, rowspan=2)
CLUSTERING = pd.DataFrame(data=MIE[::-1, :], index=np.array(clNames)[::-1], columns=np.array(clNames))
axs[0] = sns.heatmap(CLUSTERING, cmap="OrRd", annot=True, ax=axs[0], cbar=False)
axs[0].set_title('MIE matrix for ' + name + " with " +str(NClusters) + " clusters")

if NClusters > 7:
    font = {'size'   : 20}
    plt.rc('font', **font)

axs[1] = plt.subplot2grid((7,4), (0,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(quantity), "Names": clNames}
axs[1] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[1])
axs[1].set_title('Quantity of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[2] = plt.subplot2grid((7,4), (1,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(H), "Names": clNames}
axs[2] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[2])
axs[2].set_title('Informational entropy of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[3] = plt.subplot2grid((7,4), (2,0),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(intensity_act) / np.sqrt(sum(abs(np.array(intensity_act).flatten())**2)), "Names": clNames}
axs[3] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette="OrRd", dodge=False, ax=axs[3], hue="Intensity")
axs[3].legend_.remove()
axs[3].set_title('Intensity of connectivity of clusters\n with the active site for ' + name)
plt.xticks(rotation=45, ha='right');

axs[4] = plt.subplot2grid((7,4), (2,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(intensity_act_) / np.sqrt(sum(abs(np.array(intensity_act_).flatten())**2)), "Names": clNames}
axs[4] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette="OrRd", dodge=False, ax=axs[4], hue="Intensity")
axs[4].legend_.remove()
axs[4].set_title('Intensity of connectivity (normalized) of clusters\n with the active site for ' + name)
plt.xticks(rotation=45, ha='right');

axs[5] = plt.subplot2grid((7,4), (3,0),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(intensity_all) / np.sqrt(sum(abs(np.array(intensity_all).flatten())**2)), "Names": clNames}
axs[5] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette="OrRd", dodge=False, ax=axs[5], hue="Intensity")
axs[5].legend_.remove()
axs[5].set_title('Intensity of connectivity of clusters\n with the allosteric site for ' + name)
plt.xticks(rotation=45, ha='right');

axs[6] = plt.subplot2grid((7,4), (3,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(intensity_all_) / np.sqrt(sum(abs(np.array(intensity_all_).flatten())**2)), "Names": clNames}
axs[6] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette="OrRd", dodge=False, ax=axs[6], hue="Intensity")
axs[6].legend_.remove()
axs[6].set_title('Intensity of connectivity (normalized) of clusters\n with the allosteric site for ' + name)
plt.xticks(rotation=45, ha='right');

axs[7] = plt.subplot2grid((7,4), (4,0),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(MIi), "Names": clNames}
axs[7] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[7])
axs[7].set_title('Mi inside of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[8] = plt.subplot2grid((7,4), (4,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(MIi_), "Names": clNames}
axs[8] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[8])
axs[8].set_title('MI inside (normalized) of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[9] = plt.subplot2grid((7,4), (5,0),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(MIe), "Names": clNames}
axs[9] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[9])
axs[9].set_title('Mi outside of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[10] = plt.subplot2grid((7,4), (5,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(MIe_), "Names": clNames}
axs[10] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[10])
axs[10].set_title('MI outside (normalized) of clusters ' + name)
plt.xticks(rotation=45, ha='right');

axs[11] = plt.subplot2grid((7,4), (6,0),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(clAct), "Names": clNames}
axs[11] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[11])
axs[11].set_title('Number of active site residues in a cluster ' + name)
plt.xticks(rotation=45, ha='right');

axs[12] = plt.subplot2grid((7,4), (6,2),colspan=2, rowspan=1)
INTENSITY = {"Intensity": np.array(clAll), "Names": clNames}
axs[12] = sns.barplot(x="Names", y="Intensity", data=INTENSITY, palette=clColors, dodge=False, ax=axs[12])
axs[12].set_title('Number of allosteric site residues in a cluster ' + name)
plt.xticks(rotation=45, ha='right');

plt.subplots_adjust(wspace=1, hspace=1)

fig.savefig(out_path + "_"+ str(NClusters) + '_analysis.pdf')
       
    
    
