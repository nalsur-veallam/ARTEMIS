from pymol import cmd
import sys
import json
import numpy as np

if not ("-f" in sys.argv and "-n" in sys.argv and "-nclust" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f structure.pdb -n name -nclust num_of_clust")
    exit()

for i in range(1, 7) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nclust":
        NClusters = int(sys.argv[i+1])
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
        
if not (name and NClusters and path):
    print("USAGE:\n"+sys.argv[0]+" -f structure.pdb -n name -nclust num_of_clust")
    exit()

labels_path = "output/clustering/" + name + "_"+ str(NClusters) + '_clustering'
out_path = "output/clustering/" + name

with open(labels_path + '.json') as json_file:
    data = json.load(json_file)

labels = np.array(data['clustering_labels'])
names = np.array(data['names'])
NResidues = int(data['NResidues'])
real_numbers = np.array(data['real_numbers'])

groups = []
group_names = []

for i in range(NClusters):
    group = []
    group_name = []
    for k in range(NResidues):
        if labels[k] == i + 1:
            group.append(k+1)
            group_name.append(names[k])
    groups.append(group)
    group_names.append(group_name)

def transl_to_text(group, group_names):
    text = "resi " + '+'.join(str(i) for i in group) + " & resn " + '+'.join(str(i) for i in group_names)
    return text


cmd.load(path)

for i in range(NClusters):
    for j in range(len(groups[i])):
        text = "resi " + str(groups[i][j])  # + " & resn " + str(group_names[i][j])  # transl_to_text(groups[i], group_names[i])
        cmd.select("cluster_" + str(i+1), text, merge=1)

cmd.save(out_path + "_" + str(NClusters) + "_session.pse")
