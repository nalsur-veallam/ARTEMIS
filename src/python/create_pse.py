from pymol import cmd
import sys
import json
import numpy as np

print("\nSCRIPT FOR DISPLAYING CLUSTERS ON THE STRUCTURE IS LAUNCHED\n")
colors = ['blue', 'red', 'green', 'cyan', 'hotpink', 'orange', 'olive', 'aquamarine', 'yellow', 'gray', 'purple']

if not ("-f" in sys.argv and "-n" in sys.argv and "-nclust" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f structure.pdb -n name -nclust num_of_clust\n")
    exit()

for i in range(1, 7) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nclust":
        NClusters = int(sys.argv[i+1])
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]

labels_path =  "output/" + name + "/clustering/" + name + "_"+ str(NClusters) + '_clustering'
out_path =  "output/" + name + "/clustering/" + name

try:
    with open(labels_path + '.json') as json_file:
        data = json.load(json_file)

    labels = np.array(data['clustering_labels'])
except:
    print("Error reading file", labels_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f structure.pdb -n name -nclust num_of_clust\n")
    exit()

it= 0
lastResi = None

def coloring(resi, resn, aname, index, b):
    global it, lastResi
    
    resid = resi + resn
    
    if not lastResi:
        lastResi = resid
        
    if lastResi != resid:
        it += 1
        lastResi = resid
        
    clust = labels[it]
    
    clust_name = colors[clust-1] + "_cl_" + str(clust)
    
    text = "resi " + str(resi) + " & resn " + resn  + " & name " + aname + " & index " + str(index)
    cmd.select(clust_name, text, merge=1)
    
    return b

myspace = {'coloring': coloring}

cmd.load(path)
cmd.alter("all",'b = coloring(resi, resn, name, index, b)', space=myspace)
cmd.show_as("cartoon",cmd.get_object_list("all")[0])

for i in range(NClusters):
    obj = colors[i] + "_cl_" + str(i+1)
    cmd.color(colors[i], obj)
    cmd.recolor()

try:
    cmd.save(out_path + "_" + str(NClusters) + "_clusters.pse")
    print("File",out_path + "_"+ str(NClusters) + "_clustering.pse\n")
except:
    print("Error writing file",out_path + "_"+ str(NClusters) + '_clustering.pse\n')
