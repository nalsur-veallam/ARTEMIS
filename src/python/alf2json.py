import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

print("\nALPHA FOLD MATRIX CONVERTING SCRIPT IS LAUNCHED\n")

if not ("-f" in sys.argv and "-o" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -f alpha_json_file -o out_path\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f":
        i_file = sys.argv[i+1]
    if sys.argv[i] == "-o":
        o_file = sys.argv[i+1]

try:
    with open(i_file) as json_file:
        data = json.load(json_file)
    map_ = np.array(data['pae'])
except:
    print("Error reading file", i_file, ". USAGE:\n"+sys.argv[0]+" -f alpha_json_file -o out_path\n")
    exit()

labels = []
for i in range(len(map_)):
    labels.append(str(i))
    map_[i,i] = np.nan

MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="icefire")

try:
    plt.title('AlphaFold map', fontsize=20)
    fig.savefig(o_file+'.pdf')
    print("File",o_file+".pdf created")
except:
    print("Error writing file", o_file+'.pdf\n')

new_data = {}
new_data['map'] = map_.tolist()

try:
    with open(o_file+'.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",o_file+".json created\n")    
except:
    print("Error writing file", o_file+'.json\n')
    exit()
