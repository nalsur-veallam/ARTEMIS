import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json

diag = True
norm = False

def draw_map(args):

    print("\nSCRIPT FOR DRAWING MI MATRIX IS LAUNCHED\n")

    map_path =  args.files[0]
    out_path =  "./map"

    try:
        with open(map_path) as json_file:
            data = json.load(json_file)
    except:
        print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
        exit()

    try:
        map_ = np.array(data['map'])
        NResidues = len(map_)
    except:
        print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
        exit()

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

    labels = []
    for i in range(len(names)):
        labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
        if not diag:
            map_[i,i] = np.nan

    if norm:
        map_ = map_ / np.nanmax(map_)

    MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    sns.heatmap(MIE, annot=False, cmap="icefire")

    plt.title('Mutual information on residues', fontsize=20)
    try:
        fig.savefig(out_path + '.pdf')
        print("File",out_path + ".pdf created\n")
    except:
        print("Error writing file",out_path + '.pdf\n')
