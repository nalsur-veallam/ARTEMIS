import sys
import numpy as np
import pandas as pd
import json
import itertools as itt
from scipy.stats import zscore

width = .6

if not ("-allsn" in sys.argv and "-f_all" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) == 11):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name")
    exit()
    
for i in range(1, 11) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f_act":
        act_path = sys.argv[i+1]
    if sys.argv[i] == "-allsn":
        alls_name = sys.argv[i+1]
    if sys.argv[i] == "-f_all":
        all_path = sys.argv[i+1]
        
if not (name and as_name and act_path and all_path and alls_name):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name")
    exit()
    
    
map_path = "output/map/" + name + "_map"
out_path = "output/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
names = np.array(data['names'])
real_numbers = np.array(data['real_numbers'])
NResidues = data['NResidues']

with open(act_path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

with open(all_path) as json_file:
    your_data = json.load(json_file)

all_site = np.array(your_data[alls_name])


NSites = len(active_site)

Combinations = []
for i in range(NSites):
    combinations = itt.combinations(active_site, i+1)
    Combinations.append(combinations)
    

table = pd.DataFrame(columns=['Combination', 'Residues with zscore > 2', 'Zscore > 2 Residues names', 'Residues with zscore > 1.5', 'Zscore > 1.5 Residues names', 'Residues with zscore > 1', 'Zscore > 1 Residues names', 
                              'Number of Argenins with zscore > 2', 'Mean Intensity', 'STD'])

for combinations in Combinations:
    for active_site in combinations:
        intensity = []
        Allosteric = []
        Argenin = []
        for i in range(NResidues):
            
            if i+1 in all_site:
                Allosteric.append(True)
            else:
                Allosteric.append(False)
            
            if i + 1 in active_site:
                intensity.append(0)
            else:
                inten = 0
                for resid in active_site:
                    inten += map_[resid - 1][i]
                intensity.append(inten)
        
        sigma = np.std(intensity)
        mean = np.mean(intensity)
        
        intensity = zscore(intensity)
        
        new_names = []
        for i in range(NResidues):
            if names[i] == "ARG":
                Argenin.append(True)
            else:
                Argenin.append(False)
            new_names.append(names[i] + " (" + str(real_numbers[i]) +")")

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity)
        INTENSITY["Residue"] = new_names
        INTENSITY["Allosteric site"] = Allosteric
        INTENSITY["Argenin"] = Argenin
        INTENSITY = pd.DataFrame(data=INTENSITY)
        INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

        top = INTENSITY["Intensity"]

        act_s_name = []
        for i in active_site:
            act_s_name.append(names[i-1] + " (" + str(real_numbers[i-1]) +")")
        
        m2R = (INTENSITY[INTENSITY["Intensity"] > 2]["Allosteric site"] == True).sum(axis=0)
        m15R = (INTENSITY[INTENSITY["Intensity"] > 1.5]["Allosteric site"] == True).sum(axis=0)
        m1R = (INTENSITY[INTENSITY["Intensity"] > 1]["Allosteric site"] == True).sum(axis=0)
        ARGs = (INTENSITY[INTENSITY["Intensity"] > 2]["Argenin"] == True).sum(axis=0)
        INTENSITY1 = INTENSITY[INTENSITY["Intensity"] > 1]
        INTENSITY15 = INTENSITY[INTENSITY["Intensity"] > 1.5]
        INTENSITY2 = INTENSITY[INTENSITY["Intensity"] > 2]
        act_s1 = np.array(INTENSITY1[INTENSITY1["Allosteric site"] == True]["Residue"])
        act_s15 = np.array(INTENSITY15[INTENSITY15["Allosteric site"] == True]["Residue"])
        act_s2 = np.array(INTENSITY2[INTENSITY2["Allosteric site"] == True]["Residue"])

        new_row = pd.Series({'Combination': ", ".join(map(str,act_s_name)), 'Residues with zscore > 2': m2R, 'Zscore > 2 Residues names': ", ".join(map(str,act_s2)), 'Residues with zscore > 1.5': m15R, 'Zscore > 1.5 Residues names': ", ".join(map(str,act_s15)),'Residues with zscore > 1': m1R, 'Zscore > 1 Residues names': ", ".join(map(str,act_s1)), 'Number of Argenins with zscore > 2': ARGs, 'Mean Intensity': mean, 'STD': sigma})
        
        table = pd.concat([table, new_row.to_frame().T], ignore_index=True)
        
table.to_csv(out_path + "_zscore.csv")
