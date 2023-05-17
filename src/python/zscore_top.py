import sys
import numpy as np
import pandas as pd
import json
import itertools as itt
from scipy.stats import zscore
from scipy.special import comb
from tqdm import tqdm

width = .6
max_iters = 10000
noseq = 0

print("\nSCRIPT FOR ALLOSTERIC COMMUNICATION INTENSITY Z-SCORE TOP CALCULATION IS LAUNCHED\n")

if not ("-allsn" in sys.argv and "-f_all" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 11):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -noseq num_of_res(default 0)\n")
    exit()
    
for i in range(1, len(sys.argv)) :
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
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -noseq num_of_res(default 0)\n")
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


try:
    with open(act_path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", act_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", act_path,"by '" + str(as_name) + "' key\n")
    exit()
    
try:
    with open(all_path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", all_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -noseq num_of_res(default 0)\n")
    exit()
    
try:
    all_site = np.array(your_data[alls_name])
except:
    print("Error: Can't get data from file", all_path,"by '" + str(alls_name) + "' key\n")
    exit()

NSites = len(active_site)

cut = NSites
iters = 0
for i in range(NSites):
    if iters >= max_iters:
        cut = i
        print("Warning (zscore_top.py): The size of the active site is too big. Combinations of at most " + str(cut) + " remainders are investigated.\n")
        break
    iters += comb(NSites, i+1)

Combinations = []
for i in range(NSites):
    combinations = itt.combinations(active_site, i+1)
    Combinations.append(combinations)
    

table = pd.DataFrame(columns=['Combination', 'Residues with zscore > 2', 'Zscore > 2 Residues names', 'Residues with zscore > 1.5', 
                              'Zscore > 1.5 Residues names', 'Residues with zscore > 1', 'Zscore > 1 Residues names',  'Mean Intensity', 'STD',
                              'top 1 Residue', 'top 2 Residue', 'top 3 Residue'])

more2_ = []
more1_all = []
columns = ["Residue"]
leng = 0

for combinations in tqdm(Combinations, desc ="Progress", total=cut):
    leng += 1
    if leng > cut:
        break
    
    for active_site in tqdm(combinations, desc =f"ActSite size {leng} progress", total=comb(NSites, leng)):
        intensity = []
        Allosteric = []
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
                    if np.abs(resid - 1 - i) >= noseq:
                        inten += map_[resid - 1][i]
                intensity.append(inten)
        
        sigma = np.std(intensity)
        mean = np.mean(intensity)
        
        intensity = zscore(intensity)
        
        new_names = []
        for i in range(NResidues):
            new_names.append(names[i] + " (" + str(real_numbers[i]) +")")

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity)
        INTENSITY["Residue"] = new_names
        INTENSITY["Allosteric site"] = Allosteric
        INTENSITY = pd.DataFrame(data=INTENSITY)
        INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

        top1 = str(INTENSITY['Residue'][0]) + " [" +  str(np.round(INTENSITY['Intensity'][0], 2)) + "]"
        top2 = str(INTENSITY['Residue'][1]) + " [" +  str(np.round(INTENSITY['Intensity'][1], 2)) + "]"
        top3 = str(INTENSITY['Residue'][2]) + " [" +  str(np.round(INTENSITY['Intensity'][2], 2)) + "]"

        act_s_name = []
        for i in active_site:
            act_s_name.append(names[i-1] + " (" + str(real_numbers[i-1]) +")")
        
        m2R = (INTENSITY[INTENSITY["Intensity"] > 2]["Allosteric site"] == True).sum(axis=0)
        m15R = (INTENSITY[INTENSITY["Intensity"] > 1.5]["Allosteric site"] == True).sum(axis=0)
        m1R = (INTENSITY[INTENSITY["Intensity"] > 1]["Allosteric site"] == True).sum(axis=0)
        INTENSITY1 = INTENSITY[INTENSITY["Intensity"] > 1]
        INTENSITY15 = INTENSITY[INTENSITY["Intensity"] > 1.5]
        INTENSITY2 = INTENSITY[INTENSITY["Intensity"] > 2]
        act_s1 = np.array(INTENSITY1[INTENSITY1["Allosteric site"] == True]["Residue"])
        act_s15 = np.array(INTENSITY15[INTENSITY15["Allosteric site"] == True]["Residue"])
        act_s2 = np.array(INTENSITY2[INTENSITY2["Allosteric site"] == True]["Residue"])

        new_row = pd.Series({'Combination': ", ".join(map(str,act_s_name)), 'Residues with zscore > 2': m2R, 'Zscore > 2 Residues names': ", ".join(map(str,act_s2)), 
                             'Residues with zscore > 1.5': m15R, 'Zscore > 1.5 Residues names': ", ".join(map(str,act_s15)),'Residues with zscore > 1': m1R, 
                             'Zscore > 1 Residues names': ", ".join(map(str,act_s1)), 'Mean Intensity': np.round(mean, 2), 'STD': np.round(sigma, 2),
                             'top 1 Residue': top1, 'top 2 Residue': top2, 'top 3 Residue': top3})
        
        table = pd.concat([table, new_row.to_frame().T], ignore_index=True)
        
        if len(active_site) == 1:
            more2_.append(np.array(INTENSITY2["Residue"]))
            more1_all.append(np.array(INTENSITY1[INTENSITY1["Allosteric site"] == True]["Residue"]))
            columns.append(", ".join(map(str,act_s_name)))
    print('\033[F\033[K\033[F')
            
try:
    table.to_csv(out_path + "_zscore.csv", index=False)
    print("\nFile",out_path + "_zscore.csv created\n")
except:
    print("\nError writing file",out_path + '_zscore.csv\n')

overlap = pd.DataFrame(columns=columns)
overlap_all = pd.DataFrame(columns=columns)
ovAmounts = pd.DataFrame(columns=columns)
ovAmounts_all = pd.DataFrame(columns=columns)
columns = columns[1:]
overlap["Residue"] = columns
overlap_all["Residue"] = columns
ovAmounts["Residue"] = columns
ovAmounts_all["Residue"] = columns

for resid1 in range(NSites):
    lists = []
    lists_all = []
    amounts = []
    amounts_all = []
    for resid2 in range(NSites):
        list_ = list(set(more2_[resid1]) & set(more2_[resid2]))
        list_all = list(set(more1_all[resid1]) & set(more1_all[resid2]))
        lists.append(", ".join(map(str,list_)))
        lists_all.append(", ".join(map(str,list_all)))
        amounts.append(len(np.array(list_)))
        amounts_all.append(len(np.array(list_all)))
        
    overlap[columns[resid1]] = lists
    overlap_all[columns[resid1]] = lists_all
    ovAmounts[columns[resid1]] = amounts
    ovAmounts_all[columns[resid1]] = amounts_all
        
try:
    with pd.ExcelWriter(out_path + "_overlap.xlsx") as writer:
        overlap.to_excel(writer, index=False, sheet_name='overlap (lists)')
        overlap_all.to_excel(writer, index=False, sheet_name='Allosteric sites overlap (lists)')
        ovAmounts.to_excel(writer, index=False, sheet_name='overlap (amount)')
        ovAmounts_all.to_excel(writer, index=False, sheet_name='Allosteric sites overlap (amount)')
    print("\nFile",out_path + "_overlap.xlsx created\n")
except:
    print("\nError writing file",out_path + '_overlap.xlsx\n')
    
    
    
    
    
