import sys
import json
import numpy as np
from pymol import cmd
import pandas as pd

visual = "Y"
act = True
alls = True
noseq = 0

def max10(array):
    size = len(array)
    top10 = int(0.1*size)
    
    supp = np.zeros(size)
    for i in range(top10):
        supp[i] = 1
    
    d = {'data':array, 'index':np.arange(0,NResidues)}
    
    df = pd.DataFrame(data=d)
    df = df.sort_values(by=['data'], ascending=False)
    
    for i in range(NResidues):
        if np.array(df['data'])[i] < 0:
            supp[i] = 1
    
    df['data'] = df['data']*supp
    df = df.sort_values(by=['index'])
    return np.array(df['data'])

if not ("-strc" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 9):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -noseq num_of_res(default 0)")
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
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
    if sys.argv[i] == "-noseq":
        noseq = int(sys.argv[i+1])
        
if not (name and as_name and act_path and str_path):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...)")
    exit()    
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
NResidues = data['NResidues']

with open(act_path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

intensity = np.zeros(NResidues)

for resid in active_site:
    inten = []
    for i in range(NResidues):
        if i + 1 in active_site:
            inten.append(-1)
        elif np.abs(resid - 1 - i) >= noseq:
            inten.append(map_[resid - 1][i])
    intensity += max10(inten)

if ("alls_name" not in locals() and "all_path" not in locals()):
    print("ATTENTION, you did not specify an allosteric site. The program will continue without it.")
    alls = False
else:
    with open(all_path) as json_file:
        your_data = json.load(json_file)

    all_site = np.array(your_data[alls_name])
    alls_arr = []
    group_alls = []
    group_alls_names = []
    
    for i in range(NResidues):
        if i + 1 in all_site:
            alls_arr.append(1)
        else:
            alls_arr.append(0)

k = 0
resids = []
group = []
group_names = []

def renumbering(resi, resn, b):
    resid = resi + resn
    if not(resid in resids):
        resids.append(resid)
        global k
        k+=1
        
        if intensity[k-1] < 0:
            group.append(resi)
            group_names.append(resn)
        
        if alls:
            if alls_arr[k-1] == 1:
                group_alls.append(resi)
                group_alls_names.append(resn)
        
        return intensity[k-1]
    return intensity[k-1]

myspace = {'renumbering': renumbering}

cmd.load(str_path)
cmd.alter("all",'b = renumbering(resi, resn, b)', space=myspace)

if visual=="Y":
    obj = cmd.get_object_list("all")[0]
    cmd.show_as("cartoon",obj)
    cmd.cartoon("putty", obj)
    cmd.set("cartoon_putty_scale_min", min(intensity),obj)
    cmd.set("cartoon_putty_scale_max", max(intensity),obj)
    cmd.set("cartoon_putty_transform", 0,obj)
    cmd.set("cartoon_putty_radius", 0.4,obj)
    cmd.spectrum("b","rainbow", obj)
    cmd.ramp_new("count", obj, [min(intensity), max(intensity)], "rainbow")
    cmd.recolor()

for j in range(len(group)):
    text = "resi " + str(group[j]) + "& resn " + str(group_names[j])
    cmd.select("active_site", text, merge=1)
        
if act:
    obj = "active_site"
    cmd.show_as("spheres",obj)
    cmd.color("white", obj)
    cmd.recolor()
    
if alls:
    for j in range(len(group_alls)):
        text = "resi " + str(group_alls[j]) + "& resn " + str(group_alls_names[j])
        cmd.select("allosteric_site", text, merge=1)

cmd.save(out_path + "_intensity_top10.pse")
