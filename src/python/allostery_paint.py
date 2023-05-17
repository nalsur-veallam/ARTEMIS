import sys
import json
import numpy as np
from pymol import cmd

visual = "Y"
act = True
alls = True
filt = False
sasa_filt = False
noseq = 0

print("\nSCRIPT FOR DISPLAYING ALLOSTERIC COMMUNICATION INTENSITY ON THE STRUCTURE IS LAUNCHED\n")

if not ("-strc" in sys.argv and "-asn" in sys.argv and "-f_act" in sys.argv and "-n" in sys.argv and len(sys.argv) >= 9):
    print("USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
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
    if sys.argv[i] == "-filt":
        filt = True
    if sys.argv[i] == "-sasa_filt":
        sasa_filt = True
    
if sasa_filt:
    map_path =  "output/" + name + "/map/" + name + "_sasa_filt_map"
    out_path =  "output/" + name + "/analysis/" + name + "_sasa_filt"
elif filt:
    map_path =  "output/" + name + "/map/" + name + "_filt_map"
    out_path =  "output/" + name + "/analysis/" + name + "_filt"
else:
    map_path =  "output/" + name + "/map/" + name + "_map"
    out_path =  "output/" + name + "/analysis/" + name

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
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
    with open(act_path) as json_file:
        your_data = json.load(json_file)
except:
    print("Error reading file", act_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
    exit()
    
try:
    active_site = np.array(your_data[as_name])
except:
    print("Error: Can't get data from file", act_path,"by '" + str(as_name) + "' key\n")
    exit()

intensity = []
for i in range(NResidues):
    if i + 1 in active_site:
        intensity.append(-1)
    else:
        inten = 0
        for resid in active_site:
            if np.abs(resid - 1 - i) >= noseq:
                inten += map_[resid - 1][i]
        intensity.append(inten)

if ("alls_name" not in locals() and "all_path" not in locals()):
    print("ATTENTION, you did not specify an allosteric site. The program will continue without it.")
    alls = False
else:
    try:
        with open(all_path) as json_file:
            your_data = json.load(json_file)
    except:
        print("Error reading file", all_path, ". USAGE:\n"+sys.argv[0]+" -f_act active_site.json -n name -asn active_site_name -f_all allosteric_site.json -allsn allosteric_site_name -strc sctructure.pdb(.gro ...) -filt -sasa_filt -noseq num_of_res(default 0)\n")
        exit()
        
    try:
        all_site = np.array(your_data[alls_name])
    except:
        print("Error: Can't get data from file", all_path,"by '" + str(alls_name) + "' key\n")
        exit()
        
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

try:
    cmd.save(out_path + "_intensity.pse")
    print("File",out_path + "_intensity.pse created\n")
except:
    print("Error writing file",out_path + '_intensity.pse\n')

