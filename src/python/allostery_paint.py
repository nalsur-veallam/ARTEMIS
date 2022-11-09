import sys
import json
import numpy as np
from pymol import cmd

visual = "Y"
act = True

if not ("-strc" in sys.argv and "-asn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) == 9):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -strc sctructure.pdb(.gro ...) -n name -asn active_site_name")
    exit()
    
for i in range(1, 9) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-asn":
        as_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
        
if not (name and as_name and path and str_path):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name")
    exit()
    
    
map_path = "output/map/" + name + "_map"
out_path = "output/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
NResidues = data['NResidues']

with open(path) as json_file:
    your_data = json.load(json_file)

active_site = np.array(your_data[as_name])

intensity = []
for i in range(NResidues):
        if i + 1 in active_site:
            intensity.append(0)
        else:
            inten = 0
            for resid in active_site:
                inten += map_[resid - 1][i]
            intensity.append(inten)


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
        
        if intensity[k-1] == 0:
            group.append(resi)
            group_names.append(resn)
        
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

cmd.save(out_path + "_intensity.pse")

