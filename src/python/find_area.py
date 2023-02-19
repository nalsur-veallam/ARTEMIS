import sys
import json
import numpy as np
from pymol import cmd

tchain = False

if not ("-f" in sys.argv and "-sn" in sys.argv and "-cutoff" in sys.argv and "-strc" in sys.argv and "-chain" in sys.argv 
        and ("-ligname" in sys.argv or "-chain2" in sys.argv) and len(sys.argv) == 13):
    print("USAGE:\n"+sys.argv[0]+"-f source.json -sn site_name -chain protein_chain_name -ligname ligand_name -strc sctructure.pdb(.gro ...) -cutoff cutoff(Angstrom)")
    exit()
    
for i in range(1, 13) :
    if sys.argv[i] == "-cutoff":
        cutoff = float(sys.argv[i+1])
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-sn":
        sname = sys.argv[i+1]
    if sys.argv[i] == "-chain":
        sel1 = "bm. c. " + sys.argv[i+1] + " and not resn HOH"
    if sys.argv[i] == "-chain2":
        sel2 = "bm. c. " + sys.argv[i+1] + " and not resn HOH"
        tchain = True
    if sys.argv[i] == "-ligname":
        sel2 = "resn " + sys.argv[i+1]
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]

cmd.load(str_path)

m1 = cmd.get_model(sel2+" around "+str(cutoff)+" and "+sel1+" and not name c+o+n")
prot = cmd.get_model(sel1)
m2 = cmd.get_model(sel2)

resids = []
resid = ""
for i in range(len(m1.atom)):
    if m1.atom[i].chain+m1.atom[i].resn+m1.atom[i].resi != resid:
        resid = m1.atom[i].chain+m1.atom[i].resn+m1.atom[i].resi
        resids.append(resid)
area = []
resid = ""
for i in range(len(prot.atom)):
    if prot.atom[i].chain+prot.atom[i].resn+prot.atom[i].resi != resid:
        resid = prot.atom[i].chain+prot.atom[i].resn+prot.atom[i].resi
        if resid in resids and not tchain:
            area.append(True)
        else:
            area.append(False)

if tchain:
    resid = ""
    counter = -1
    for c1 in range(len(prot.atom)):
        if prot.atom[c1].chain+prot.atom[c1].resn+prot.atom[c1].resi != resid:
            resid = prot.atom[c1].chain+prot.atom[c1].resn+prot.atom[c1].resi
            counter += 1
            
        if not area[counter]:
            for c2 in range(len(m2.atom)):
                distance=np.sqrt(sum(map(lambda f: (f[0]-f[1])**2, zip(prot.atom[c1].coord,m2.atom[c2].coord))))
                if distance<float(cutoff):
                    area[counter] = True

site = []
for i in range(len(area)):
    if area[i]:
        site.append(i+1)
try:        
    with open(path) as json_file:
        data = json.load(json_file)
        
    data[sname] = site
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
except:
    data = {}
    data[sname] = site
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
    

