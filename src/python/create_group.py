import sys
import json
import numpy as np
from pymol import cmd

print("\nSCRIPT FOR CREATING GROUP FILE FROM PYMOL SELECTION IS LAUNCHED\n")

if not ("-ps" in sys.argv and "-f" in sys.argv and "-sel" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f source.json -ps pymol_session.pse -sel selection_name\n")
    exit()
    
for i in range(1, 7) :
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-ps":
        ses_path = sys.argv[i+1]
    if sys.argv[i] == "-sel":
        selname = sys.argv[i+1]

k = 0
resids = []
group = []
jgroup = []

def write(resi, chain, resn, name, index, b):
    group.append(chain+resi+resn+name+str(index))
    
    return b

def check(resi, chain, resn, name, index, b):
    resid = chain+resi+resn
    
    if not(resid in resids):
        resids.append(resid)
        global k
        k+=1
        
        atom = chain+resi+resn+name+str(index)
        if atom in group:
            jgroup.append(k)
            
    return b

myspace = {'write': write, 'check': check}

cmd.load(ses_path)
cmd.alter(selname,'b = write(resi, chain, resn, name, index, b)', space=myspace)
cmd.alter("all",'b = check(resi, chain, resn, name, index, b)', space=myspace)

try:        
    with open(path) as json_file:
        data = json.load(json_file)
        
    data[selname] = jgroup
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
    print("\nFile",path,"wrote\n")
except:
    data = {}
    data[selname] = jgroup
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
    print("\nFile",path," created\n")

