import sys
import json
import numpy as np
from pymol import cmd

K = 5
K_1 = 4
N = 15

if not ("-strc" in sys.argv and "-grn" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) == 9):
    print("USAGE:\n"+sys.argv[0]+" -f group.json -strc sctructure.pdb(.gro ...) -n name -grn group_name")
    exit()
    
for i in range(1, 9) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-grn":
        gr_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-strc":
        str_path = sys.argv[i+1]
        
if not (name and gr_name and path and str_path):
    print("USAGE:\n"+sys.argv[0]+" -f active_site.json -n name -asn active_site_name")
    exit()
    
    
out_path = "output/analysis/" + name

with open(path) as json_file:
    your_data = json.load(json_file)

group_site = np.array(your_data[gr_name])

group = []

def selection(resi, idx):
    if not (idx in group) and (int(resi) in group_site):
        group.append(idx)
        

myspace = {'selection': selection}

cmd.load(str_path)
cmd.alter("all",'selection(resi, index)', space=myspace)


def do_line(arr):
    line = []
    first = True
    for idx in arr:
        if first:
            k = len(str(idx))
            for i in range(K_1 - k):
                line.append(' ')
            line.append(str(idx))
            first = False
        else:
            k = len(str(idx))
            for i in range(K - k):
                line.append(' ')
            line.append(str(idx))
    return line
    
n = len(group) // N
group = np.array(group)
print(len(group))

f = open(out_path + '.ndx', 'w')
f.write('['+name+']\n')
for i in range(n):
    line = do_line(group[N*i : N*(i+1)])
    f.write("".join(map(str,line)) + '\n')
line = do_line(group[(n)*N :-1])
f.write("".join(map(str,line)))
f.close()
    
    



    
