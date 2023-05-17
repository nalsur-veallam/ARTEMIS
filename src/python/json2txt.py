import sys
import numpy as np
import json

print("\nSCRIPT FOR CALCULATING COMPACTNESS IS LAUNCHED\n")

p_charge = ["ARG", "LYS"]
n_charge = ["ASP", "GLU"]

if not ("-n" in sys.argv and "-o" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -n name -o out_path\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-o":
        o_file = sys.argv[i+1]
        
        
map_path =  "output/" + name + "/map/" + name + "_map"

try:
    with open(map_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", map_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
    exit()

try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", map_path + '.json',"by 'map' key\n")
    exit()
     
try:
    names = np.array(data['names'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key\n")
    exit()
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'real_numbers' key\n")
    exit()
    
NResidues = int(np.shape(map_)[0])

def do_line(res1, res2):
    line = []
    line.append(names[res1])
    line.append(real_numbers[res1])
    line.append(names[res2])
    line.append(real_numbers[res2])
    line.append(map_[res1, res2])
    return line

try:
    f = open(o_file + '.txt', 'w')
    for i in range(NResidues-1):
        for j in range(i+1, NResidues):
            line = do_line(i, j)
            f.write(" ".join(map(str,line)) + '\n')
    line = do_line(NResidues-1, NResidues-1)
    f.write(" ".join(map(str,line)))
    f.close()
    print("File",o_file + ".txt created\n")
except:
    print("Error writing file", o_file + ".txt\n")
    exit()
