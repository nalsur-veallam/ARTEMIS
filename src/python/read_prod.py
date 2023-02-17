import numpy as np
import json
import sys

prot_prot = False

if not ("-prod" in sys.argv and "-f" in sys.argv and "-sn" in sys.argv and "-cutoff" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+"-prod prodigy_file -f source_json -sn site_name -cutoff cutoff(Angstrom) -pp")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-cutoff":
        cutoff = float(sys.argv[i+1])
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
    if sys.argv[i] == "-prod":
        prod_path = sys.argv[i+1]
    if sys.argv[i] == "-sn":
        sname = sys.argv[i+1]
    if sys.argv[i] == "-pp":
        prot_prot = True

data = np.loadtxt(prod_path, dtype=np.unicode_)

if not prot_prot:
    site = data[np.argwhere(data[:, 8].astype('float64') < cutoff).reshape(-1)]
else:
    site = data[np.argwhere(data[:, 8].astype('float64') < cutoff).reshape(-1)]

try:        
    with open(path) as json_file:
        data = json.load(json_file)
        
    data[sname] = np.unique((site[:, 2]).astype('int32') - 13).tolist()
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
except:
    data = {}
    data[sname] = np.unique((site[:, 2]).astype('int32') - 13).tolist()
    with open(path, 'w') as outfile:
        json.dump(data, outfile)
