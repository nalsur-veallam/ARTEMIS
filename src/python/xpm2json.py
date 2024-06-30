import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
import re
import json
import sys

print("\nXPM GROMACS MATRIX CONVERTING SCRIPT IS LAUNCHED\n")

if not ("-f" in sys.argv and "-o" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -f xpm_file -o out_path\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f":
        i_file = sys.argv[i+1]
    if sys.argv[i] == "-o":
        o_file = sys.argv[i+1]

symbols = []
counts = []

try:
    f = open(i_file, 'r')
    r_code = True
    r_mat = False
    for row in f:
        if r_code:
            match = re.search(r'\w\s+\w\s#\w{6}\s"', row)
            if match:
                symbol = match[0][0]
                match = re.search(r'"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"', row)
                count = float(re.search(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?', match[0])[0])
                symbols.append(symbol)
                counts.append(count)
            else:
                match = re.search(r'"\d+\s\d+\s+\d+\s\d"', row)
                if match:
                    match = re.search(r'"\d+', match[0])
                    N = int(re.search(r'\d+', match[0])[0])
                    Mat = np.empty(shape=(N,N))
                    
            
        if r_code:
            if re.search(r'/* x-axis', row):
                r_code = False
                r_mat = True
                nrow = 0
                symbols = np.array(symbols)
                counts = np.array(counts)
                
        if r_mat:
            if row[0] == '"':
                n = len(row)
                for i in range(0, N):
                    count = counts[symbols == row[i+1]]
                    Mat[nrow, i] = float(count)
                nrow += 1
    f.close()
except:
    print("Error reading file", i_file, ". USAGE:\n"+sys.argv[0]+" -f xpm_file -o out_path\n")
    exit()

DM = pd.DataFrame(data=(-1*Mat + np.max(Mat))/np.nanmax(-1*Mat + np.max(Mat)), index=np.arange(1, Mat.shape[0] + 1)[::-1], columns=np.arange(1, Mat.shape[0] + 1))


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(DM, annot=False, cmap="icefire")

plt.title('Distance matrix', fontsize=20)
try:
    fig.savefig(o_file + '.pdf')
    print("File",o_file+".pdf created")
except:
    print("Error writing file", o_file+'.pdf')

new_data = {}
new_data['map'] = Mat.tolist()
try:
    with open(o_file+'.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",o_file+".json created\n")    
except:
    print("Error writing file", o_file+'.json\n')
    exit()

