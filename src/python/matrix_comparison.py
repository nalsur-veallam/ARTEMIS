import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

diag = True
print("\nMATRIX COMPARISON SCRIPT IS LAUNCHED\n")

if not ("-f1" in sys.argv and "-f2" in sys.argv and "-o" in sys.argv and len(sys.argv) >= 7):
    print("USAGE:\n"+sys.argv[0]+" -f1 matrix1.json -f2 matrix2.json -o output -nodiag\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f1":
        f1 = sys.argv[i+1]
    if sys.argv[i] == "-o":
        out_path = sys.argv[i+1]
    if sys.argv[i] == "-f2":
        f2 = sys.argv[i+1]
    if sys.argv[i] == "-nodiag":
        diag = False

try:
    with open(f1) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f1, ". USAGE:\n"+sys.argv[0]+" -f1 matrix1.json -f2 matrix2.json -o output -nodiag\n")
    exit()
    
try:
    map_ = np.array(data['map'])
except:
    print("Error: Can't get data from file", f1 ,"by 'map' key\n")
    exit()
    
map_[np.isnan(map_)] = 0

names = []
for i in range(len(map_)):
    names.append(str(i))

if not diag:
    map_ = map_ - np.diag(np.diag(map_))

try:
    with open(f2) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f2, ". USAGE:\n"+sys.argv[0]+" -f1 matrix1.json -f2 matrix2.json -o output -nodiag\n")
    exit()
    
try:
    your_map = np.array(data['map'])
except:
    print("Error: Can't get data from file", f2 ,"by 'map' key\n")
    exit()
    
your_map[np.isnan(your_map)] = 0    

if not diag:
    your_map = your_map - np.diag(np.diag(your_map))

if (np.shape(your_map) != np.shape(map_)):
    print('Error: matrices have different size (' + str(np.shape(map_)) +' and ' +str(np.shape(your_map)) + ')\n')
    exit()

frob1 = np.sqrt(sum(abs(map_.flatten())**2)) #Frobenius norm of a matrix
frob2 = np.sqrt(sum(abs(your_map.flatten())**2))

if frob1 != 0:
    mat1 = map_/frob1
if frob2 != 0:
    mat2 = your_map/frob2
mat = mat1-mat2

frob = np.sqrt(sum(abs(mat.flatten())**2))

print('The original matrices have Frobenius norms equal to', frob1, '(First matrix) and', frob2, '(Second matrix);')
print('The Frobenius norm of the difference of matrices is', frob, '.')

MAT = pd.DataFrame(data=mat[::-1, :], index=names[::-1], columns=names)

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MAT, annot=False, cmap="bwr", center=0)

plt.title('Difference matrix of normalized matrices with norm ' + str(round(frob,3)), fontsize=20)
try:
    fig.savefig(out_path + '.pdf')
    print("File",out_path + ".pdf created")
except:
    print("Error writing file",out_path + '.pdf')

new_data = {}
new_data['first_matrix_norm'] = frob1
new_data['second_matrix_norm'] = frob2
new_data['diff_of_matrix_norm'] = frob
new_data['difference'] = mat.tolist()
try:
    with open(out_path + '.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + ".json created\n")
except:
    print("Error writing file",out_path + '.json\n')
