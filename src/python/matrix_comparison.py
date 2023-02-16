import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

if not ("-matname" in sys.argv and "-f" in sys.argv and "-n" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f matrix.json -n name -matname matrix_name")
    exit()

for i in range(1, 7) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-matname":
        matrix_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
        
if not (name and matrix_name and path):
    print("USAGE:\n"+sys.argv[0]+" -f matrix.json -n name -matname matrix_name")
    exit()
    
    
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/analysis/" + name

with open(map_path + '.json') as json_file:
    data = json.load(json_file)

map_ = np.array(data['map'])
map_ = map_ - np.diag(np.diag(map_))
names = np.array(data['names'])
NResidues = data['NResidues']

with open(path) as json_file:
    data = json.load(json_file)

your_map = np.array(data[matrix_name])
your_map = your_map - np.diag(np.diag(your_map))

if (np.shape(your_map) != np.shape(map_)):
    print('Error: matrices have different size (' + str(np.shape(map_)) +' and ' +str(np.shape(your_map)) + ')')
    exit()

frob1 = np.sqrt(sum(abs(map_.flatten())**2)) #Frobenius norm of a matrix
frob2 = np.sqrt(sum(abs(your_map.flatten())**2))

mat1 = map_/frob1
mat2 = your_map/frob2
mat = mat1-mat2

frob = np.sqrt(sum(abs(mat.flatten())**2))

print('The original matrices have Frobenius norms equal to', frob1, '(PARENT matrix) and', frob2, '(custom matrix);')
print('The Frobenius norm of the difference of matrices is', frob, '.')

MAT = pd.DataFrame(data=mat[::-1, :], index=names[::-1], columns=names)

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MAT, annot=False, cmap="bwr", center=0)

plt.title('Difference matrix of normalized matrices for ' + name, fontsize=20)
fig.savefig(out_path + '_matrix_comparison.pdf')

new_data = {}
new_data['name'] = name
new_data['PARENT_matrix_norm'] = frob1
new_data['custom_matrix_norm'] = frob2
new_data['diff_of_matrix_norm'] = frob
with open(out_path + '_matrix_comparison.json', 'w') as outfile:
    json.dump(new_data, outfile)
