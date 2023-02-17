import json
import numpy as np
import sys

if not ("-matname" in sys.argv and "-f" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -f matrix.json -matname matrix_name")
    exit()

for i in range(1, 5) :
    if sys.argv[i] == "-matname":
        matrix_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]
        
if not (matrix_name and path):
    print("USAGE:\n"+sys.argv[0]+" -f matrix.json -n name -matname matrix_name")
    exit()
    
    
with open(path) as json_file:
    data = json.load(json_file)

your_map = np.array(data[matrix_name])

frob = np.sqrt(sum(abs(your_map.flatten())**2)) #Frobenius norm of a matrix
print('The Frobenius norm of your matrix is', frob)
