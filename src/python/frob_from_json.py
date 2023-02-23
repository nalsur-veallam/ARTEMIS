import json
import numpy as np
import sys

print("\nSCRIPT FOR CALCULATION OF THE FROBENIUS NORM IS LAUNCHED\n")

if not ("-matname" in sys.argv and "-f" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -f matrix.json -matname matrix_name\n")
    exit()

for i in range(1, 5) :
    if sys.argv[i] == "-matname":
        matrix_name = sys.argv[i+1]
    if sys.argv[i] == "-f":
        path = sys.argv[i+1]    
    
try:
    with open(path) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", path, ". USAGE:\n"+sys.argv[0]+" -f matrix.json -matname matrix_name\n")
    exit()

try:
    your_map = np.array(data[matrix_name])
except:
    print("Error: Can't get data from file", path,"by", matrix_name, "key\n")
    exit()

frob = np.sqrt(sum(abs(your_map.flatten())**2)) #Frobenius norm of a matrix
print('The Frobenius norm of your matrix is', frob, '\n')
