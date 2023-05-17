import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json
from sklearn.linear_model import LinearRegression

diag = True
covar = False

rnames = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
masses = [71.0788, 156.1875, 114.1038, 115.0886, 103.1388, 103.1388, 128.1307, 57.0519, 137.1411, 113.1594, 113.1594, 128.1741, 131.1926, 147.1766, 97.1167, 87.0782, 101.1051, 186.2132, 163.1760, 99.1326]
charges = [0, +3, 0, -2, 0, -2, 0, 0, +2, 0, 0, +1, 0, 0, 0, 0, 0, 0, 0, 0]

print("\nSCRIPT FOR CALCULATING LINEAR REGESSION FOR MI MATRIX IS LAUNCHED\n")

if not ("-n" in sys.argv and len(sys.argv) >= 5):
    print("USAGE:\n"+sys.argv[0]+" -n name -nodiag -dist dist_matrix -covar covar_matrix\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n":
        name = sys.argv[i+1]
    if sys.argv[i] == "-nodiag":
        diag = False
    if sys.argv[i] == "-dist":
        dist_path = sys.argv[i+1]
    if sys.argv[i] == "-covar":
        covar = True
        covar_path = sys.argv[i+1]
        
        
map_path =  "output/" + name + "/map/" + name + "_map"
out_path =  "output/" + name + "/map/" + name + "_LinReg_map"

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
    NResidues = len(names)
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'names' key. Continues without this data.")
    exit()
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", map_path + '.json',"by 'real_numbers' key. Continues without this data.")
    real_numbers = np.empty(NResidues)
    
try:
    with open(dist_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", dist_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
    exit()

try:
    dist_map = np.array(data['map'])
except:
    print("Error: Can't get data from file", dist_path + '.json',"by 'map' key\n")
    exit()
    
if covar:
    try:
        with open(covar_path + '.json') as json_file:
            data = json.load(json_file)
    except:
        print("Error reading file", covar_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n name -nodiag -norm\n")
        exit()

    try:
        covar_map = np.array(data['map'])
    except:
        print("Error: Can't get data from file", covar_path + '.json',"by 'map' key\n")
        exit()

mass_map = np.zeros((NResidues, NResidues))

labels = []
for i in range(len(names)):
    labels.append(names[i] + " (" + str(real_numbers[i]) + ")")
    
Masses = np.zeros(NResidues)
for i in range(20):
    resname = rnames[i]
    Masses[np.argwhere(names == resname).reshape(-1)] = masses[i]
    
i, j = np.indices(map_.shape)
mass_map[i, j] = Masses[i] * Masses[j]

charge_map = np.zeros((NResidues, NResidues))

Charges = np.zeros(NResidues)
for i in range(20):
    resname = rnames[i]
    Charges[np.argwhere(names == resname).reshape(-1)] = charges[i]
    
i, j = np.indices(map_.shape)
charge_map[i, j] = np.abs(Charges[i] * Charges[j])

if not diag:
        map_[i,i] = 0
        charge_map[i,i] = 0
        dist_map[i,i] = 0
        mass_map[i,i] = 0
        if covar:
            covar_map[i,i] = 0

map_ = map_ / np.sqrt(sum(abs(map_.flatten())**2))
dist_map = dist_map / np.sqrt(sum(abs(dist_map.flatten())**2))
mass_map = mass_map / np.sqrt(sum(abs(mass_map.flatten())**2))
charge_map = charge_map / np.sqrt(sum(abs(charge_map.flatten())**2))
if covar:
    covar_map = covar_map / np.sqrt(sum(abs(covar_map.flatten())**2))

m = map_.flatten()
m1 = dist_map.flatten().reshape(-1, 1)
m2 = mass_map.flatten().reshape(-1, 1)
m3 = charge_map.flatten().reshape(-1, 1)
if covar:
    m4 = covar_map.flatten().reshape(-1, 1)

    X = np.concatenate([m1, m2, m3, m4], axis=1)
else:
    X = np.concatenate([m1, m2, m3], axis=1)
Y = m
model = LinearRegression(fit_intercept=True)
model.fit(X, Y)

coefs = model.coef_
frob = np.sqrt(sum(abs((Y - model.predict(X)).flatten())**2))

print("The coefficient for Distance matrix is", round(coefs[0], 4))
print("The coefficient for Mass matrix is", round(coefs[1], 4))
print("The coefficient for Charge matrix is", round(coefs[2], 4))
print("The intercept is", round(model.intercept_, 4))
if covar:
    print("The coefficient for Covar matrix is", round(coefs[3], 4))
print("Frobenius norm of difference is", round(frob, 4), '\n')

map_ = (Y - model.predict(X)).reshape(NResidues, NResidues)
MIE = pd.DataFrame(data=map_[::-1, :], index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MIE, annot=False, cmap="bwr")

plt.title('Difference matrix in Linear Regression for ' + name + 'with norm ' + str(round(frob, 4)), fontsize=20)
try:
    fig.savefig(out_path + '.pdf')
    print("File",out_path + ".pdf created")
except:
    print("Error writing file",out_path + '.pdf\n')
    
new_data = {}
new_data['names'] = names.tolist()
new_data['NResidues'] = NResidues
new_data['name'] = name
new_data['map'] = map_.tolist()
new_data['real_numbers'] = real_numbers.tolist()
new_data['coefficients'] = coefs.tolist()
new_data['norm'] = frob
new_data['covar'] = covar

try:
    with open(out_path + '.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + ".json created\n")
except:
    print("Error writing file",out_path + '.json\n')
