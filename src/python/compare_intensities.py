import sys
import json
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd

width = .6
print("\nINTENSITIES COMPARISON SCRIPT IS LAUNCHED\n")

if not ("-f1" in sys.argv and "-f2" in sys.argv and "-o" in sys.argv and len(sys.argv) == 7):
    print("USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -o output\n")
    exit()
    
for i in range(1, 7) :
    if sys.argv[i] == "-f1":
        f1 = sys.argv[i+1]
    if sys.argv[i] == "-o":
        out_path = sys.argv[i+1]
    if sys.argv[i] == "-f2":
        f2 = sys.argv[i+1]

try:
    with open(f2) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f2, ". USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -o output\n")
    exit()

try:
    frob2 = np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
    intensity2 = np.array(data['intensity'])/np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
except:
    print("Error: Can't get data from file", f2,"by 'intensity' key\n")
    exit()

try:
    with open(f1) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f1, ". USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -o output\n")
    exit()

try:
    frob1 = np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
    intensity1 = np.array(data['intensity'])/np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
except:
    print("Error: Can't get data from file", f1,"by 'intensity' key\n")
    exit()

try:
    NResidues = int(data['NResidues'])
except:
    print("Caution:, Can't get data from file",f1 ,"by 'NResidues' key. Continues without this data.")
    NResidues = len(intensity1)
    
try:
    names = np.array(data['names'])
except:
    print("Caution:, Can't get data from file", f1,"by 'names' key. Continues without this data.")
    names = np.arange(1, NResidues+1)
    
try:
    real_numbers = np.array(data['real_numbers'])
except:
    print("Caution:, Can't get data from file", f1,"by 'real_numbers' key. Continues without this data.")
    real_numbers = np.empty(NResidues)

if (np.shape(intensity1) != np.shape(intensity2)):
    print('Error: Intensities have different size (' + str(np.shape(intensity1)) +' and ' +str(np.shape(intensity2)) + ')')
    exit()

new_names = []
for i in range(NResidues):
    new_names.append(names[i] + "\n(" + str(real_numbers[i]) +")")

norm = np.sqrt(sum(abs((intensity1 - intensity2).flatten())**2))

print('The original intensities have Frobenius norms equal to', frob1, '(First) and', frob2, '(Second);')
print('The Frobenius norm of the difference of intensities is', norm, '.')

INTENSITY = {}
INTENSITY["Intensity"] = (intensity1 - intensity2) / norm
INTENSITY["Residue"] = new_names

colarr = INTENSITY["Intensity"] - np.min(INTENSITY["Intensity"])
colors = plt.cm.viridis(colarr/np.max(colarr))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title('Compare intensity with norm ' + str(norm), fontsize=40)
try:
    fig.savefig(out_path+'.pdf')
    print("File",out_path + ".pdf created")
except:
    print("Error writing file",out_path + '.pdf')

new_data = {}
new_data['first_inten_norm'] = frob1
new_data['second_inten_norm'] = frob2
new_data['diff_of_inten_norm'] = norm
new_data['difference'] =(intensity1 - intensity2).tolist()
try:
    with open(out_path + '.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + ".json created\n")
except:
    print("Error writing file",out_path + '.json\n')
