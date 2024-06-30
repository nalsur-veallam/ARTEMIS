import sys
import numpy as np
import pylab as plt
import seaborn as sns
import pandas as pd
import json
from scipy.spatial import distance

width = .6
noseq = 0

print("\nCOMPARISON 2 SYSTEMS SCRIPT IS LAUNCHED\n")

if not ("-n1" in sys.argv and "-n2" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -n1 name1 -n2 name2\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-n1":
        name1 = sys.argv[i+1]
    if sys.argv[i] == "-n2":
        name2 = sys.argv[i+1]
        
map1_path =  "output/" + name1 + "/map/" + name1 + "_map"
map2_path =  "output/" + name2 + "/map/" + name2 + "_map"
intensity1_path = "output/" + name1 + "/analysis/" + name1 + "_intensity"
intensity2_path = "output/" + name2 + "/analysis/" + name2 + "_intensity"
out_path =  "output/" + name1 + "/analysis/" + name1 + "_with_" + name2

try:
    with open(map1_path + '.json') as json_file:
        data1 = json.load(json_file)
except:
    print("Error reading file", map1_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n1 name1 -n2 name2\n")
    exit()

try:
    map1 = np.array(data1['map'])
except:
    print("Error: Can't get data from file", map1_path + '.json',"by 'map' key\n")
    exit()
    
try:
    NResidues = int(data1['NResidues'])
except:
    print("Caution:, Can't get data from file", map1_path + '.json',"by 'NResidues' key. Continues without this data.")
    NResidues = len(map1)
    
try:
    names1 = np.array(data1['names'])
except:
    print("Caution:, Can't get data from file", map1_path + '.json',"by 'names' key. Continues without this data.")
    names1 = np.arange(1, NResidues+1)
    
try:
    real_numbers1 = np.array(data1['real_numbers'])
except:
    print("Caution:, Can't get data from file", map1_path + '.json',"by 'real_numbers' key. Continues without this data.")
    real_numbers1 = np.empty(NResidues)
    
try:
    with open(map2_path + '.json') as json_file:
        data2 = json.load(json_file)
except:
    print("Error reading file", map2_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n1 name1 -n2 name2\n")
    exit()

try:
    map2 = np.array(data2['map'])
except:
    print("Error: Can't get data from file", map2_path + '.json',"by 'map' key\n")
    exit()
    
if (np.shape(map1) != np.shape(map2)):
    print('Error: matrices have different size (' + str(np.shape(map1)) +' and ' +str(np.shape(map2)) + ')\n')
    exit()
    
try:
    names2 = np.array(data2['names'])
except:
    print("Caution:, Can't get data from file", map2_path + '.json',"by 'names' key. Continues without this data.")
    names2 = np.arange(1, NResidues+1)
    
try:
    real_numbers2 = np.array(data2['real_numbers'])
except:
    print("Caution:, Can't get data from file", map2_path + '.json',"by 'real_numbers' key. Continues without this data.")
    real_numbers2 = np.empty(NResidues)
    
frob1 = np.sqrt(sum(abs(map1.flatten())**2)) #Frobenius norm of a matrix
frob2 = np.sqrt(sum(abs(map2.flatten())**2))

if frob1 != 0:
    mat1 = map1/frob1
    i, j = np.indices(mat1.shape)
    mat1[i==j] = 0
if frob2 != 0:
    mat2 = map2/frob2
    i, j = np.indices(mat2.shape)
    mat2[i==j] = 0
mat = mat1-mat2

frob = np.sqrt(sum(abs(mat.flatten())**2))

labels = []
Entropy1 = []
Entropy2 = []
res_diff = []
res_jsd = []
res1_mean = np.mean(map1, axis=0)
res2_mean = np.mean(map2, axis=0)
res1_std = np.std(map1, axis=0)
res2_std = np.std(map2, axis=0)
for i in range(len(names1)):
    res1_distr = map1[i, :]/np.sqrt(sum(abs(map1[i, :].flatten())**2))
    res2_distr = map2[i, :]/np.sqrt(sum(abs(map2[i, :].flatten())**2))
    res_diff.append(np.sqrt(sum(abs((res1_distr - res2_distr).flatten())**2)))
    res_jsd.append(distance.jensenshannon(res1_distr, res2_distr))
    Entropy1.append(map1[i, i])
    Entropy2.append(map2[i, i])
    labels.append(names1[i] + "(" + str(real_numbers1[i]) + ")/" + names2[i] + "(" + str(real_numbers2[i]) + ")")

efrob1 = np.sqrt(sum(np.abs(Entropy1)**2)) #Frobenius norm of a matrix
efrob2 = np.sqrt(sum(np.abs(Entropy2)**2))

if efrob1 != 0:
    entropy1 = Entropy1/efrob1
if efrob2 != 0:
    entropy2 = Entropy2/efrob2
entropy = entropy1 - entropy2
efrob = np.sqrt(sum(abs(entropy.flatten())**2))

print('The original MI matrices have Frobenius norms equal to', frob1, '(First matrix) and', frob2, '(Second matrix);')
print('The diagonal elements of original MI matrices have Frobenius norms equal to', efrob1, '(First matrix) and', efrob2, '(Second matrix);')
print('The Frobenius norm of the difference of matrices (without the diagonals) is', frob, ';')
print('The Frobenius norm of the difference of diagonals of matrices is', efrob, '.\n')

MAT = pd.DataFrame(data=mat[::-1, :] , index=np.array(labels)[::-1], columns=np.array(labels))

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

sns.heatmap(MAT, annot=False, cmap="bwr", center=0)

plt.title('Difference matrix of '+ name1 + ' and '+ name2 + ' matrices with norm ' + str(round(frob,3)), fontsize=20)
try:
    fig.savefig(out_path + '_mapdiff.pdf')
    print("File",out_path + "_mapdiff.pdf created")
except:
    print("Error writing file",out_path + '_mapdiff.pdf')

INTENSITY = {}
INTENSITY["Residue (" + name1 + "/" + name2 + ")"] = labels
INTENSITY["Frobenius norm of difference per residue"] = res_diff
INTENSITY["JSD for residue between systems"] = res_jsd
INTENSITY["Mean MI for 1st system"] = res1_mean.tolist()
INTENSITY["Mean MI for 2nd system"] = res2_mean.tolist()
INTENSITY["Entropy for 1st system"] = Entropy1
INTENSITY["Entropy for 2nd system"] = Entropy2
INTENSITY["STD MI for 1st system"] = res1_std.tolist()
INTENSITY["STD MI for 2nd system"] = res2_std.tolist()
INTENSITY = pd.DataFrame(data=INTENSITY)

mutations = pd.DataFrame(data=names1!=names2, columns=['Mutation'])
INTENSITY = pd.concat([INTENSITY, mutations], axis=1)

INTENSITY = INTENSITY.sort_values(by=['Frobenius norm of difference per residue'], ascending=False, ignore_index=True)

INTENSITY = np.round(INTENSITY, 3)        

try:
    INTENSITY.to_csv(out_path + "_diff_table.csv")
    print("File",out_path + "_diff_table.csv created")
except:
    print("Error writing file",out_path + '_diff_table.csv')
    
INTENSITY = {}
INTENSITY["Residue (" + name1 + "/" + name2 + ")"] = labels
INTENSITY["Entropy difference"] = np.abs(np.array(Entropy1) - np.array(Entropy2))
INTENSITY["Frobenius norm of difference per residue"] = res_diff
INTENSITY["JSD for residue between systems"] = res_jsd
INTENSITY["Mean MI for 1st system"] = res1_mean.tolist()
INTENSITY["Mean MI for 2nd system"] = res2_mean.tolist()
INTENSITY["Entropy for 1st system"] = Entropy1
INTENSITY["Entropy for 2nd system"] = Entropy2
INTENSITY["STD MI for 1st system"] = res1_std.tolist()
INTENSITY["STD MI for 2nd system"] = res2_std.tolist()
INTENSITY = pd.DataFrame(data=INTENSITY)

mutations = pd.DataFrame(data=names1!=names2, columns=['Mutation'])
INTENSITY = pd.concat([INTENSITY, mutations], axis=1)

INTENSITY = INTENSITY.sort_values(by=['Entropy difference'], ascending=False, ignore_index=True)

INTENSITY = np.round(INTENSITY, 3)        

try:
    INTENSITY.to_csv(out_path + "_entropydiff_table.csv")
    print("File",out_path + "_entropydiff_table.csv created")
except:
    print("Error writing file",out_path + '_entropydiff_table.csv')

new_data = {}
new_data['name1'] = name1
new_data['name2'] = name2
new_data['first_matrix_norm'] = frob1
new_data['second_matrix_norm'] = frob2
new_data['diff_of_matrices_norm'] = frob
new_data['difference'] = mat.tolist()
new_data['names1'] = names1.tolist()
new_data['names2'] = names2.tolist()
new_data['NResidues'] = NResidues
new_data['real_numbers1'] = real_numbers1.tolist()
new_data['real_numbers2'] = real_numbers2.tolist()
new_data['first_entropy_norm'] = efrob1
new_data['second_entropy_norm'] = efrob2
new_data['diff_of_entropies_norm'] = efrob
new_data['mutations']=(names1!=names2).tolist()
try:
    with open(out_path + '_mapdiff.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + "_mapdiff.json created\n")
except:
    print("Error writing file",out_path + '_mapdiff.json\n')
    
try:
    with open(intensity2_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", intensity2_path, ". USAGE:\n"+sys.argv[0]+" -n1 name1 -n2 name2\n")
    exit()

try:
    frob2 = np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
    intensity2 = np.array(data['intensity'])/np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
except:
    print("Error: Can't get data from file", intensity2_path + '.json',"by 'intensity' key\n")
    exit()

try:
    with open(intensity1_path + '.json') as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", intensity1_path + '.json', ". USAGE:\n"+sys.argv[0]+" -n1 name1 -n2 name2\n")
    exit()

try:
    frob1 = np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
    intensity1 = np.array(data['intensity'])/np.sqrt(sum(abs(np.array(data['intensity']).flatten())**2))
except:
    print("Error: Can't get data from file", intensity1_path + '.json',"by 'intensity' key\n")
    exit()
    
if (np.shape(intensity1) != np.shape(intensity2) != NResidues):
    print('Error: Intensities have different size (' + str(np.shape(intensity1)) +' and ' +str(np.shape(intensity2)) + ')')
    exit()
    
norm = np.sqrt(sum(abs((intensity1 - intensity2).flatten())**2))

print('The original intensities have Frobenius norms equal to', frob1, '(First) and', frob2, '(Second);')
print('The Frobenius norm of the difference of intensities is', norm, '.')
print("Jensen-Shannon distance between intensities is", distance.jensenshannon(intensity1, intensity2), '\n')

INTENSITY = {}
INTENSITY["Residue (" + name1 + "/" + name2 + ")"] = labels
INTENSITY["Intensity difference"] = np.abs(intensity1 - intensity2)
INTENSITY["Entropy difference"] = np.abs(np.array(Entropy1) - np.array(Entropy2))
INTENSITY["Frobenius norm of difference per residue"] = res_diff
INTENSITY["JSD for residue between systems"] = res_jsd
INTENSITY["Intensity for 1st system"] = intensity1
INTENSITY["Intensity for 2nd system"] = intensity2
INTENSITY["Mean MI for 1st system"] = res1_mean.tolist()
INTENSITY["Mean MI for 2nd system"] = res2_mean.tolist()
INTENSITY["Entropy for 1st system"] = Entropy1
INTENSITY["Entropy for 2nd system"] = Entropy2
INTENSITY["STD MI for 1st system"] = res1_std.tolist()
INTENSITY["STD MI for 2nd system"] = res2_std.tolist()
INTENSITY = pd.DataFrame(data=INTENSITY)

mutations = pd.DataFrame(data=names1!=names2, columns=['Mutation'])
INTENSITY = pd.concat([INTENSITY, mutations], axis=1)

INTENSITY = INTENSITY.sort_values(by=['Intensity difference'], ascending=False, ignore_index=True)

INTENSITY = np.round(INTENSITY, 3)        

try:
    INTENSITY.to_csv(out_path + "_intensitydiff_table.csv")
    print("File",out_path + "_intensitydiff_table.csv created")
except:
    print("Error writing file",out_path + '_intensitydiff_table.csv')

labels = []
for i in range(len(names1)):
    labels.append(names1[i] + "(" + str(real_numbers1[i]) + ")\nand\n" + names2[i] + "(" + str(real_numbers2[i]) + ")")

INTENSITY = {}
INTENSITY["Intensity"] = (intensity1 - intensity2) / norm
INTENSITY["Residue"] = labels

colarr = INTENSITY["Intensity"] - np.min(INTENSITY["Intensity"])
colors = plt.cm.viridis(colarr/np.max(colarr))
colors = sns.color_palette(colors, as_cmap=True)

fig, axs = plt.subplots(figsize=(NResidues*width, 20), constrained_layout=True)
axs = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.title("Compare intensities of " + name1 + " and " + name2 + "with norm " + str(norm), fontsize=40)
try:
    fig.savefig(out_path+'_intensitydiff.pdf')
    print("File",out_path + "_intensitydiff.pdf created")
except:
    print("Error writing file",out_path + '_intensitydiff.pdf')

new_data = {}
new_data['first_inten_norm'] = frob1
new_data['second_inten_norm'] = frob2
new_data['diff_of_inten_norm'] = norm
new_data['difference'] =(intensity1 - intensity2).tolist()
new_data['names1'] = names1.tolist()
new_data['names2'] = names2.tolist()
new_data['NResidues'] = NResidues
new_data['real_numbers1'] = real_numbers1.tolist()
new_data['real_numbers2'] = real_numbers2.tolist()
try:
    with open(out_path + '_intensitydiff.json', 'w') as outfile:
        json.dump(new_data, outfile)
    print("File",out_path + "_intensitydiff.json created\n")
except:
    print("Error writing file",out_path + '_intensitydiff.json\n')
