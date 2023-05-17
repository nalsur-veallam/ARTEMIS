import sys
import json
import numpy as np
import pylab as plt
from scipy.spatial import distance

hist=False
prob=False
width = .6
nbins=10

print("\nJENSEN SHANNON DISTANCE CALCULATION SCRIPT IS LAUNCHED\n")

if not ("-f1" in sys.argv and "-arrn1" in sys.argv and "-arrn2" in sys.argv and "-f2" in sys.argv and len(sys.argv) >= 9):
    print("USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -hist output -prob -arrn1 arr1name -arrn2 arr2name\n")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f1":
        f1 = sys.argv[i+1]
    if sys.argv[i] == "-hist":
        hist = True
        out_path = sys.argv[i+1]
    if sys.argv[i] == "-prob":
        prob = True
    if sys.argv[i] == "-f2":
        f2 = sys.argv[i+1]
    if sys.argv[i] == "-arrn2":
        arr2name = sys.argv[i+1]
    if sys.argv[i] == "-arrn1":
        arr1name = sys.argv[i+1]

try:
    with open(f1) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f1, ". USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -hist output -prob -arrn1 arr1name -arrn2 arr2name\n")
    exit()

try:
    intensity1 = np.array(data[arr1name]).reshape(-1)
except:
    print("Error: Can't get data from file", f1,"by",arr1name,"key\n")
    exit()

try:
    with open(f2) as json_file:
        data = json.load(json_file)
except:
    print("Error reading file", f2, ". USAGE:\n"+sys.argv[0]+" -f1 intensity1.json -f2 intensity2.json -hist output -prob -arrn1 arr1name -arrn2 arr2name\n")
    exit()

try:
    intensity2 = np.array(data[arr2name]).reshape(-1)
except:
    print("Error: Can't get data from file", f2,"by",arr2name,"key\n")
    exit()

if prob:
    if (np.shape(intensity1) != np.shape(intensity2)):
        print('Error: Probability vectors have different size (' + str(np.shape(intensity1)) +' and ' +str(np.shape(intensity2)) + ')\n')
        exit()
    print(" Jensen-Shannon distance is", distance.jensenshannon(intensity1, intensity2), '\n')
    
else:
    if hist:
        fig, axs = plt.subplots(figsize=(20, 15), constrained_layout=True)
        axs.hist(intensity1, bins=50)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.title('First histogram', fontsize=40)
        try:
            fig.savefig(out_path + '_hist1.pdf')
            print("File",out_path + "_hist1.pdf created")
        except:
            print("Error writing file",out_path + '_hist1.pdf')
        
        fig, axs = plt.subplots(figsize=(20, 15), constrained_layout=True)
        axs.hist(intensity2, bins=50)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.title('Second histogram', fontsize=40)
        try:
            fig.savefig(out_path + '_hist2.pdf')
            print("File",out_path + "_hist2.pdf created")
        except:
            print("Error writing file",out_path + '_hist2.pdf')
    
    min_ = np.min([intensity1, intensity2])
    max_ = np.max([intensity1, intensity2])
    
    probs1 = np.histogram(intensity1, bins=nbins, range=[min_, max_])[0]
    probs2 = np.histogram(intensity2, bins=nbins, range=[min_, max_])[0]
    print(" Jensen-Shannon distance is", distance.jensenshannon(probs1, probs2), '\n')
    
    
