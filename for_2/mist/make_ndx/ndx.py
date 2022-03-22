import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
from sklearn.cluster import AgglomerativeClustering

def show_number(line, i):
    numbers = []
    while line[i] != '\n':
        k = []
        number = 0
        while(line[i] != ' ' and line[i] != '\n' and line[i] != '-'):
            k.append(line[i])
            i+=1
        for j in range(len(k)):
            number+=(10**(len(k)-j - 1))*int(k[j])
        numbers.append(number)
        if line[i] != '\n':
            i+=1
    return numbers

B = True
A = False
D = False


fig, ax = plt.subplots()

data_mi = np.loadtxt('MIE.txt')
N = 965
clust = 6

BOND = data_mi [:, 0]
B_1 = data_mi [:, 1]
B_2 = data_mi [:, 2]
S_2 = data_mi [:, 3]


MAP1 = np.zeros((N, N))
MAP2 = np.zeros((N, N))
for i in range(len(BOND)):
    if(int(B_1[i]) <= int(B_2[i]) and BOND[i] == 33):
        MAP2[int(B_1[i])-1][int(B_2[i])-1] = S_2[i]
        MAP2[int(B_2[i])-1][int(B_1[i])-1] = S_2[i]

cl2 = AgglomerativeClustering(n_clusters = clust).fit(MAP2)
with open('clust.ndx', 'w') as outfile:
    for i in range(clust):
        outfile.write('[' + str(i + 1) + ']\n')
        B = True
        A = False
        D = False
        BAT = []
        f3 = open('PARENT_MIE_topology.txt', 'r')
        for j in f3:
            if(j[0] == '#'):
                    if(j[1] == 'b'):
                        B = True
                        A = False
                        D = False
                    elif(j[1] == 'a'):
                        B = False
                        A = True
                        D = False
                    elif(j[1] == 'd'):
                        B = False
                        A = False
                        D = True
            else:
                coord = show_number(j, 0)
                if D:
                    if cl2.labels_[coord[0] - 1] == i:
                        for k in coord[1:5]:
                            BAT.append(k)
                            
        f3.close()
        outfile.write(" ".join(map(str,np.unique(BAT))) + '\n')
        
        
        
        
        
        
        
        
        
