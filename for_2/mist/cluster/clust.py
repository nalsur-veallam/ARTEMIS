import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
from sklearn.cluster import AgglomerativeClustering
import json


with open('data.txt') as json_file:
    data = json.load(json_file)

n_clusters = 16

fig, ax = plt.subplots()

data_mi = np.loadtxt('MIE.txt')
N = 2898

BOND = data_mi [:, 0]
B_1 = data_mi [:, 1]
B_2 = data_mi [:, 2]
S_2 = data_mi [:, 3]

for i in range(len(BOND)):
    if (BOND[i] == 12):
        B_2[i] += N/3 + 1
    elif (BOND[i] == 22):
        B_2[i] += N/3 + 1
        B_1[i] += N/3 + 1
    elif (BOND[i] == 13):
        B_2[i] += 2*N/3 + 1
    elif (BOND[i] == 23):
        B_1[i] += N/3 + 1
        B_2[i] += 2*N/3 + 1
    elif (BOND[i] == 33):
        B_1[i] += 2*N/3 + 1
        B_2[i] += 2*N/3 + 1



MAP1 = np.zeros((N, N))
MAP2 = np.zeros((N, N))
for i in range(len(BOND)):
    if(int(B_1[i]) <= int(B_2[i])):
        MAP2[int(B_1[i])-1][int(B_2[i])-1] = S_2[i]
        MAP2[int(B_2[i])-1][int(B_1[i])-1] = S_2[i]

for n_clusters in range(2, 20):
    cl2 = AgglomerativeClustering(n_clusters).fit(MAP2)
    print('N = ' + str(cl2.n_clusters_))
    for i in range(N):
        for j in range(N):
            if cl2.labels_[i] == cl2.labels_[j]:
                MAP2[i][j] = cl2.labels_[i]


    fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

    p1 = axs.imshow(MAP2, extent=[0, N, 0, N], cmap = 'tab20b')


    fig.colorbar(p1, ax=axs)
    plt.title('MI2 - AgglomerativeClustering')

    fig.set_figwidth(10)
    fig.set_figheight(10)

    fig.savefig(str(n_clusters) + 'clustagglmi2all')
