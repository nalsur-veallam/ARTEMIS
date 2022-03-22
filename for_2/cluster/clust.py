import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
from sklearn.cluster import AgglomerativeClustering
import json


with open('data.txt') as json_file:
    data = json.load(json_file)



fig, ax = plt.subplots()

data_mi = np.loadtxt('MIE.txt', skiprows = 3735521)
N = 965

BOND = data_mi [:, 0]
B_1 = data_mi [:, 1]
B_2 = data_mi [:, 2]
S_1 = data_mi [:, 3]
S_2 = data_mi [:, 4]


MAP1 = np.zeros((N, N))
MAP2 = np.zeros((N, N))
for i in range(len(BOND)):
    if(int(B_1[i]) <= int(B_2[i])):
        MAP1[int(B_1[i])-1][int(B_2[i])-1] = S_1[i]
        MAP1[int(B_2[i])-1][int(B_1[i])-1] = S_1[i]
        MAP2[int(B_1[i])-1][int(B_2[i])-1] = S_2[i]
        MAP2[int(B_2[i])-1][int(B_1[i])-1] = S_2[i]

cl2 = AgglomerativeClustering(n_clusters = 35).fit(MAP2)
print('N = ' + str(cl2.n_clusters_))
for i in range(N):
    for j in range(N):
        if cl2.labels_[i] == cl2.labels_[j]:
            MAP2[i][j] = cl2.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(MAP2, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI2 - AgglomerativeClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('35clustagglmi2')
