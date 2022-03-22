import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans

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

MAP = MAP1
M = MAP2
M1 = MAP1
M2 = MAP2
MAP1 = abs(MAP1)
sc1 = SpectralClustering(10, affinity='precomputed', n_init=100, assign_labels='discretize')
sc1.fit_predict(MAP1)
for i in range(N):
    for j in range(N):
        if sc1.labels_[i] == sc1.labels_[j]:
            MAP1[i][j] = sc1.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(MAP1, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI1 - SpectralClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustspeecmi1')


sc2 = SpectralClustering(10, affinity='precomputed', n_init=100, assign_labels='discretize')
sc2.fit_predict(MAP2)
for i in range(N):
    for j in range(N):
        if sc2.labels_[i] == sc2.labels_[j]:
            MAP2[i][j] = sc2.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(MAP2, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI2 - SpectralClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustspeecmi2')


cl = AgglomerativeClustering(n_clusters = 10).fit(MAP)
for i in range(N):
    for j in range(N):
        if cl.labels_[i] == cl.labels_[j]:
            MAP[i][j] = cl.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(MAP, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI1 - AgglomerativeClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustagglmi1')



cl2 = AgglomerativeClustering(n_clusters = 10).fit(M)
for i in range(N):
    for j in range(N):
        if cl2.labels_[i] == cl2.labels_[j]:
            M[i][j] = cl2.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(M, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI2 - AgglomerativeClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustagglmi2')



km1 = KMeans(n_clusters=10, random_state=0).fit(M1)
for i in range(N):
    for j in range(N):
        if km1.labels_[i] == km1.labels_[j]:
            M1[i][j] = km1.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(M1, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI1 - KmeansClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustkmmi1')


km2 = KMeans(n_clusters=10, random_state=0).fit(M2)
for i in range(N):
    for j in range(N):
        if km2.labels_[i] == km2.labels_[j]:
            M2[i][j] = km2.labels_[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(M2, extent=[0, N, 0, N], cmap = 'tab20b')


fig.colorbar(p1, ax=axs)
plt.title('MI2 - KmeansClustering')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('clustkmmi2')

