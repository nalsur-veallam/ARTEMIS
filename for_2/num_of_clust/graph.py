import numpy as np
import pylab as plt
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

nums = []
for i in range(250, 350):
    cl2 = AgglomerativeClustering(n_clusters = None, distance_threshold=i*data['params'][1], compute_distances=True).fit(MAP2)
    nums.append(cl2.n_clusters_)
    print(i)
    
fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.plot(np.linspace(250,350, 100),nums)


plt.title('Nums clusters from distance')
plt.grid()
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('numclust_3')
