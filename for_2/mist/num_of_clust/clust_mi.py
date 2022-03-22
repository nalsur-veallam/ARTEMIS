import numpy as np
import pylab as plt
from sklearn.cluster import AgglomerativeClustering

data_mi = np.loadtxt('MIE.txt')
N = 965

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


MI = []
for i in range (10, 21):
    mi = 0;
    cl = AgglomerativeClustering(n_clusters = i).fit(MAP2)
    for j in range(N):
        for k in range(N):
            if cl.labels_[k] != cl.labels_[j]:
                mi += MAP2[k][j]
    MI.append(mi/np.log(i))
    print(i)

    
fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.plot(np.linspace(10,20, 11),MI)

plt.title('MI/log(N) from num of clusters')
plt.xlabel('N - Num of clusters')
plt.ylabel('MI/log(N)')
plt.grid()
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('MI:log(N)_2')
