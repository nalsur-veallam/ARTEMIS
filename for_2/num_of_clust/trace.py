import numpy as np
import pylab as plt
from sklearn.cluster import AgglomerativeClustering

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


TR = []
for i in range(20, 30):
    tr = 0
    A = 0
    cl = AgglomerativeClustering(n_clusters = i).fit(MAP2)
    for j in range(i):
        a = 0
        for l in range(N):
            for k in range(N):
                if(cl.labels_[k] == cl.labels_[l] and cl.labels_[l] == j):
                    a += MAP2[k][l]
        B = []
        for h in range(i):
            b = 0
            for l in range(N):
                for k in range(N):
                    if(cl.labels_[k] != cl.labels_[l] and cl.labels_[l] == j and cl.labels_[k] == h):
                        b += MAP2[k][l]
            B.append(b)
        b = max(B)
        print(b, a)
        a = (b - a)/max(b, a)
        A += a
    TR.append(A/np.log(i))
    print(i)
    
fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.plot(np.linspace(20,30, 10),TR)


plt.title('Trace of clusters')
plt.xlabel('N - Num of clusters')
plt.ylabel('TR/log(N)')
plt.grid()
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('TRwLog')
