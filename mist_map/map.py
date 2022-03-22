import numpy as np
import pylab as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots()

data_mi = np.loadtxt('MIE.txt', skiprows = 2898)
N = 2898

BOND = data_mi [:, 0]
B_1 = data_mi [:, 1]
B_2 = data_mi [:, 2]
S_1 = data_mi [:, 3]
S_2 = data_mi [:, 4]

for i in range(len(BOND)):
    if (BOND[i] == 12):
        B_2[i] += 967
    elif (BOND[i] == 22):
        B_2[i] += 967
        B_1[i] += 967
    elif (BOND[i] == 13):
        B_2[i] += 1933
    elif (BOND[i] == 23):
        B_1[i] += 967
        B_2[i] += 1933
    elif (BOND[i] == 33):
        B_1[i] += 1933
        B_2[i] += 1933


MAP1 = np.zeros((N, N))
MAP2 = np.zeros((N, N))
for i in range(len(BOND)):
    if(int(B_1[i]) <= int(B_2[i])):
        MAP1[N - int(B_1[i])][int(B_2[i]) - 1] = S_1[i]
        MAP2[N - int(B_1[i])][int(B_2[i]) - 1] = S_2[i]
        MAP1[N - int(B_2[i])][int(B_1[i]) - 1] = S_1[i]
        MAP2[N - int(B_2[i])][int(B_1[i]) - 1] = S_2[i]


fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.imshow(MAP1, extent=[0, N, 0, N], cmap = 'PuOr')


fig.colorbar(p1, ax=axs)
plt.title('MIE for all degrees of freedom')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('mie_map')
