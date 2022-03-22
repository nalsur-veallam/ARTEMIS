import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
import scipy.stats as stats
import json

data_mi = np.loadtxt('MIE.txt')

BOND = data_mi [:, 0]
B_1 = data_mi [:, 1]
B_2 = data_mi [:, 2]
S_2 = data_mi [:, 3]

fig, axs = plt.subplots(figsize=(10,10), constrained_layout=True)

p1 = axs.hist(S_2, bins= 50)


plt.title('MIE hist')
plt.show()

fig.set_figwidth(10)
fig.set_figheight(10)

fig.savefig('hist')

normal = stats.norm
params = normal.fit(S_2)
data = {}
data['params'] = params
data['variance'] = np.var(S_2)
data['mean'] = np.mean(S_2)
with open('data.txt', 'w') as outfile:
    json.dump(data, outfile)
