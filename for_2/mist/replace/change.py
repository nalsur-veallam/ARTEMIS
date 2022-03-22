import numpy as np

data = np.loadtxt('MIST.txt')


for i in range(len(data[:,0])):
    if data[i][0] == 32:
        data[i][0] = 23
        data[i][1], data[i][2] = data[i][2], data[i][1]
    if data[i][0] == 31:
        data[i][0] = 13
        data[i][1], data[i][2] = data[i][2], data[i][1]
    if data[i][0] == 21:
        data[i][0] = 12
        data[i][1], data[i][2] = data[i][2], data[i][1]


with open('MIE.txt', 'w') as f:
    for row in data:
        f.write(str(int(row[0])) + ' ' + str(int(row[1])) + ' ' + str(int(row[2])) + ' ' + str(row[3]) + '\n')
