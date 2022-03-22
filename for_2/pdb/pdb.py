import numpy as np

def show_number(line, i):
    numbers = []
    while line[i] != ' ':
        k = []
        number = 0
        while(line[i] != ' '):
            k.append(line[i])
            i-=1
        for j in range(len(k)):
            number+=(10**(j))*int(k[j])
        numbers.append(number)
        if line[i] != ' ':
            i-=1
    return number

groups = []
group = []
first = True
f1 = open('clust.ndx', 'r')
for line in f1:
    if first:
        first = False
        ng = 1
    else:
        if line[0] == '[':
            groups.append(group)
            group = []
            ng += 1
        else:
            group = line.split()
groups.append(group)


f3 = open('new_1.pdb', 'w')

with open('ev1.pdb', 'r') as f2:
    for lin in f2:
        if lin[0] == 'A':
            num = show_number(lin, 10)
            d = True
            n = 0
            while d and n < ng:
                n += 1
                for k in range(len(groups[n - 1])):
                    if int(num) == int(groups[n - 1][k]):
                        d = False
                        lin = lin[:25] + str(n) + lin[26:]
                        lin = lin[:24] + ' ' + lin[25:]
                        break
        f3.write(lin)
f3.close()
f1.close()
