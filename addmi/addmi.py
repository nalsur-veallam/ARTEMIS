import numpy as np


def show_number(line, i):
    numbers = []
    while line[i] != '\n':
        k = []
        number = 0
        while(line[i] != ' ' and line[i] != '\n' and line[i] != '-'):
            k.append(line[i])
            i+=1
        for j in range(len(k)):
            number+=(10**(len(k)-j - 1))*int(k[j])
        numbers.append(number)
        if line[i] != '\n':
            i+=1
    return numbers

count = 0
B = True
A = False
D = False

f1 = open('groups.txt', 'r')
for line in f1:
    if(line[0] == '['):
        count += 1


f1.close()
f1 = open('groups.txt', 'r')

addmid = []
addmib = []
lines = f1.readlines()
for i in range(1, count):
    B = True
    A = False
    D = False
    f3 = open('PARENT_MIE_topology.txt', 'r')
    for j in f3:
        if(j[0] == '#'):
                if(j[1] == 'b'):
                    B = True
                    D = False
                elif(j[1] == 'a'):
                    D = False
                    B = False
                elif(j[1] == 'd'):
                    B = False
                    D = True
        else:
            atoms1 = lines[2*i-1].split()
            atoms2 = lines[2*i+1].split()
            coord = show_number(j, 0)
            if D:
                if((str(coord[1]) in atoms1) and (str(coord[2]) in atoms1) and (str(coord[3]) in atoms2) and (str(coord[4]) in atoms2)):
                    addmid.append(coord[0] + 1933)
            elif B:
                if((str(coord[1]) in atoms1) and (str(coord[2]) in atoms2)):
                    addmib.append(coord[0])
    
    f3.close()
print(addmib)
print(addmid)
f1.close()
