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


B = True
A = False
D = False
BAT = []
first = True
f1 = open('groups_2.txt', 'r')
f2 = open('newgroups.txt', 'w')
for line in f1:
    if(line[0] == '['):
        if first:
            f2.write(line)
            first = False
        else:
            f2.write(" ".join(map(str,BAT)) + '\n')
            f2.write(line)
            BAT = []
    else:
        B = True
        A = False
        D = False
        f3 = open('PARENT_MIE_topology.txt', 'r')
        for j in f3:
            if(j[0] == '#'):
                    if(j[1] == 'b'):
                        B = True
                        A = False
                        D = False
                    elif(j[1] == 'a'):
                        B = False
                        A = True
                        D = False
                    elif(j[1] == 'd'):
                        B = False
                        A = False
                        D = True
            else:
                atoms = line.split()
                coord = show_number(j, 0)
                if B:
                    if((str(coord[1]) in atoms) and (str(coord[2]) in atoms)):
                        BAT.append(coord[0])
                elif A:
                    if((str(coord[1]) in atoms) and (str(coord[2]) in atoms) and (str(coord[3]) in atoms)):
                        BAT.append(coord[0] + 967)
                elif D:
                    if((str(coord[1]) in atoms) and (str(coord[2]) in atoms) and (str(coord[3]) in atoms) and (str(coord[4]) in atoms)):
                        BAT.append(coord[0] + 1933)
        
        f3.close()
                
                
f2.write(" ".join(map(str,BAT)))
f1.close()
f2.close()
