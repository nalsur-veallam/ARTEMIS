import numpy as np

def show_number(line, i):
    numbers = []
    while(line[i] != '\n' and i != len(line) - 1):
        k = []
        number = 0
        while(line[i] != ' ' and line[i] != '\n' and i != len(line) - 1):
            k.append(line[i])
            i+=1
        for j in range(len(k)):
            number+=(10**(len(k)-j - 1))*int(k[j])
        numbers.append(number)
        if(line[i] != '\n' and i != len(line) - 1):
            i+=1
    return numbers

N = 484
first = True

f1 = open('groups.txt', 'r')
f2 = open('groups_2.txt', 'w')
for line in f1:
    if(line[0] == '['):
        if first:
            f2.write(line)
            first = False
            name = line
        else:
            f2.write('\n' + line)
            name = line
    else:
        atoms = np.array(show_number(line, 0)) + N
        f2.write(line[:-1])
        f2.write('\n' + name)
        f2.write(" ".join(map(str,atoms)))
f1.close()
f2.close()
