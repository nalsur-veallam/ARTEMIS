mass = []
first = True
f1 = open('atoms.txt', 'r')
f2 = open('groups.txt', 'w')
for line in f1:
    if(line[0] == ';'):
        if first:
            f2.write('[' + line[:-1] + ']\n')
            first = False
        else:
            f2.write(" ".join(map(str,mass)) + '\n')
            f2.write('[' + line[:-1] + ']\n')
            mass = []
    else:
        i = 0
        atom = []
        while line[i] == ' ':
            i+=1
        while line[i] != ' ':
            atom.append(line[i])
            i+=1
        k = 0
        for j in range(len(atom)):
            k+=(10**(len(atom)-j - 1))*int(atom[j])
        mass.append(k)
                
f2.write(" ".join(map(str,mass)))
f1.close()
f2.close()
