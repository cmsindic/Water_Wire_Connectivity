import os

dirs = [d for d in os.scandir() if 'ases' in d.name]

n = 0
n0 = 0

for d in dirs:

    for file in os.scandir(d):
        if 'wire' in file.name:
            n += 1
            lines = []
            with open(file,'r') as f:
                for line in f:
                    line = line.split(',')
                    lines.append(line)

            if len(lines) == 0:
                continue

            ncontacted = int(lines[6][0])

            if ncontacted == 0:
                n0 += 1

print(n,n0,100*n0/n)
