#!/usr/bin/env python3

f1 = open('target_sequence.txt', 'r') #file with a sequence
f2 = open('iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3', 'r') #gff3 file with TE coordinates
f3 = open('target_sequence_masked.txt', 'w') #masked sequence
f4 = open('target_position.txt', 'r')
f5 = open('SSRs.txt', 'r')
chr = ''
begin_reg = 0
end_reg = 0

for line in f4:
    if line.strip() != '':
        lss = line.strip().split()
        begin_reg = int(lss[2])
        end_reg = int(lss[3])
        chr = lss[1]

intervals = []

for line in f2:
    if line[0:2] != '##':
        ls = line.strip()
        lss = line.split()
        if lss[0] == chr:
            if begin_reg <= int(lss[3]) <= end_reg or begin_reg <= int(lss[4]) <= end_reg:
                intervals.append(lss[3] + '\t' + lss[4])

intervals_new = []
for interval in intervals:
    te = interval.split()
    lte = int(te[0]) - begin_reg
    rte = int(te[1]) - begin_reg
    interval = str(lte) + '\t' + str(rte)
    intervals_new.append(interval)

start = '0'

for line in f5:
    if line[0] == '>':
        start = line.strip().split('=')[1]
    else:
        stop = str(int(start) + len(line.strip()) + 1)
        interval = start + '\t' + stop
        intervals_new.append(interval)

stroka = 0
n = 0
for line in f1:
    stroka += 1
    if stroka == 1:
        f3.write(line)
        continue
    out = ''
    for i in line:
        if i == 'A' or i == 'T' or i == 'G' or i == 'C' or i == 'N':
            n += 1
            inside = False
            for interval in intervals_new:
                te = interval.split()
                lte = int(te[0])
                rte = int(te[1])
                if lte < n < rte:
                    inside = True
            #print(n, lte, rte, inside)
            if inside:
                L = i.lower()
            else:
                L = i
            out = out + L
        else:
            out = out + i

    f3.write(out)

f1.close()
f2.close()
f3.close()
