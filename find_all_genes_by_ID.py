#!/usr/bin/env python3
# Finds gene positions in gff3 file by gene ID

f1 = open('IWGSC_v1.1_HC_20170706.gff3', 'r')
f2 = open('all_gene_IDs.txt', 'r')
f3 = open('all_positions.txt', 'w')

genes = []
found = False
wider = 500

for line in f2:
    ln = line.strip().split()
    wider = int(ln[1]) + 500
    genes.append(ln[0].split('.')[0])

for line in f1:
    if line[0] != '#':
        ln = line.strip()
        li = ln.split()
        if li[2] == 'gene':
            id = li[8][li[8].find('ID=')+3:li[8].find(';')]
            for i in genes:
                if id == i:
                    found = True
                    break
                else:
                    found = False
            if found:
                chrom = li[0]
                start = str(int(li[3]) - wider)
                end = str(int(li[4]) + wider)
                strand = li[6]
                res = id + '\t' + chrom + '\t' + start + '\t' + end + '\t' + strand
                print(res)
                f3.write(res + '\n')
f1.close()
f2.close()
f3.close()

# Retrieves gene sequences from genome file by gene positions


def revcomp(q):
    body = q  # тело записи
    body = body.replace('\n', '')  # удаляем переводы строки
    rev = body[len(body):0:-1] + body[0]
    trt = rev.maketrans('ACGTMRWSYKVHDBN', 'TGCAKYWSRMBDHVN')
    comp = rev.translate(trt)
    comps = ''
    ll = len(comp)
    if ll > 60:
        nos = int(ll / 60)
        if round(ll / 60) < ll / 60:
            nos = round(ll / 60) + 1
        for k in range(nos):
            if k < nos:
                comps = comps + comp[k * 60: k * 60 + 60] + '\n'
            if k == nos:
                comps = comps + comp[k * 60:]
    out = comps
    return out


fasta = '161010_Chinese_Spring_v1.0_pseudomolecules.fasta'
fai = fasta + '.fai'

f4 = open('all_positions.txt', 'r')
f5 = open(genes[0]+'_homeologs.txt', 'w')

for gene in f4:
    gs = gene.strip().split()
    seq = gs[1]
    start = int(gs[2])
    end = int(gs[3])
    out = gs[0] + '.fasta'
    f1 = open(fai, 'r')
    f2 = open(fasta, 'r')
    res = ''
    for line in f1:
        ln = line.strip()
        li = ln.split()
        if li[0] == seq:
            p = int(li[2]) + start + int(start / int(li[3])) * (int(li[4]) - int(li[3])) - 1
            end1 = end + int((end - start) / int(li[3])) * (int(li[4]) - int(li[3])) - 1
            f2.seek(p)
            res = f2.read(end1 - start)
            if gs[4] == '-':
                res = revcomp(res)
            if gs[4] == '+':
                res = revcomp(revcomp(res))
    res = '>' + gs[1] + ' ' + gs[0] + '\n' + res
    f5.write(res)
    f1.close()
    f2.close()

f4.close()
f5.close()

