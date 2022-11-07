#!/usr/bin/env python3
# Finds gene positions in gff3 file by gene ID

f1 = open('IWGSC_v1.1_HC_20170706.gff3', 'r')
f2 = open('TARGET_GENE_ID.txt', 'r')
f3 = open('genes_positions.txt', 'w')

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
                start = li[3]
                end = li[4]
                strand = li[6]
                res = id + '\t' + chrom + '\t' + start + '\t' + end + '\t' + strand
                # print(res)
                f3.write(res + '\n')
f1.close()
f2.close()
f3.close()

# Finds gene positions in gff3 file by gene ID


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
        if nos < ll / 60:
            nos += 1
        for k in range(nos):
            if k < nos:
                comps = comps + comp[k * 60: k * 60 + 60] + '\n'
            if k == nos:
                comps = comps + comp[k * 60:]
    else:
        comps = comp
    out = comps
    return out


f3 = open('IWGSC_v1.1_HC_20170706.gff3', 'r')
f4 = open('genes_positions.txt', 'r')

fasta = '161010_Chinese_Spring_v1.0_pseudomolecules.fasta'
fai = fasta + '.fai'

genes = []
found = False

for line in f4:
    ln = line.strip().split()
    genes.append(ln[0])

id = '0'
chromseq = ''
geneseq = ''
exonseq = ''
cdsseq = ''
strand = '+'
for line in f3:
    if line[0] == '#':
        if found:
            f5 = open(id + '.fasta', 'w')
            f5.write(chromseq)
            f5.write(geneseq)
            if strand == '+':
              exonseq = revcomp(revcomp(exonseq))
            else:
                exonseq = revcomp(exonseq)
            f5.write('>exons ' + id + '\n' + exonseq)
            if strand == '+':
                cdsseq = revcomp(revcomp(cdsseq))
            else:
                cdsseq = revcomp(cdsseq)
            f5.write('>CDS ' + id + '\n' + cdsseq)
            geneseq = ''
            exonseq = ''
            cdsseq = ''
            f5.close()
    else:
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
                seq = li[0]
                start = int(li[3])
                end = int(li[4])
                strand = li[6]
                f1 = open(fai, 'r')
                f2 = open(fasta, 'r')
                res = ''
                for line in f1:
                    ln = line.strip()
                    li = ln.split()
                    if li[0] == seq:
                        p = int(li[2]) + start + int(start / int(li[3])) * (int(li[4]) - int(li[3])) - 1
                        f2.seek(p)
                        while len(res) < end - start + 1:
                            res = res + f2.read(60).replace('\n', '')
                        res = res[0:end - start + 1]
                        if strand == '-':
                            res = revcomp(res)
                        if strand == '+':
                            res = revcomp(revcomp(res))
                geneseq = '>' + id + '\n' + res
                f1.close()
                f2.close()

                start = start - wider
                end = end + wider
                f1 = open(fai, 'r')
                f2 = open(fasta, 'r')
                res = ''
                for line in f1:
                    ln = line.strip()
                    li = ln.split()
                    if li[0] == seq:
                        p = int(li[2]) + start + int(start / int(li[3])) * (int(li[4]) - int(li[3])) - 1
                        f2.seek(p)
                        while len(res) < end - start + 1:
                            res = res + f2.read(60).replace('\n', '')
                        res = res[0:end - start + 1]
                        if strand == '-':
                            res = revcomp(res)
                        if strand == '+':
                            res = revcomp(revcomp(res))
                chromseq = '>' + seq + '\n' + res
                f1.close()
                f2.close()

        if li[2] == 'exon':
            fid = li[8][li[8].find('ID=')+3:li[8].find(';')]
            id = fid[0:fid.find('.')]
            rna = int(fid.replace(':', '.').split('.')[1])
            for i in genes:
                if id == i:
                    found = True
                    break
                else:
                    found = False
            if found:
                seq = li[0]
                start = int(li[3])
                end = int(li[4])
                strand = li[6]
                f1 = open(fai, 'r')
                f2 = open(fasta, 'r')
                res = ''
                for line in f1:
                    ln = line.strip()
                    li = ln.split()
                    if li[0] == seq:
                        p = int(li[2]) + start + int(start / int(li[3])) * (int(li[4]) - int(li[3])) - 1
                        f2.seek(p)
                        while len(res) < end - start + 1:
                            res = res + f2.read(60).replace('\n', '')
                        res = res[0:end - start + 1]
                if rna == 1:
                    exonseq = exonseq + res
                f1.close()
                f2.close()
        if li[2] == 'CDS':
            fid = li[8][li[8].find('ID=') + 3:li[8].find(';')]
            id = fid[0:fid.find('.')]
            rna = int(fid.replace(':', '.').split('.')[1])
            for i in genes:
                if id == i:
                    found = True
                    break
                else:
                    found = False
            if found:
                seq = li[0]
                start = int(li[3])
                end = int(li[4])
                strand = li[6]
                f1 = open(fai, 'r')
                f2 = open(fasta, 'r')
                res = ''
                for line in f1:
                    ln = line.strip()
                    li = ln.split()
                    if li[0] == seq:
                        p = int(li[2]) + start + int(start / int(li[3])) * (int(li[4]) - int(li[3])) - 1
                        f2.seek(p)
                        while len(res) < end - start + 1:
                            res = res + f2.read(60).replace('\n', '')
                        res = res[0:end - start + 1]
                if rna == 1:
                    cdsseq = cdsseq + res
                f1.close()
                f2.close()
f3.close()
f4.close()

