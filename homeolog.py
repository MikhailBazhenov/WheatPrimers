f2 = open('TARGET_GENE_ID.txt')
f3 = open('offtarget_genes_IDs.txt', 'w')
f4 = open('all_gene_IDs.txt', 'w')

for line in f2:
    lss = line.strip().split()
    gene = lss[0].replace('02G', '01G')
    wider = lss[1]
    f1 = open('Triticum_aestivum_V1_PGSB.homeologous_gene_groups.txt', 'r')
    f4.write(gene.replace('01G', '02G') + '\t' + wider + '\n')
    for line in f1:
        lss = line.strip().split()
        cat = lss[1].replace('"', '').replace(',', '\t') + '\t' + lss[2].replace('"', '').replace(',', '\t') + '\t' + \
              lss[3].replace('"', '').replace(',', '\t')
        list = cat.split()
        if gene in list:
            for l in list:
                if l != gene:
                    f3.write(l.replace('01G', '02G') + '\t' + wider +'\n')
                    f4.write(l.replace('01G', '02G') + '\t' + wider + '\n')