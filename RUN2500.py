import subprocess
f3 = open('LIST.txt', 'r')
f4 = open('MARGINS.txt', 'r')
margins = f4.read().strip()

for item in f3:
    print('Starting primer search for ', item.strip())
    f5 = open('TARGET_GENE_ID.txt', 'w')
    f5.write(item.strip() + '\t' + margins)
    f5.close()
    print('Search for homeologs')
    subprocess.run(['python3', 'homeolog.py'])
    print('Search for genes coordinates and sequences')
    subprocess.run(['python3', 'find_target_gene_by_ID.py'])
    subprocess.run(['python3', 'find_offtarget_genes_by_ID.py'])
    subprocess.run(['python3', 'find_all_genes_by_ID.py'])
    subprocess.run(['python3', 'find_gene+features.py'])
    
    print('Search for SSRs in target sequence')
    result = subprocess.run(['python3', 'microsat.py', 'target_sequence.txt', '10'], capture_output=True,
                            encoding='utf-8')
    f1 = open('SSRs.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    print('Masking TEs and SSRs')
    subprocess.run(['python3', 'masker_new.py'])
    
    print('Making primers for 500-1000 bp PCR products')
    result = subprocess.run(['python3', 'taskmaker.py', 'target_sequence_masked.txt', '50', '500-1000', 'offtarget_sequences.txt'],
                            capture_output=True, encoding='utf-8')
    out = result.stdout
    result = subprocess.run(['primer3_core'], input=out, capture_output=True, encoding='utf-8')
    f1 = open('primer3_out.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    result = subprocess.run(['python3', 'out-format.py', 'primer3_out.txt'], capture_output=True, encoding='utf-8')
    f1 = open('primers_500.txt', 'w')
    f1.write(result.stdout)
    f1.close()
     
    print('Making primers for 1000-1500 bp PCR products')
    result = subprocess.run(
        ['python3', 'taskmaker.py', 'target_sequence_masked.txt', '50', '1000-1500', 'offtarget_sequences.txt'],
        capture_output=True, encoding='utf-8')
    out = result.stdout
    result = subprocess.run(['primer3_core'], input=out, capture_output=True, encoding='utf-8')
    f1 = open('primer3_out.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    result = subprocess.run(['python3', 'out-format.py', 'primer3_out.txt'], capture_output=True, encoding='utf-8')
    f1 = open('primers_1000.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    print('Making primers for 1500-2000 bp PCR products')
    result = subprocess.run(
        ['python3', 'taskmaker.py', 'target_sequence_masked.txt', '50', '1500-2000', 'offtarget_sequences.txt'],
        capture_output=True, encoding='utf-8')
    out = result.stdout
    result = subprocess.run(['primer3_core'], input=out, capture_output=True, encoding='utf-8')
    f1 = open('primer3_out.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    result = subprocess.run(['python3', 'out-format.py', 'primer3_out.txt'], capture_output=True, encoding='utf-8')
    f1 = open('primers_1500.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    print('Making primers for 2000-2500 bp PCR products')
    result = subprocess.run(
        ['python3', 'taskmaker.py', 'target_sequence_masked.txt', '50', '2000-2500', 'offtarget_sequences.txt'],
        capture_output=True, encoding='utf-8')
    out = result.stdout
    result = subprocess.run(['primer3_core'], input=out, capture_output=True, encoding='utf-8')
    f1 = open('primer3_out.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    result = subprocess.run(['python3', 'out-format.py', 'primer3_out.txt'], capture_output=True, encoding='utf-8')
    f1 = open('primers_2000.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    result = subprocess.run(['cat', 'primers_500.txt', 'primers_1000.txt', 'primers_1500.txt', 'primers_2000.txt'], 
                             capture_output=True, encoding='utf-8')
    f1 = open('primers.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    
    f2 = open('target_position.txt', 'r')
    gene = f2.read().strip().split()
    print('Finding optimal combination of primers')
    result = subprocess.run(
        ['python3', 'primer-sorter.py', 'primers.txt', '500-' + str(int(gene[3]) - int(gene[2]) - 500),
         '500'], capture_output=True, encoding='utf-8')
    f1 = open(gene[0] + '_primers.txt', 'w')
    f1.write(result.stdout)
    f1.close()
    print('Completed')
    print('')

#cd /home/linux/TaSCL14
#python3 masker.py chr4B.txt chr4B_RM.txt
#python3 microsat.py chr4B.txt 10 > SSR_chr4B.txt
#cat chr4A-4D.txt SSR_chr4B.txt > Unwanted.txt

#python3 taskmaker.py masked_chr4B.txt 50 500-1000 Unwanted.txt > primer3_task.txt
#primer3_core < primer3_task.txt > primer3_out.txt
#python3 out-format.py primer3_out.txt > primers_chr4B_500.txt

#python3 taskmaker.py masked_chr4B.txt 50 1000-1500 Unwanted.txt > primer3_task.txt
#primer3_core < primer3_task.txt > primer3_out.txt
#python3 out-format.py primer3_out.txt > primers_chr4B_1000.txt

#cat primers_chr4B_500.txt primers_chr4B_1000.txt > primers_chr4B.txt
#python3 primer-sorter.py primers_chr4B.txt 200-3210 500 > best_primers_chr4B.txt
