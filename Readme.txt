This program package allows to design the homeolog-specific primes
for gene cloning (by fragments) in allopolyploid bread wheat (Triticum aestivum). 
The package uses the information from IWGSC RefSeq1.0 wheat genome and its annotations,
and primer3 program package (on Linux).

INSTALLATION

Please, download all the files from a GitHub repository to a folder on your computer.

Then download the wheat genome RefSeq1.0 archive as well as genome annotation archives using the links below:

https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/iwgsc_refseqv1.0_all_chromosomes.zip
https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3.zip
https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/iwgsc_refseqv1.1_genes_2017July06.zip
https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_PGSB_annotation_files.zip

Unzip the following files to the folder, containing the files of the repository (the links and the files are in correspondance):

161010_Chinese_Spring_v1.0_pseudomolecules.fasta
iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3
IWGSC_v1.1_HC_20170706.gff3
Triticum_aestivum_V1_PGSB.homeologous_gene_groups.txt

You will need unzip and add exactly the files listed above, nothing more. 
These files were not added to the GitHub repository due to their great volume.

Also, you will need to install primer3. 
On Ubuntu you may use the following command:

sudo apt install primer3


OPERATION

Put the genes IDs to the LIST.txt file.

For example:
TraesCS2B02G383600

Put a single number to the MARGINS.txt file, which will indicate the size in base pairs
of the flanking sequences around the gene of interest to be covered by PCR amplicons

For example:
700


Decide the largest size of the amplicon you wish to get from the following list and run the corresponding script:

bp	Script
1500	RUN1500.py
2000	RUN2000.py
2500	RUN2500.py
3500	RUN3500.py

For example:
python3 RUN2000.py

The primers will be found in the file named 'TraesCS2B02G383600_primers.txt'