#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_meta.txt --time 12:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1

# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p

# 2) Annotating the genome with KOfam hits
for file in $(ls *.db); do anvi-run-hmms -c $file  -H HMM_AAB ; done

# 3) Estimating metabolism
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                -o pfam-proteins.fa \
                                --hmm-source HMM_AAB \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate

# Compute the phylogenomic tree
anvi-gen-phylogenomic-tree -f pfam-proteins.fa \
                           -o pfam-phylogenomic-tree.txt

# Optional: Compute the phylogenomic tree with IQ-TREE
module load iqtree
iqtree -s pfam-proteins.fa