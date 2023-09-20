#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_sccg.txt --time 16:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1

# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p

# 2) Run HMMs for 71 Bacteria copy core genes
for file in $(ls *.db); do anvi-run-hmms -c $file -T 4 -I Bacteria_71 ; done

# 3) Estimating metabolism
anvi-get-sequences-for-hmm-hits -e external-genomes.txt --hmm-sources Bacteria_71 --get-aa-sequences --concatenate-genes -o ALIGN --return-best-hit

# 4) Constructing the tree with IQTree
module load iqtree
iqtree -s ALIGN
