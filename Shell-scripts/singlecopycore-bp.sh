#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_sccg.txt --time 16:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1

# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p

# 4) Constructing the tree with IQTree
module load iqtree
iqtree -s ALIGN-bp -bb 1000 -msub nuclear
