#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_iqtree.txt --time 24:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load iqtree

#construct iq-trees
for seq in $(ls *.fa); do iqtree -s $seq -bb 1000 ; done
