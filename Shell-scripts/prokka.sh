#!/bin/sh
#SBATCH -c 8 --mem 1G --output=prokka_output.txt --time 12:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load prokka/1.14

p=/projects/mjolnir1/people/vhp327/Acetic
cd $p
mkdir Acetic_annotated

# 1) Prokka annotation
for file in *.fna; do tag=${file%.fna}; prokka --prefix "$tag" --outdir Acetic_annotated --force $file ; done

# 2) Delete unnecessary files
for file in ${p}/Acetic_annotated/*; do
    find ${file} ! -name '*.tsv' -type f -delete
done
