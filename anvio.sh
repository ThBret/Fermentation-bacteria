#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_anvio.txt --time 12:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1

for i in $(ls *.fna); do
    mv $i "${i%.*}".fa
done

for i in $(ls *.fa); do
    anvi-script-reformat-fasta $i -l 0 --simplify-names --seq-type NT -o "${i}.fixed"
    mv ${i}.fixed $i
    if [ -f $i.db ]; then
        echo "$i.db already exists."
    else
         anvi-gen-contigs-database -f $i -o ${i%.*}_CONTIGS.db -T 4
    fi
done

# Make an external genomes database
for i in $(ls *fa); do
    f=${i%%.*}
    echo -e ${f}"\t"${f}_CONTIGS.db
done | grep -v 'FIXED' > external-genomes.txt

