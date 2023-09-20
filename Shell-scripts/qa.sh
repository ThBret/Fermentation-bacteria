#!/bin/sh
#SBATCH -c 8 --mem 4G --output=qa_output.txt --time 16:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load checkm-genome/1.1.3 
p=/projects/mjolnir1/people/vhp327/new/
cd $p
mkdir output
for i in $(ls $p)
do if [ -d $i ]
    then 
        checkm lineage_wf -x fna $i ${i}_output 
        cd ${i}_output
        checkm qa lineage.ms ./ -f qa_${i}.txt
        mv qa_${i}.txt ../output
        cd $p
    fi
done


