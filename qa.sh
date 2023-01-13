module load checkm-genome/1.1.3 
p=/projects/mjolnir1/people/vhp327/
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