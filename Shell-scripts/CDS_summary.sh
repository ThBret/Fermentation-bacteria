#!/bin/sh
#SBATCH -c 8 --mem 1G --output=CDS_output.txt --time 12:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL

# 1) Set the directory and output file name
p=/projects/mjolnir1/people/vhp327/Acetic/Acetic_annotated
cd $p
output_file="CDS_summary.txt"

# 2) Loop over all .tsv files in the directory
for file in $(ls *.tsv); do
  # Compute the sum of all elements in the third column
  total_sum=$(awk '{sum+=$3} END {print sum}' "$file")
  # Compute the sum of all elements in the third column when the second column is 'CDS'
  cds_sum=$(awk '$2=="CDS" {sum+=$3} END {print sum}' "$file")
  # Get the filename without the directory path or extension
  filename=$(basename -- "$file")
  filename="${filename%.*}"
  # Append the filename and sums to the output file
  echo "$filename,$total_sum,$cds_sum" >> "$output_file"
done

