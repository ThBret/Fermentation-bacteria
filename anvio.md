# Anvio pangenomics:
**1.** Generate a genomes storage using **anvi-gen-genomes-storage**. The genomes storage is generated from the 'external_genomes.txt' file and is a needed input file for later steps.
~~~
anvi-gen-genomes-storage -e external_genomes.txt -o STORAGE-GENOMES.db 
~~~

**2.** With the genomes storage ready, we can use the program **anvi-pan-genome** to run the actual pangenomic analysis.
~~~
anvi-pan-genome -g STORAGE-GENOMES.db -n PANGENOME
~~~

**3.** After the analysis is done, we can use the program **anvi-display-pan** to display your results.
~~~
anvi-display-pan -p PANGENOME-PAN.db -g STORAGE-GENOMES.db
~~~

**4.** Get a FASTA file with aligned and concatenated amino acid sequences corresponding to gene clusters found in the pangenome. This will be used to perform a phylogenomic analysis.
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o concatenated-proteins.fa --max-num-genes-from-each-genome 1
~~~

*Optional*: Get the same FASTA file with non-concatenated amino acids.
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db -o genes.fa
~~~
  
**5.** Compute both the geometric homogeneity and functional homogeneity for the gene clusters in a pangenome database and add this information to the database. 
This is done because the phylogenemic inference cannot be performed on the entire pangenome, instead we will only use gene clusters with significant variation (combined homogeneity < 0.75).
~~~
anvi-compute-gene-cluster-homogeneity -p PANGENOME-PAN.db -g STORAGE-GENOMES.db -o homogeneity_output.txt --store-in-db
~~~

**6.** Get a FASTA file with aligned and concatenated amino acid sequences corresponding to the selected gene clusters.
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o filtered-concatenated-proteins.fa --max-num-genes-from-each-genome 1 --max-combined-homogeneity-index 0.75
~~~

**7.** Perform phylogenetic inference based on the previously generated FASTA file containing aligned and concatenated gene clusters of interest.
~~~
anvi-gen-phylogenomic-tree -f filtered-concatenated-proteins.fa -o tree.newick
~~~

**8.** Display the phylogenetic tree in the Anvi'o interactive interface.
~~~
anvi-interactive -p phylogenomic-profile.db -t tree.newick --title "Pangenome tree" --manual
~~~
