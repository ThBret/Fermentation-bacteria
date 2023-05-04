# Anvio pangenomics:
/!\ These commands are run locally after running the **anvio.sh** script on the server (which generates contigs databases for every genome and the 'external_genomes.txt' file) and after downloading the folder containing those files from the server. /!\

**1.** Activate the Anvi'o Conda environment.
~~~
conda activate anvio-7.1
~~~

**2.** Generate a genomes storage using **anvi-gen-genomes-storage**. The genomes storage is generated from contigs databases corresponding to every genome sequence which are accessed using the 'external_genomes.txt' file. The genomes storage is a needed input file for later steps.
~~~
anvi-gen-genomes-storage -e external_genomes.txt -o STORAGE-GENOMES.db 
~~~

**3.** With the genomes storage ready, we can use the program **anvi-pan-genome** to run the actual pangenomic analysis.
~~~
anvi-pan-genome -g STORAGE-GENOMES.db -n PANGENOME
~~~

**4.** Move the genomes storage to the newly created Pangenome folder and change the working directory to this same folder. The rest of the analysis will be performed within that directory.
~~~
mv STORAGE-GENOMES.db PANGENOME
cd PANGENOME
~~~

**5.** After the analysis is done, we can use the program **anvi-display-pan** to display your results.
~~~
anvi-display-pan -p PANGENOME-PAN.db -g STORAGE-GENOMES.db
~~~

*Optional*: Get a FASTA file with aligned and concatenated amino acid sequences corresponding to all gene clusters found in the pangenome. This is optional since we cannot build the phylogeny using all gene clusters (as the alignment size would be too enormous).
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o concatenated-proteins.fa --max-num-genes-from-each-genome 1
~~~

*Optional*: Get the same FASTA file with non-concatenated amino acids.
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db -o genes.fa
~~~
  
**6.** Compute both the geometric homogeneity and functional homogeneity for the gene clusters in a pangenome database and add this information to the database. Since the phylogenemic inference cannot be performed on the entire pangenome, we will instead only use gene clusters with significant variation (combined homogeneity < 0.75).
~~~
anvi-compute-gene-cluster-homogeneity -p PANGENOME-PAN.db -g STORAGE-GENOMES.db -o homogeneity_output.txt --store-in-db
~~~

**7.** Get a FASTA file with aligned and concatenated amino acid sequences corresponding to the selected gene clusters. This will be used to perform the phylogenomic analysis. We set the *--max-combined-homogeneity-index* to 0.75 to limit our selection to highly variable gene clusters (as they will have a bigger impact on the phylogeny than gene clusters with low variability). We also set the *genomes-gene-cluster-occurs* parameter to 32 as we have 32 genomes in the analysis and we want the gene clusters to be present in every genome.
~~~
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o filtered-concatenated-proteins.fa --max-num-genes-from-each-genome 1 --min-num-genes-from-each-genome 1 --min-num-genomes-gene-cluster-occurs 32 --max-combined-homogeneity-index 0.75
~~~

**8.** Perform phylogenetic inference based on the previously generated FASTA file containing aligned and concatenated gene clusters of interest.
~~~
anvi-gen-phylogenomic-tree -f filtered-concatenated-proteins.fa -o tree.newick
~~~

**9.** Display the phylogenetic tree in the Anvi'o interactive interface.
~~~
anvi-interactive -p phylogenomic-profile.db -t tree.newick --title "Pangenome tree" --manual
~~~
