# Fermentation-bacteria

The scripts used as part of this project are here described in order:

## initialiser.py
**Local**: *Used after downloading whole-genomes sequences from NCBI.*

- Sets up the working directory.
- Fills a **Pandas** dataframe with information on the genomes using the "data_summary.tsv" file present in the folder downloaded from **NCBI**.
- The dataframe is filled with informations such as the organism qualifier, taxonomic ID, assembly name/accession, size, gene count, BioProject ID, BioSample ID, etc.
- Removes unnecessary files.
- Moves all the genome sequences to the bacterium-specific folder (removes unnecessary nested folders).
- Saves the full bacteria log data as an **Excel** sheet --> "bacteria_log.xlsx"


## [CheckM](https://github.com/Ecogenomics/CheckM)
**Mjolnir**: *Used after exporting the whole-genome sequence files.*

~~~
#!/bin/sh
#SBATCH -c 8 --mem 40G --output=Acetic.xmfa --time 14-0:00
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
~~~

**What does CheckM do?**
- Assesses the quality of microbial genomes recovered from isolates, single cells, and metagenomes.
- Uses **Prodigal** to identify genes in a given sequence.
- Places genome bins into a reference genome tree and relies on that phylogenetic placement to identify sets of genes that should be in the genome.
- Looks for those genes using **Hmmer** to search the genome. 
- The completeness is an estimate of the fraction of genes that are expected to be there which were actually found. 
- The contamination is based on identifying the number of single-copy genes, that should only be there once.
- Returns quality assessment files include the percentage of completeness and the degree of contamination for each given genome.
- The output files are saved in a directory named "output".


## qa_filtering.py
**Local**: *Used after importing the "output" directory locally.*

- Reads the quality assessment output files generated by **CheckM**.
- Converts all the qa data into a **Pandas** dataframe that is then converted into an **Excel** sheet --> "total_qa_df.xlsx".

- Reads the "total_qa_df.xlsx" sheet.
- Based on the completeness and the percentage of contamination of each genome, given the predetermined threshold parameters (completeness threshold = 0.9; contamination threshold = 0.05), filters out genomes that are not fit for the analysis. 
- This is done by creating a new directory called "Genomes_filtered" where only the valid whole-genome sequences are present.


## acetic_or_lactic.py
**Mjolnir**: *Used after exporting the "Genomes_filtered" directory to the Mjolnir server.*

- Writes a file ("sort_file.txt") that lists all the genomes and indicates whether each genome is part of the Lactic or the Acetic bacteria family.

- Reads the "sort_file.txt" file that assigns every genome to either the Lactic or the Acetic bacteria family.
- Correspondingly moves the genome sequence files in the "Acetic" folder or the "Lactic" folder. 
- That way, lactic and acetic bacteria can be treated separately for later steps of the analysis.


## first_tree.py
**Local**: *Used after separating the acetic bacteria from the lactic bacteria and filtering out contaminated and incomplete genomes.*

- Retrieves the percentage of completeness of all genomes from the "bacteria log v2.xlsx" Excel sheet (the updated version of "bacteria_log.xlsx").
- For each species of interest (acetic or lactic), finds the file corresponding to the most complete genome of that species.
- Creates a new directory ("Acetic_unique"/"Lactic_unique") containing only the most complete genome per species.

This directory is then exported to the server so that **Prokka** can be used to annotate all the selected genomes.


## [Prokka](https://github.com/tseemann/prokka)
**Mjolnir**: *Used after exporting the "Acetic_unique"/"Lactic_unique" directory (containing only the most complete genomes) to the server.*

~~~
#!/bin/sh
#SBATCH -c 8 --mem 40G --output=Acetic.xmfa --time 14-0:0
module load prokka/1.14
p=/projects/mjolnir1/people/vhp327/Acetic_unique
cd $p
mkdir Acetic_unique_annotated
for file in *.fna; do tag=${file%.fna}; prokka --prefix "$tag" --outdir Acetic_unique_annotated/"$tag"_prokka "$file"; done
~~~

**What is Prokka?**
- Command line software tool that can be installed on any Unix system. Prokka coordinates a suite of existing software tools to achieve a rich and reliable annotation of genomic bacterial sequences.
- Prokka expects preassembled genomic DNA sequences in FASTA format. Finished sequences without gaps are the ideal input, but it is expected that the typical input will be a set of scaffold sequences produced by de novo assembly software.
- The tools used are: **Prodigal** (Hyatt 2010) to identify coding sequence (CDS), **RNAmmer** (Lagesen et al., 2007) to identify ribosomal RNA genes (rRNA), **Aragorn** (Laslett and Canback, 2004) to identify transfer RNA genes, **SignalP** (Petersen et al., 2011) to identify signal leader peptides and **Infernal** (Kolbe and Eddy, 2011) for non-coding RNA.
- The traditional way to predict what a gene codes for is to compare it with a large database of known sequences, usually at a protein sequence level, and transfer the annotation of the best significant match.
- Prokka uses this method, but in a hierarchical manner, starting with a smaller trustworthy database, moving to medium-sized but domain-specific databases, and finally to curated models of protein families.
- Prokka produces 10 files in the specified output directory, all with a common prefix. The GFF v3 file (.gff) containing sequences and annotations is the one that will be used later in the pipeline.


## [Roary](https://sanger-pathogens.github.io/Roary/)
**Mjolnir**: *Used after annotating the genomes located in the "Acetic_unique"/"Lactic_unique" directory with Prokka.*

Once the genomes have been annotated, they can be aligned using **Roary**.

~~~
#!/bin/sh
#SBATCH -c 8 --mem 40G --output=Acetic.xmfa --time 14-0:00
module load roary/3.13.0
p=/projects/mjolnir1/people/vhp327/Acetic_unique
cd $p/Acetic_unique_annotated
mkdir ../Acetic_unique_roary
for dir in $(ls $p)
do
  cp $dir/*.gff ../../Acetic_unique_roary
done
cd $p/Acetic_unique_roary
roary –e –mafft *.gff
~~~

**What does Roary do?**
- From the ".gff" files given by **Prokka**, converts coding sequences into protein sequences.
- Cluster these protein sequences by several methods.
- Further refines clusters into orthologous genes.
- For each sample, determines if gene is present/absent: produces "gene_presence_absence.csv".
- Uses this gene p/a information to build a tree, using **FastTree**: produces "accessory_binary_genes.fa.newick".
- Overall, calculates number of genes that are shared, and unique: produces "summary_statistics.txt".
- Aligns the core genes using **Mafft** (if option used, as above) for downstream analyses: produces "core_gene_alignment.aln".


## [Peppan](https://github.com/zheminzhou/PEPPAN)
**Mjolnir**: *Alternative to Roary.*

**PEPPAN (Phylogeny Enhanced Pipeline for PAN-genome)**: Pipeline that can construct a pan-genome from thousands of genetically diversified bacterial genomes. PEPPAN implements a combination of tree- and synteny-based approaches to identify and exclude paralogous genes, as well as similarity-based gene predictions that support consistent annotations of genes and pseudogenes in individual genomes.

PEPPAN's workflow consists of the following five successive groups of operations:
- Identifying representative gene sequences: The inputs for PEPPAN consist of GFF3 formatted genome assemblies. To reduce the number of genes used in downstream analyses, PEPPAN iteratively clusters genes using Linclust (Steinegger and Söding 2017), resulting in a single representative gene sequence per 90% nucleotide homology cluster.
- Identifying gene candidates: Each representative gene is aligned to all genomes using both BLASTN (Altschul et al. 1990), which accurately locates short inserts and deletions (indels), and DIAMOND (Buchfink et al. 2015), which generates amino acid alignments and has greater sensitivity with divergent sequences than BLASTN. Alignments are rescored, and all sequences with homology ≥50% across ≥50% of the representative sequence are clustered in a neighbour-joining tree using RapidNJ (Simonsen et al. 2011).
- Identifying clusters of orthologous genes: PEPPAN identifies putative orthologs by calculating a paralogous score for each branch in a gene cluster tree based on the ratio of the pairwise genetic distances of candidate genes within each cluster to the average genetic distances of their host genomes. Using average genetic distances avoids potential errors that can be introduced by using a “species” tree to reconcile individual gene cluster trees. Branches with a paralogous score of greater than one are iteratively pruned until none remain. The remaining monophyletic subtrees are treated as putative orthologs.
- The genomic locations of multiple putative orthologs may overlap in some genomes owing to either inconsistent genome annotations or a failure to cluster divergent orthologous sequences in the first stage. These conflicts are resolved by retaining the ortholog with the greatest information score and eliminating all other gene candidates for that region.
- The remaining gene candidates from each genome are ordered according to their genomic coordinates, and the final set of orthologous genes is identified based on synteny.
- Pangenome outputs: Each gene candidate in each genome is categorised as either an intact coding sequence (CDS) or a pseudogene, depending on the size of the aligned reading frame relative to its representative gene. It is also possible to predict pseudogenes that are disrupted in all genomes by importing their intact analog into PEPPAN as an external representative gene. Finally, the evaluations of all genes, as well as their genomic coordinates and orthologous group, are output in GFF3 format. The extent of the regions that match their representative genes is saved in FASTA format.
- Pangenome analysis: A separate tool, PEPPAN_parser, generates analyses of the estimated pangenome based on the GFF3 outputs from PEPPAN. Similar to Roary (Page et al. 2015) and PIRATE (Bayliss et al. 2019), these include rarefaction curves, gene presence matrices, and gene presence trees. In addition, PEPPAN_parser can also calculate a core genome tree based on allelic differences of genes that are conserved in most genomes. These core genome trees can scale to 10,000s genomes and provide the basis for all core genome MLST schemes in EnteroBase (Supplemental Text 3; Zhou et al. 2020).


## [FastTree](http://www.microbesonline.org/fasttree/)
**Mjolnir**: *Used after constructing the pan-genome alignment with Roary.*

Once the genomes have been aligned, we can construct the tree using **FastTree**.

~~~
#!/bin/sh
#SBATCH -c 8 --mem 40G --output=Acetic.xmfa --time 14-0:00
module load roary/3.13.0
FastTree -nt -gtr output_with_alignment/core_gene_alignment.aln \
    > tree.newick
~~~

**What does FastTree do?**
- Infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequence.
- Generates the tree file in **Newick format** ("tree.newick") based on the core alignment file ("core_gene_alignment.aln").
- Faster than **PhyML 3.0** and **RAxML 7** for large alignments.
- (-gtr) --> generalized time-reversible model.


## Visualisation
### With Phandango
- Go to http://phandango.net
- Drag and drop the "tree.newick" and "gene_presence_absence.csv" files onto the landing page.
- View the tree of samples and their core and pan genomes

### With Roary plots
~~~
python3 roary_plots.py tree.newick gene_presence_absence.csv
~~~
