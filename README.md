# Fermentation-bacteria

The entire pipeline of this year-long project is unrolled here in detail.

# November-December 2022 - Preprocessing
## 1) Downloading the whole-genome sequences
The first step consists of downloading whole-genome sequences from the **National Center for Biotechnology Information (NCBI) database** given a list of species and genera of interest. We limit our selection to sequences provided by GenBank only in order to avoid inconsistencies.

## 2) Initialisation
This step takes place right after downloading the genomes from NCBI; the code is run locally. Before proceeding to quality assessment, we must first move all the whole-genome sequences to a single directory, save all the information we have on each genome to an Excel log, and remove unnecessary files. The custom Python script below (*initialiser.py*) was written to do all of this.

### initialiser.py
```python
import pandas as pd
from tqdm import tqdm
import os

#prepare the working directory
bacteria_log = pd.DataFrame()
path = '/Users/thibaultbret/Genomes/'
os.chdir(path)
bacteria = [file for file in sorted(os.listdir(path)) if not file.startswith('.')]

for bacterium in bacteria:
    #select the working directory
    path = '/Users/thibaultbret/Genomes/' + bacterium + '/'
    os.chdir(path)
    #fill the bacteria log dataframe using the "data_summary.tsv" file present in the folder downloaded from NCBI
    info_df = pd.read_table('data_summary.tsv', usecols = [0,2,3,4,5,7,8,9,10,11,12,13,14]).rename(columns = {'Organism Scientific Name':'Organism'}).set_index('Organism')
    bacteria_log = bacteria_log.append(info_df)
    #remove unnecessary files
    os.remove('dataset_catalog.json')
    os.remove('assembly_data_report.jsonl')
    #move all the genome sequences to the bacterium-specific folder
    for elem in tqdm(os.listdir(path)):
        if os.path.isdir(elem):
            os.system('mv ' + elem + '/*.fna ' + path)
            if os.path.exists(elem + '/sequence_report.jsonl'): os.remove(elem + '/sequence_report.jsonl')
            os.rmdir(elem)

#Save the full bacteria log data as an Excel sheet
bacteria_log.index = [idx[0] + ' ' + idx[1] for idx in bacteria_log.index.str.split(' ')]
bacteria_log.to_excel('/Users/thibaultbret/bacteria_log00.xlsx')
```

This script creates a directory named *Genomes* to store all the whole-genome sequences. Unnecessary files are removed. Additionally, the sequences are grouped in subfolders corresponding to their species. The script then fills a **Pandas** data frame with information on the genomes extracted from the "data_summary.tsv" file present in the original folders downloaded from **NCBI**. This includes the organism qualifier, the taxonomy ID, the assembly name, the assembly accession, the annotation level, the genome size, the submission date, the gene count, the corresponding BioProject and BioSample IDs, whether it is part of the lactic or the acetic bacterial family, and information on the isolation source (species, body part and country/region). Finally, the the full Pandas log data frame is exported as an **Excel** file named *bacteria_log.xlsx*

## 3) Quality assessment using [CheckM](https://github.com/Ecogenomics/CheckM)

Once all the genome sequences have been downloaded and their information saved to the *bacteria_log.xlsx* Excel sheet, we can export them to the **Mjolnir** server with the following command:

`scp -r /Users/thibaultbret/Genomes vhp327@mjolnirhead01fl.unicph.domain:/projects/mjolnir1/people/vhp327/`

Now that the genome sequences are present on the server, we can proceed to the quality assessment. We will use **CheckM**, a set of tools already installed on the server that provides robust estimates of genome completeness and contamination. CheckM assumes by default that the genomes consist of contigs/scaffolds in FASTA format, which corresponds to the format of our files.

I wrote a custom script to run a CheckM lineage analysis

*qa.sh*
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


## 4) Quality filtering
Custom script --> qa_filtering.py
**Local**: *Used after importing the "output" directory locally.*

- Reads the quality assessment output files generated by **CheckM**.
- Converts all the qa data into a **Pandas** dataframe that is then converted into an **Excel** sheet --> "total_qa_df.xlsx".

- Reads the "total_qa_df.xlsx" sheet.
- Based on the completeness and the percentage of contamination of each genome, given the predetermined threshold parameters (completeness threshold = 0.9; contamination threshold = 0.05), filters out genomes that are not fit for the analysis. 
- This is done by creating a new directory called "Genomes_filtered" where only the valid whole-genome sequences are present.


## 5) Sorting
Custom script --> acetic_or_lactic.py
**Mjolnir**: *Used after exporting the "Genomes_filtered" directory to the Mjolnir server.*

- Writes a file ("sort_file.txt") that lists all the genomes and indicates whether each genome is part of the Lactic or the Acetic bacteria family.

- Reads the "sort_file.txt" file that assigns every genome to either the Lactic or the Acetic bacteria family.
- Correspondingly moves the genome sequence files in the "Acetic" folder or the "Lactic" folder. 
- That way, lactic and acetic bacteria can be treated separately for later steps of the analysis.


## 6) Further quality filtering
Custom script --> first_tree.py
**Local**: *Used after separating the acetic bacteria from the lactic bacteria and filtering out contaminated and incomplete genomes.*

- Retrieves the percentage of completeness of all genomes from the "bacteria log v2.xlsx" Excel sheet (the updated version of "bacteria_log.xlsx").
- For each species of interest (acetic or lactic), finds the file corresponding to the most complete genome of that species.
- Creates a new directory ("Acetic_unique"/"Lactic_unique") containing only the most complete genome per species.

This directory is then exported to the server so that we can proceed to the following stages with **Anvi'o**.

# Complementary plots

G-C contents:
for file in *.fa ; do gc=$(awk '/^>/ {next;} {gc+=gsub(/[GCgc]/,""); at+=gsub(/[ATat]/,"")} END {print (gc/(gc+at))*100}' $file) ; echo "${file%.fa}:$gc" >> GC_contents.txt ; done

# Anvi'o pangenomics

conda activate anvio-7.1

anvi-gen-genomes-storage -e external-genomes.txt -o STORAGE-GENOMES.db

anvi-pan-genome -g STORAGE-GENOMES.db -n PANGENOME
  
mv STORAGE-GENOMES.db PANGENOME

cd PANGENOME

anvi-compute-gene-cluster-homogeneity -p PANGENOME-PAN.db -g STORAGE-GENOMES.db -o homogeneity_output.txt --store-in-db

anvi-display-pan -p PANGENOME-PAN.db -g STORAGE-GENOMES.db

anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o filtered-concatenated-proteins.fa --max-num-genes-from-each-genome 1 --min-num-genes-from-each-genome 1 --min-num-genomes-gene-cluster-occurs 32 --max-combined-homogeneity-index 0.75

#anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db -o genes.fa#

anvi-gen-phylogenomic-tree -f filtered-concatenated-proteins.fa -o tree.newick

anvi-interactive -p phylogenomic-profile.db -t tree.newick --title "Pangenome tree" --manual


# Single-copy core genes tree

...



# July-August 2023 - New trees based on selected genes and pathways specific to acetic acid bacteria

## General tree
### Preface
- **Contigs database**: Generated by **Anvi'o**, a contigs database is equivalent to a FASTA file with the added advantage of being able to store additional information about the sequences it contains. In this project, each genome sequence was converted to a corresponding contigs database file ending with the suffix "*.CONTIGS.db*".
- **Pfam Accessions**: Pfam is a database of protein families that helps in identifying common motifs in protein families and domains. Each protein family in Pfam is assigned a unique identifier called a "Pfam accession".
- **Hidden Markov Models (HMMs)**: HMMs are probabilistic  models that are commonly used in statistical pattern recognition and classification. In this context, HMMs are generated to represent the characteristics of a given Pfam protein family and are based on patterns in their amino acid sequences. Generating HMMs will later allow us to scan through the genomes and identify matching sequences.

### Steps
1. List pathways/proteins/enzymes of interest and find the corresponding Pfam accessions by browsing the database. This is the only manual step of this pipeline.
*Note: The following code can only be run with **anvio-dev**, a corresponding module exists on the server but I could not make it work so I installed a version of **anvio-dev** locally following this tutorial:* https://anvio.org/install/#5-follow-the-active-development-youre-a-wizard-arry
2. Once everything is properly set up, activate the **anvio-dev** environment with:
```bash
conda activate anvio-dev
 ```
3. Given the list of Pfam accessions of interest, run the **anvi-script-pfam-accessions-to-hmms-directory** command to generate **Hidden Markov Models (HMMs)** corresponding to each Pfam sequence. This will later allow us to scan through any contigs database and detect the presence of matching sequences. An output folder for the HMMs also has to be provided, here we name it *HMM_AAB* (for Hidden Markov Model - Acetic Acid Bacteria).
```bash
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list PF03070 PF08042 PF13360 PF00171 PF00465 PF02317 PF00005 PF02887 PF00923 PF03971 PF00168 PF01161 PF13243 PF00330 PF17327 PF00196 PF00958 PF08240 PF00118  -O HMM_AAB
 ```
*Note: The rest of the code was run from the server in a file called metabolics.sh (shown below).*
![](pfam-general-tree.png)

4. Annotate the genomes with HMMs hits using the **anvi-run-hmms** command.
5. Now we can go through every contigs database (listed in the file *external-genomes.txt*) and retrieve sequences that got a hit for any given HMM using the **anvi-get-sequences-for-hmm-hits** command. The sequences are concatenated so that a phylogenetic tree can be drawn from them (as aligned sequences are a prerequisite). Furthermore, we select on the best hits using the *--return-best-hit* flag. From the Anvi'o manual: "This flag is most appropriate if one wishes to perform phylogenomic analyses, which ensures that for any given protein family, there will be only one gene reported from a given genome". The concatenated best hits are returned in the *pfam-proteins.fa* file.
6. Now a phylogenetic tree can be drawn from the *pfam-proteins.fa* file using the **anvi-gen-phylogenomic-tree** command. However, this will generate a rough tree which is not entirely reliable. Alternatively, a tree can drawn using IQ-TREE simply by using the following command:
```bash
iqtree -s pfam-proteins.fa
```
*Note: However, after trying to run this code, I noticed that the model selection process runs indefinitely (possibly because of the number of sequences the tree is based on). Therefore, I decided to stick to the Anvi'o-generated tree.*


## Individual trees
The idea here is to see if the tree structure would vary depending on the genes we base it on. To test this hypothesis, we build individual phylogenies, still using all 32 genomes, but each being based on a different Pfam identifier.
### 1) Set up HMMs
For each individual Pfam ID, generate a corresponding Hidden Markov Model.
```bash
pfams=("PF03070" "PF08042" "PF13360" "PF00171" "PF00465" "PF02317" "PF00005" "PF02887" "PF00923" "PF03971" "PF00168" "PF01161" "PF13243" "PF00330" "PF17327" "PF00196" "PF00958" "PF08240" "PF00118")
for pfam in "${pfams[@]}"; do anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list $pfam -O HMMs/${pfam}_HMM ; done
```

### 2) Annotating the genome with HMM hits
```bash
for hmm in $(ls HMMs); do for file in *CONTIGS.db; do anvi-run-hmms -c "$file" -H HMMs/$hmm -T 4; done ; done
```

### 3) Estimating metabolism
```bash
cd HMMs
for hmm in $(ls); do anvi-get-sequences-for-hmm-hits --external-genomes ../external-genomes.txt -o  ../PFAM_seqs/${hmm%%_HMM}-proteins.fa --hmm-source $hmm --return-best-hit --get-aa-sequences --concatenate ; done
```
### 4) Building the trees with IQ-TREE
![](iqtree.png)
The resulting trees are shown in the *trees.pdf* file.


