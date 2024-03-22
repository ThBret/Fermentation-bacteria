# Fermentation-bacteria

> The entire pipeline of this year-long project is unrolled here in detail.

# November-December 2022 - Preprocessing
## 1) Downloading the whole-genome sequences
The first step consists of downloading whole-genome sequences from the **National Center for Biotechnology Information (NCBI) database** given a list of species and genera of interest. We limit our selection to sequences provided by GenBank only in order to avoid inconsistencies.

## 2) Initialisation
This step takes place right after downloading the genomes from NCBI; the code is run locally. Before proceeding to quality assessment, we must first move all the whole-genome sequences to a single directory, save all the information we have on each genome to an Excel log, and remove unnecessary files. The custom Python script below (*initialiser.py*) was written to do all of this.

<details>
  <summary><b>initialiser.py</b> <i>(see code)</i></summary>

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

</details>

This script creates a directory named *Genomes* to store all the whole-genome sequences. Unnecessary files are removed. Additionally, the sequences are grouped in subfolders corresponding to their species. The script then fills a **Pandas** data frame with information on the genomes extracted from the "data_summary.tsv" file present in the original folders downloaded from **NCBI**. This includes the organism qualifier, the taxonomy ID, the assembly name, the assembly accession, the annotation level, the genome size, the submission date, the gene count, the corresponding BioProject and BioSample IDs, whether it is part of the lactic or the acetic bacterial family, and information on the isolation source (species, body part and country/region). Finally, the the full Pandas log data frame is exported as an **Excel** file named *bacteria_log.xlsx*

## 3) Quality assessment using [CheckM](https://github.com/Ecogenomics/CheckM)
Once all the genome sequences have been downloaded and their information saved to the *bacteria_log.xlsx* Excel sheet, we can export them to the **Mjolnir** server with the following command:

```bash
scp -r /Users/thibaultbret/Genomes vhp327@mjolnirhead01fl.unicph.domain:/projects/mjolnir1/people/vhp327/
```

Now that the genome sequences are present on the server, we can proceed to the quality assessment. We will use **CheckM**, a set of tools already installed on the server that provides robust estimates of genome completeness and contamination. CheckM assumes by default that the genomes consist of contigs/scaffolds in FASTA format, which corresponds to the format of our files. I wrote a custom Slurm script to run a CheckM lineage analysis on the server: *qa.sh*.

<details>
  <summary><b>qa.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=qa_output.txt --time 16:00:00
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
```

</details>

**What does CheckM do?**
- Assesses the quality of microbial genomes recovered from isolates, single cells, and metagenomes.
- Uses **Prodigal** to identify genes in a given sequence.
- Places genome bins into a reference genome tree and relies on that phylogenetic placement to identify sets of genes that should be in the genome.
- Looks for those genes using **Hmmer** to search the genome. 
- The completeness is an estimate of the fraction of genes that are expected to be there which were actually found. 
- The contamination is based on identifying the number of single-copy genes, that should only be there once.
- Returns quality assessment files include the percentage of completeness and the degree of contamination for each given genome.
- The output files are saved in a directory named *output*.

## 4) Quality filtering
After running the CheckM analysis on the server, output files are returned and saved in a directory named *output*. This folder is imported locally to filter the list of genomes based on the results of the completeness and contamination assessments. The filtering is performed by another custom Python script: *qa_filtering.py*.

<details>
  <summary><b>qa_filtering.py</b> <i>(see code)</i></summary>

```python
import pandas as pd
from tqdm import tqdm
import os

def qa_reader(path = '/Users/thibaultbret/output'):
    '''
    Reads the quality assessment output files generated by CheckM and converts all the qa data into a Pandas dataframe
    that is then converted into an Excel sheet.
    '''
    #find the qa files
    total_qa_df = pd.DataFrame()
    os.chdir(path)
    qa_files = [file for file in sorted(os.listdir(path)) if not file.startswith('.')]
    #gather the data
    for qa_file in tqdm(qa_files):
        names = ['Bin Id','Marker lineage','#genomes','#markers','#marker sets','0','1','2','3','4','5+','Completeness','Contamination','Strain heterogeneity']
        data = []
        with open(qa_file) as fl:
            for line in fl.readlines()[2:]:
                if not line.startswith('-'):
                    data.append([l for l in line.split(' ') if not l in ['','\n'] and not l.startswith('(')])
        qa_df = pd.DataFrame(data, columns = names)
        dtypes_dict = {'Bin Id':str, 'Marker lineage':str, '#genomes':int, '#markers':int, '#marker sets':int, '0':int, '1':int, '2':int, '3':int, '4':int, '5+':int, 'Completeness':float, 'Contamination':float, 'Strain heterogeneity':float}
        qa_df = qa_df.astype(dtypes_dict)
        qa_df['Organism'] = qa_file.replace('qa_','').replace('.txt','')
        total_qa_df = total_qa_df.append(qa_df, ignore_index=True)
    #export excel sheet 
    total_qa_df.to_excel('/Users/thibaultbret/total_qa_df.xlsx', index = False)

def qa_filterer(path = '/Users/thibaultbret', total_qa_xlsx = '/Users/thibaultbret/total_qa_df0.xlsx', completeness_th = 0.9, contamination_th = 0.05):
    '''
    Reads the "total_qa_df.xlsx" Excel sheet and, based on the completeness and the percentage of contamination of each genome, given
    the predetermined threshold parameters, filters out genomes that are not fit for the analysis.
    '''
    #Extract the excel qa sheet and filter out uninformative genomes
    total_qa_df = pd.read_excel(total_qa_xlsx)
    completeness_mask = total_qa_df['Completeness'] > completeness_th * 100
    contamination_mask = total_qa_df['Contamination'] < contamination_th * 100
    filtered_qa_df = total_qa_df[~(completeness_mask & contamination_mask)]
    genomes_to_remove = list(filtered_qa_df['Bin Id'])
    #create a copy of the Genome folder without the uninformative genomes
    os.chdir(path)
    os.system('cp -r Genomes Genomes_filtered')
    new_path = path + '/Genomes_filtered'
    os.chdir(new_path)
    genome_folders = [folder for folder in sorted(os.listdir(new_path)) if not folder.startswith('.')]
    #move all .fna files to the newly created directory
    for genome_folder in genome_folders:
        os.system('mv ' + genome_folder + '/*.fna ' + new_path)
    #remove all nested directories
    for genome_folder in genome_folders:
        os.system('rm -r ' + genome_folder)
    #delete the .fna files that correspond to genomes that should be filtered out
    fna_files = [file for file in sorted(os.listdir(new_path)) if not file.startswith('.')]
    for fna_file in fna_files:
        if fna_file in [gen + '.fna' for gen in genomes_to_remove]:
            os.remove(fna_file)

if __name__ == "__main__":
    qa_reader()
    qa_filterer()
```

</details>

This script reads the quality assessment output files generated by **CheckM** and converts all the qa data into a **Pandas** data frame that is then exported to an **Excel** sheet named *total_qa_df.xlsx*. The second part of the script then reads the *total_qa_df.xlsx* file. Based on the completeness and the percentage of contamination of each genome, and given predetermined threshold parameters (**completeness threshold = 0.9**; **contamination threshold = 0.05**), genomes that do not meet the threshold are filtered out and excluded from the rest of the analysis. This is done by creating a new directory called *Genomes_filtered* where only the whole-genome sequences that passed the completeness and contaminations tests are present.

## 5) Sorting
Once created, the *Genomes_filtered* directory, containing only whole-genome sequences that passed the completeness and contamination quality assessment tests, is exported to the **Mjolnir** server. The following custom Python script (*acetic_or_lactic.py*) is then executed on the server to sort the files based on whether they belong to the Acetic or Lactic bacterial family. The bacterial genomes are divided as such because we will only work with Acetic Acid bacteria from now on.

<details>
  <summary><b>acetic_or_lactic.py</b> <i>(see code)</i></summary>

```python
import pandas as pd
import shutil
import os

def write_sort_file(path = os.path.expanduser('~')):
    '''Writes a file ("sort_file.txt") that lists all the genomes and indicates whether each genome is part of the Lactic
    or the Acetic bacteria family.
    '''
    #gather the data 
    bacteria_log = pd.read_excel(os.path.join(path, 'Bacteria log v2.xlsx'))
    sort_dict = bacteria_log.set_index('Assembly Accession').to_dict()['Lactic/acetic']
    sort_dict = {key.split('.')[0]: value for key, value in sort_dict.items()}
    #write the file
    with open(os.path.join(path, 'sort_file.txt'), 'w') as f:
        for k, v in sort_dict.items(): 
            f.write('%s:%s\n' % (k, v))

def reorganiser(path = os.path.expanduser('~'), directory = 'Genomes'):
    '''
    Moves all the ".fna" sequence files from nested directories to the "Genomes" directory.
    '''
    species = [f for f in sorted(os.listdir(os.path.join(path, directory))) if not f.startswith('.')]
    for sp in species:
        genomes = [gen for gen in sorted(os.listdir(os.path.join(path, directory, sp))) if gen.endswith('.fna')]
        for gen in genomes:
            shutil.move(os.path.join(path, directory, sp, gen), os.path.join(path, directory))

def assigner(path = os.path.expanduser('~'), directory = 'Genomes'):
    '''
    Reads the "sort_file.txt" file that assigns every genome to either the Lactic or the Acetic bacteria family and correspondingly
    moves the genome sequence files in the "Acetic" folder or the "Lactic" folder. That way, lactic and acetic bacteria can be treated separately
    for later steps of the analysis.
    '''
    genomes = [file for file in sorted(os.listdir(os.path.join(path, directory))) if not file.startswith('.')]
    if not os.path.exists(os.path.join(path, 'Acetic')): os.makedirs(os.path.join(path, 'Acetic'))
    if not os.path.exists(os.path.join(path, 'Lactic')): os.makedirs(os.path.join(path, 'Lactic'))
    with open(os.path.join(path, 'sort_file.txt')) as f:
        lines = f.readlines()
        for line in lines:
            id, family = line.strip().split(':')
            for gen in genomes:
                if id == gen.split('.')[0]:
                    shutil.copy(os.path.join(path, directory, gen), os.path.join(path, family))

if __name__ == "__main__":
    write_sort_file()
    reorganiser()
    assigned()
```

</details>

This custom script writes a file (*sort_file.txt*) that lists all the genomes and indicates whether a genome is part of the Lactic or the Acetic bacterial family. Then, nested species folders in the *Genomes* directory are removed and all the bacterial genomes are moved to either of the newly created *Acetic* and *Lactic* directory. Genomes are assigned to a directory based on information present in the sort file. That way, lactic and acetic bacteria can be processed separately for later steps of the analysis.

## 6) Further filtering (one sequence per species)
After discussing the next stages of the project, we agreed to start by only using one genome per species in the phylogeny. Since we are comparing a multitude of species whose genomes vary quite significantly, we concluded that using multiple genomes per species would be unnecessary and would only add to the computational complexity of the task. Therefore, I imported the *Acetic* directory from the server and ran another custom Python script locally (first_tree.py) which finds the most complete genome for each species and copies it to a new directory named *Acetic_unique*.

<details>
  <summary><b>first_tree.py</b> <i>(see code)</i></summary>

```python
import pandas as pd
import shutil
import os

def get_most_complete_genomes(Acetic_or_Lactic = 'Acetic', data = 'Bacteria log v2.xlsx', path = '/Users/thibaultbret/', Strict = True, Complete = True):
    global most_complete_genomes
    global rows_to_drop
    '''
    Get the most complete genome for each species selected to construct the tree.
    '''
    #make sure the path is the correct one
    if 'Python' in path: path = path.replace('Python','')
    #retrieve data to find the completeness of every genome
    bacteria_log = pd.read_excel(path + data)
    Acetic = bacteria_log[bacteria_log['Lactic/acetic'] == Acetic_or_Lactic]
    most_complete_genomes = Acetic.sort_values('Completeness', ascending = False).drop_duplicates('Unnamed: 0', keep = 'first')[['Unnamed: 0','Assembly Accession']].set_index('Unnamed: 0')
    #if Strict is set to True, remove certain species for the final tree
    if Strict == True:
        to_drop = ['Gluconobacter albidus','Gluconobacter cerevisiae','Gluconobacter japonicus','Gluconobacter oxydans','Neokomagataea thailandica','Neokomagataea tanensis','Neokomagataea anthophila']
        rows_to_drop = most_complete_genomes.loc[to_drop].index
        most_complete_genomes.drop(rows_to_drop, inplace = True)
    #return the assembly accession of the most complete genomes without moving the files
    if Complete != True: return most_complete_genomes
    #get most complete genome files from the Acetic/Lactic directory
    bacteria = [file for file in sorted(os.listdir(path + Acetic_or_Lactic)) if not file.startswith('.')]
    #create new directory containing only the most complete genome per species
    if not os.path.exists(path + Acetic_or_Lactic + '_unique'): os.makedirs(path + Acetic_or_Lactic + '_unique')
    #move the most complete genome files to the new directory
    for bact in bacteria:
        for assemb in most_complete_genomes['Assembly Accession']:
            if bact.startswith(assemb):
                shutil.copyfile(path + Acetic_or_Lactic + '/' + bact, path + '/' + Acetic_or_Lactic + '_unique/' + bact)

if __name__ == "__main__":
    get_most_complete_genomes(data = 'bacteria_log2.xlsx', Complete = False, Strict = False)
```

</details>

This script retrieves the percentage of completeness of all genomes from the "bacteria log v2.xlsx" Excel sheet (the updated version of "bacteria_log.xlsx"). For each species of interest, it finds the sequence file corresponding to the most complete genome of that species. Then it creates a new directory (*Acetic_unique*) containing only the most complete genome per species.

This directory is then exported to the server so that we can proceed with the construction of the first phylogeny.

## 7) Genome annotation using [Prokka](https://github.com/tseemann/prokka)
After exporting the *Acetic_unique* directory to the **Mjolnir** server and while looking into programs that could be used to run a pangenomic analysis, we performed genome annotation using a command line tool called **Prokka**. Genome annotation was initially meant to generate files that could be used to run a pangenomic analysis with **Roary**. However, when we decided against using Roary (as it is meant to be used on sequences coming from a single species), the Prokka annotation still retained its significance by generating files containing data that would be used later in the pipeline to produce complementary plots.

<details>
  <summary><b>prokka.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 1G --output=prokka_output.txt --time 12:00:00
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
```

</details>

**What is Prokka?**
- Command line software tool that can be installed on any Unix system. Prokka coordinates a suite of existing software tools to achieve a rich and reliable annotation of genomic bacterial sequences.
- Prokka expects preassembled genomic DNA sequences in FASTA format. Finished sequences without gaps are the ideal input, but it is expected that the typical input will be a set of scaffold sequences produced by de novo assembly software.
- The tools used are: **Prodigal** (Hyatt 2010) to identify coding sequence (CDS), **RNAmmer** (Lagesen et al., 2007) to identify ribosomal RNA genes (rRNA), **Aragorn** (Laslett and Canback, 2004) to identify transfer RNA genes, **SignalP** (Petersen et al., 2011) to identify signal leader peptides and **Infernal** (Kolbe and Eddy, 2011) for non-coding RNA.
- The traditional way to predict what a gene codes for is to compare it with a large database of known sequences, usually at a protein sequence level, and transfer the annotation of the best significant match.
- Prokka uses this method, but in a hierarchical manner, starting with a smaller trustworthy database, moving to medium-sized but domain-specific databases, and finally to curated models of protein families.
- Prokka produces 10 files in the specified output directory, all with a common prefix. The GFF v3 file (.gff) containing sequences and annotations is the one that will be used later in the pipeline.

# January-February 2023 - First attempts at Pangenomics
Making the first tree was a long process as a lot of options were explored including **MUSCLE**, **MAUVE**, **Roary** and **PEPPAN**. We finally decided to use **Anvi'o** after assessing that it was the most suitable tool in this situation.

# March-May 2023 - Anvi'o Pangenomics
> These commands are run locally after running the **anvio.sh** script on the server (which generates contigs databases for every genome and the 'external-genomes.txt' file) and after downloading the folder containing those files from the server.

**1.** Activate the Anvi'o Conda environment.
```bash
conda activate anvio-7.1
```

**2.** Generate a genomes storage using **anvi-gen-genomes-storage**. The genomes storage is generated from contigs databases corresponding to every genome sequence which are accessed using the 'external-genomes.txt' file. The genomes storage is a needed input file for later steps.
```bash
anvi-gen-genomes-storage -e external-genomes.txt -o STORAGE-GENOMES.db 
```

**3.** With the genomes storage ready, we can use the program **anvi-pan-genome** to run the actual pangenomic analysis.
```bash
anvi-pan-genome -g STORAGE-GENOMES.db -n PANGENOME
```

**4.** Move the genomes storage to the newly created Pangenome folder and change the working directory to this same folder. The rest of the analysis will be performed within that directory.
```bash
mv STORAGE-GENOMES.db PANGENOME
cd PANGENOME
```

**5.** Compute both the geometric homogeneity and functional homogeneity for the gene clusters in a pangenome database and add this information to the database. Since the phylogenemic inference cannot be performed on the entire pangenome, we will instead only use gene clusters with significant variation (combined homogeneity < 0.75).
```bash
anvi-compute-gene-cluster-homogeneity -p PANGENOME-PAN.db -g STORAGE-GENOMES.db -o homogeneity_output.txt --store-in-db
```

**6.** After the pangenomic analysis is done and the homogeneity values have been computed, we can use the program **anvi-display-pan** to display the results.
```bash
anvi-display-pan -p PANGENOME-PAN.db -g STORAGE-GENOMES.db
```

*Optional*: Get a FASTA file with aligned and concatenated amino acid sequences corresponding to all gene clusters found in the pangenome. This is optional since we cannot build the phylogeny using all gene clusters (as the alignment size would be too enormous).
```bash
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o concatenated-proteins.fa --max-num-genes-from-each-genome 1
```

*Optional*: Get the same FASTA file with non-concatenated amino acids.
```bash
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db -o genes.fa
```

**7.** Get a FASTA file with aligned and concatenated amino acid sequences corresponding to the selected gene clusters. This will be used to perform the phylogenomic analysis. We set the *--max-combined-homogeneity-index* to 0.75 to limit our selection to highly variable gene clusters (as they will have a bigger impact on the phylogeny than gene clusters with low variability). We also set the *genomes-gene-cluster-occurs* parameter to 32 as we have 32 genomes in the analysis and we want the gene clusters to be present in every genome.
```bash
anvi-get-sequences-for-gene-clusters -g STORAGE-GENOMES.db -p PANGENOME-PAN.db --concatenate-gene-clusters -o filtered-concatenated-proteins.fa --max-num-genes-from-each-genome 1 --min-num-genes-from-each-genome 1 --min-num-genomes-gene-cluster-occurs 32 --max-combined-homogeneity-index 0.75
```

**8.** Perform phylogenetic inference based on the previously generated FASTA file containing aligned and concatenated gene clusters of interest.
```bash
anvi-gen-phylogenomic-tree -f filtered-concatenated-proteins.fa -o tree.newick
```

**9.** Display the phylogenetic tree in the Anvi'o interactive interface.
```bash
anvi-interactive -p phylogenomic-profile.db -t tree.newick --title "Pangenome tree" --manual
```

# April-May 2023 - Complementary plots
## CDS content
Using the ".tsv" file generated by the **Prokka** annotation earlier in the pipeline, we can assess the proportion of coding sequences (CDS) i.e. regions of a gene's DNA that code for proteins. This is done using the following custom Slurm script (*CDS_summary.sh*).

<details>
  <summary><b>CDS_summary.sh</b> <i>(see code)</i></summary>

```bash
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
```

</details>

This script returns a text file named *CDS_summary.txt* which contains 3 columns. The first one includes the accession ID of a genome, the second corresponds to the total size of the genome and the third column shows the length of CDS in this same genome.

## G-C content
Next, we can assess the G-C content for each genome using the following command, executed on the **Mjolnir** server in the *Acetic_unique* directory:
```bash
for file in *.fa ; do gc=$(awk '/^>/ {next;} {gc+=gsub(/[GCgc]/,""); at+=gsub(/[ATat]/,"")} END {print (gc/(gc+at))*100}' $file) ; echo "${file%.fa}:$gc" >> GC_contents.txt ; done
```

# March-May 2023 - Single-copy core genes tree
After discussing our current methods and their limitations, we decided to try an alternative to the pangenome approach: a tree based on a predetermined set of 71 bacteria-specific single-copy core genes (SCCGs) i.e. genes that are present in exactly one copy in the genome. These SCCGs are included in *Bacteria_71*, an HMM profile provided by **Anvi'o**. SCCGs provide a less biased basis for tree construction as they minimise the risk that the genes being considered are not under similar evolutionary pressures (one of the built-in assumptions in phylogenetics). If we were to use genes regardless of the presence of multiple gene copies within the same genome, as we did with the pangenome tree construction, we would be more likely to violate that assumption because different copies of the same gene may experience different evolutionary pressures. 

We can reuse some of our previous code here, in this custom Slurm script (*singlecopycore.sh*).

<details>
  <summary><b>singlecopycore.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_sccg.txt --time 16:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1
# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p
# Optionally migrate the Contigs databases if they are not up to date
# anvi-migrate --migrate-safely *.db
# 2) Run HMMs for 71 Bacteria copy core genes
for file in $(ls *.db); do anvi-run-hmms -c $file -T 4 -I Bacteria_71 ; done
# 3) Estimating metabolism
anvi-get-sequences-for-hmm-hits -e external-genomes.txt --hmm-sources Bacteria_71 --get-aa-sequences --concatenate-genes -o ALIGN --return-best-hit
# 4) Constructing the tree with IQTree
module load iqtree
iqtree -s ALIGN
```

</details>

Repeat the tree construction with **IQTree** using bootstrapping (1,000 bootstraps).

<details>
  <summary><b>singlecopycore-bs.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_sccg.txt --time 16:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1
# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p
# 2) Constructing the tree with IQTree
module load iqtree
iqtree -s ALIGN-bp -bb 1000 -msub nuclear
```

</details>

Below is the single-copy core genes tree with support values:
![sccgtree-bs](https://github.com/THibaultBret/Fermentation-bacteria/assets/90853477/1da84c37-e22c-40a9-b702-c5010dc300a9)


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
*Note: The rest of the code was run from the server in a custom Slurm script called metabolics.sh (shown below).*

<details>
  <summary><b>metabolics.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_meta.txt --time 12:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio/7.1

# 1) Setup
p=/projects/mjolnir1/people/vhp327/Acetic_all
cd $p
# 2) Annotating the genome with KOfam hits
for file in $(ls *.db); do anvi-run-hmms -c $file  -H HMM_AAB ; done
# 3) Estimating metabolism
anvi-get-sequences-for-hmm-hits --external-genomes external-genomes.txt \
                                -o pfam-proteins.fa \
                                --hmm-source HMM_AAB \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate
# 4) Compute the phylogenomic tree
anvi-gen-phylogenomic-tree -f pfam-proteins.fa \
                           -o pfam-phylogenomic-tree.txt
# 5) build iq-tree
module load iqtree
iqtree -s pfam-proteins.fa
```

</details>

4. Annotate the genomes with HMMs hits using the **anvi-run-hmms** command.
5. Now we can go through every contigs database (listed in the file *external-genomes.txt*) and retrieve sequences that got a hit for any given HMM using the **anvi-get-sequences-for-hmm-hits** command. The sequences are concatenated so that a phylogenetic tree can be drawn from them (as aligned sequences are a prerequisite). Furthermore, we select on the best hits using the *--return-best-hit* flag. From the Anvi'o manual: "This flag is most appropriate if one wishes to perform phylogenomic analyses, which ensures that for any given protein family, there will be only one gene reported from a given genome". The concatenated best hits are returned in the *pfam-proteins.fa* file.
6. Now a phylogenetic tree can be drawn from the *pfam-proteins.fa* file using the **anvi-gen-phylogenomic-tree** command. However, this will generate a rough tree which is not entirely reliable. Alternatively, a tree can drawn using IQ-TREE simply by using the following command:
```bash
iqtree -s pfam-proteins.fa
```
*Note: However, after trying to run this code, I noticed that the model selection process runs indefinitely (possibly because of the number of sequences the tree is based on). Therefore, I decided to stick to the Anvi'o-generated tree.*

Below is the general Pfam tree:
<img width="659" alt="pfam-general-tree" src="https://github.com/THibaultBret/Fermentation-bacteria/assets/90853477/251cda9d-4bac-4905-b49f-f7e5f09f9b7f">

And a comparison of the general Pfam tree (on the left) and the SCCG tree (on the right):
<img width="1361" alt="Comparison of trees" src="https://github.com/THibaultBret/Fermentation-bacteria/assets/90853477/a9828d33-b7e4-4254-be3f-6578613130f9">


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
The following custom Slurm script (*pfamtrees.sh*) is executed in the *PFAM_seqs* directory, located on the **Mjolnir** server, which contains all the Pfam HMMs previously generated. This script generates phylogenies based on each individual Pfams with bootstrapping.

<details>
  <summary><b>pfamtrees.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_iqtree.txt --time 24:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load iqtree

# Construct iq-trees
for seq in $(ls *.fa); do iqtree -s $seq -bb 1000 ; done
```

</details>

The resulting trees are shown in the *pfamtrees-bs.pdf* file.


## Automatically download new genomes using [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)

Using the NCBI datasets command line tool, new genomes can be downloaded automatically in a matter of minutes.

```bash
datasets download genome accession GCF_000010825.1 --include gff3,rna,cds,protein,genome,seq-report
```

Simply replace the accession with the ones corresponding to the genome you want to download.


If you wanna easily download genome sequences only:
```bash
datasets download genome accession <accessions> --include genome
```

Then unzip the resulting zip file, rename the "data" folder to "FinalGenomes" and run

```bash
# Specify the source directory
source_dir="FinalGenomes"

# Move files ending with .fna to FinalGenomes directory
for dir in "$source_dir"/*/; do
    mv "$dir"*.fna "$source_dir"
done

# Remove empty directories
find FinalGenomes -type d -empty -delete
```

If you need to rename the fasta files, run:
```bash
for file in *.fna; do mv "$file" "${file%.fna}.fasta"; done
```

# February 2024 - Metabolic analysis

## Method 1 - Restricted metabolic analysis

The aim here is to run the metabolic analysis with only a restricted set of enzymes of interest in order to optimise the speed of the task (which won't be slowed down by factoring irrelevant enzymes). The disadvantage here is that we don't necessarily know which enzymes might be relevant so we might be missing out on important data, but it's also a way to experiment with the  metabolism suite of programs in Anvi’o. All these commands can be run locally using the **anvio-dev** conda working environment. Our working directory is named **FinalTree** and contains our tree files ("sccg-tree-noSRR.treefile" and "sccg-tree.treefile"), our external-genomes files ("external-genomes.txt" and "external-genomes-noSRR.txt") as well as a folder named **ContigsDB** containing all the contigs databases we have been working with so far.

### 1) Setting up KEGG data
The first step consists of defining the working directory and setting up a directory containing all the KEGG data.
```bash
p=/Users/thibaultbret/All-db
cd $p
anvi-setup-kegg-kofams --kegg-data-dir $p/KEGG
```

### 2) Generate Hidden Markov Models (HMMs) corresponding to each Pfam sequence
Then we can create Hidden Markov Model directories corresponding to each Pfam sequence we are interested in. This Anvi'o script will generate this for us, we only need to provide it with a list of the Pfam accessions and a name for the output directory.
```bash
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list "PF00005" "PF00196" "PF00923" "PF02887" "PF08240" "PF00118" "PF00330" "PF00958" "PF03070" "PF13243" "PF00171" "PF00465" "PF01161" "PF08042" "PF13360" -O HMM_AAB
```

### Optional - update databases if you have to
It is sometimes necessary to migrate the CONTIGS databases if Anvi'o has been updated to a newer version.
```bash
anvi-migrate --migrate-safely ContigsDB/*_CONTIGS.db
```

### 3) Annotate the genomes with Pfam hits
We can now annotate the genomes with Pfam hits based on the HMM directories we created earlier.
```bash
for file in $(ls ContigsDB/*_CONTIGS.db); do anvi-run-hmms -c $file -H HMM_AAB --add-to-functions-table --just-do-it ; done
```

### 4) Generate a user-defined module file
We need to generate a module file outlining the enzymes we are conducting this analysis with, their HMM source and orthology. It should look as follows:
enzymes-list-for-module.txt
```bash
enzyme	source	orthology
PF00005	HMM_AAB	"ABC transporter"
PF00196	HMM_AAB	"Bacterial regulatory proteins, luxR family"
PF00923	HMM_AAB	"Transaldolase/Fructose-6-phosphate aldolase"
PF02887	HMM_AAB	"Pyruvate kinase, alpha/beta domain"
PF08240	HMM_AAB	"Alcohol dehydrogenase GroES-like domain"
PF00118	HMM_AAB	"TCP-1/cpn60/GroEL chaperonin family"
PF00330	HMM_AAB	"Aconitase family (aconitate hydratase)"
PF00958	HMM_AAB	"GMP synthase C terminal domain"
PF03070	HMM_AAB	"TENA/THI-4/PQQC family"
PF13243	HMM_AAB	"Squalene-hopene cyclase C-terminal domain"
PF00171	HMM_AAB	"Aldehyde dehydrogenase family"
PF00465	HMM_AAB	"Iron-containing alcohol dehydrogenase"
PF01161	HMM_AAB	"Phosphatidylethanolamine-binding protein"
PF08042	HMM_AAB	"PqqA family"
PF13360	HMM_AAB	"PQQ-like domain"
```

```bash
anvi-script-gen-user-module-file -I "AAB.txt" \
                  -n "Collection of 15 enzyme families present in AAB" \
                  -c "User modules; Test set; AAB enzymes" \
                  -e enzymes-list-for-module.txt
```


### 5) Make modules database
Necessary steps before estimating metabolism.
```bash
mkdir AAB_METABOLISM && cd AAB_METABOLISM
mkdir modules
mv $p/AAB.txt modules
cd ..
anvi-setup-user-modules --user-modules AAB_METABOLISM --kegg-data-dir ~/FinalTree/KEGG
```


### 6) Estimate metabolism
```bash
anvi-estimate-metabolism  -e external-genomes.txt \
                         --user-modules AAB_metabolism/ \
                         --only-user-modules \
                         -O AAB-metabolism-results \
                         --matrix-format
```

### 7) Visualise the result
The previous step will generate a file named "AAB-metabolism-results-enzyme_hits-MATRIX.txt". We could directly feed it into the next command and obtain a rough visual representation but to make sure that the heatmap includes as much information as possible, it is preferable to copy that data and transpose it into an Excel sheet. Then using previously calculated parameters such as GC content, CDS proportion, genome size and isolation source, this Excel file can serve as the basis for the heatmap. To feed it back into Anvi'o we first need to convert it back to a .txt file that we will name "AAB-metabolics-data.txt".

Initiate the graph:
```bash
anvi-interactive -d AAB-metabolics-data.txt -p AAB-metabolism-heatmap.db --manual
```

### 8) Add the phylogenetic tree
```bash
anvi-import-items-order -i sccg-tree-noSRR.treefile \
                        -p AAB-metabolism-heatmap.db \
                        --name taxonomy

anvi-interactive -d AAB-metabolics-data.txt -p AAB-metabolism-heatmap.db --manual
```

The result:
![Metabolics-heatmap-noSRR](https://github.com/THibaultBret/Fermentation-bacteria/assets/90853477/852ed5d9-400f-4879-848f-76371dcabe89)


## Method 2 - Full metabolic analysis

The aim here is to run the metabolic analysis without first specifying enzymes or pathways of interest so that every factor can be taken into consideration. This script has to be run from the server as it is a very computationally heavy task. After setting the working directory, optionally migrating outdated databases and setting up a directory with the KEGG data, the first step is to annotate the genomes with KOfam hits, which takes between 5 and 15 minutes for each contigs database. Once this is done for all 154 genomes, we can run the **anvi-estimate-metabolism** Anvi'o command which can either return a "matrix-format" output which is useful if we want to determine which pathways/enzymes are associated with a given genome, but it's a lot of information, or we can leave it to the default output format, which may be less interesting to look at but is necessary to later run the enrichment analysis.

<details>
  <summary><b>anvio-metabolics.sh</b> <i>(see code)</i></summary>

```bash
#!/bin/sh
#SBATCH -c 8 --mem 4G --output=output_anvio.txt --time 2-00:00:00 --mail-user=thibault.bret@sund.ku.dk --mail-type=ALL
module load anvio
module load anvio-dev

# Setup
p=/projects/mjolnir1/people/vhp327/FinalTree/
cd $p

# Optionally migrate outdated databases
#anvi-migrate --migrate-safely *_CONTIGS.db
#anvi-migrate --migrate-safely KEGG/MODULES.db

# 1) Set up KEGG data if not done already
#anvi-setup-kegg-kofams --kegg-data-dir $p/KEGG

# 2) Annotating the genome with KOfam hits
for file in $(ls ContigsDB/*_CONTIGS.db); do
    anvi-run-kegg-kofams -c ${file} -T 4 --kegg-data-dir $p/KEGG
done

# 3) Estimating metabolism (pick preferred output format)
#anvi-estimate-metabolism -e external-genomes.txt --kegg-data-dir $p/KEGG --matrix-format --include-metadata
anvi-estimate-metabolism -e external-genomes.txt --kegg-data-dir $p/KEGG
```

</details>

Running the **anvi-estimate-metabolism** with default setting will return a file named "kegg-metabolism_modules.txt". This file can be used as the input for a very import command called **anvi-compute-metabolic-enrichment** which will compute an enrichment score and a list of associated groups for each module that is present in at least one genome (modules are considered ‘present’ in a genome if they have a high enough completeness score in that genome). It also requires a text file that will describe the groups between which we want to separate the genomes. It should look like the following:

```bash
item	group
Acetobacter_ascendens	Ferment
Acetobacter_cibinongensis	Plant
Acetobacter_conturbans	Ferment
Acetobacter_fallax	Ferment
Acetobacter_farinalis	Ferment
Acetobacter_garciniae	Ferment
Acetobacter_ghanensis	Ferment
Acetobacter_lambici	Ferment
Acetobacter_lovaniensis	Plant
Acetobacter_malorum	Ferment
...
```

We can finally run the following command, either on the server or in the **anvio-dev** conda environment:

Run the enrichment analysis based on isolation source groups:

```bash
anvi-compute-metabolic-enrichment -M kegg-metabolism_modules.txt \
                                  -G groups-isolation.txt \
                                  -o functional-enrichment-isolation-groups.txt
```

Run the enrichment analysis based on genus groups:

```bash
anvi-compute-metabolic-enrichment -M kegg-metabolism_modules.txt \
                                  -G groups-genus.txt \
                                  -o functional-enrichment-genus-groups.txt
```

### Visualise the result
Now based on the pathways that are significantly enriched in a group of set of groups compared to others, we can generate a new heatmap that will be much more informative.
We need the "AAB-metabolism-results-enzyme_hits-MATRIX.txt" returned by the "anvi-estimate-metabolism" command with the "--matrix-format" and "include-metadata" settings. Again for better visualisation, it is preferable to copy that data and transpose it into an Excel sheet. Then using previously calculated parameters such as GC content, CDS proportion, genome size and isolation source, this Excel file can serve as the basis for the heatmap. To feed it back into Anvi'o we first need to convert it back to a .txt file that we will name "AAB-full-KEGG-data.txt".

Initiate the graph:
```bash
anvi-interactive -d kegg-metabolism-module_pathwise_completeness-MATRIX.txt -p AAB-full-KEGG-heatmap.db --manual
```

### Sort metabolic pathways & add custom bins
```bash
# Cluster the pathways so that metabolisms with similar distributions across MAGs will be closer together
anvi-matrix-to-newick kegg-metabolism-module_pathwise_completeness-MATRIX.txt -o module_organization.newick

# Cluster the genomes according to their metabolic capacity, so that MAGs with similar capacities will be closer together
anvi-script-transpose-matrix kegg-metabolism-module_pathwise_completeness-MATRIX.txt -o kegg-metabolism-module_pathwise_completeness-MATRIX-TRANSPOSED.txt
anvi-matrix-to-newick kegg-metabolism-module_pathwise_completeness-MATRIX-TRANSPOSED.txt -o mag_organization.newick

# Now we have two dendrograms, one for each side of the heatmap. Let’s add them to the display by importing them into the profile database
anvi-import-items-order -i module_organization.newick \
                        -p AAB-full-KEGG-heatmap.db \
                        --name module_organization

# Layer orders have to be imported into the database using anvi-import-misc-data, which requires copying the tree into a tab-delimited file that also describes the name and type of the ordering.
TREE=$(cat mag_organization.newick)  # copy the dendrogram into a variable called 'TREE'
echo -e "item_name\tdata_type\tdata_value\nmag_organization\tnewick\t$TREE" > layer_order.txt # put it into a file
# import into the database
anvi-import-misc-data -p AAB-full-KEGG-heatmap.db \
                      -t layer_orders \
                      layer_order.txt

# Import custom bins
anvi-import-collection gengroups.txt -p AAB-full-KEGG-heatmap.db -C gengroups --bins-info bins-info.txt

# Show the figure
anvi-interactive -d kegg-metabolism-module_pathwise_completeness-MATRIX.txt -p AAB-full-KEGG-heatmap.db --manual
```



### Only looking at pathways enriched in certain genera
```bash
anvi-interactive -d AAB-genus-KEGG-data.txt -p AAB-genus-KEGG-heatmap.db --manual
```

```bash
anvi-import-items-order -i sccg-tree-noSRR.treefile \
                        -p AAB-genus-KEGG-heatmap.db \
                        --name taxonomy

anvi-import-collection gengroups.txt -p AAB-genus-KEGG-heatmap.db -C gengroups --bins-info bins-info.txt

anvi-import-state -p AAB-genus-KEGG-heatmap.db -s state -n default

anvi-interactive -d AAB-genus-KEGG-data.txt -p AAB-genus-KEGG-heatmap.db --manual
```

### Only looking at pathways enriched in certain isolation groups
```bash
anvi-interactive -d AAB-isolation-KEGG-data.txt -p AAB-isolation-KEGG-heatmap.db --manual
```

```bash
anvi-import-items-order -i sccg-tree-noSRR.treefile \
                        -p AAB-isolation-KEGG-heatmap.db \
                        --name taxonomy

anvi-import-state -p AAB-isolation-KEGG-heatmap.db -s state -n default

anvi-import-collection gengroups.txt -p AAB-isolation-KEGG-heatmap.db -C gengroups --bins-info bins-info.txt

anvi-interactive -d AAB-isolation-KEGG-data.txt -p AAB-isolation-KEGG-heatmap.db --manual
```


# February-March 2024 - Time scaled phylogeny & concordance factors

Our single-copy core gene trees are already informative but one thing we can add is a time scale. A time scale can provide a rough chronology of the splits in the tree and potentially corroborate with existing hypotheses. To generate a time scale we first have to calibrate it based on pre-existing knowledge. We can give it two nodes of the tree .


