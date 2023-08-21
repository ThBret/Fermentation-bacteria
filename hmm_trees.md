## July-August 2023 - New trees based on selected genes and pathways specific to acetic acid bacteria:

# 1) General tree
1. List pathways/proteins/enzymes of interest and find the corresponding Pfam accessions.
2. *Note: The following code can only be run with **anvio-dev**, a corresponding module exists on the server but I could not make it work so I installed a version of **anvio-dev** locally following this tutorial:* https://anvio.org/install/#5-follow-the-active-development-youre-a-wizard-arry
3. Once everything is properly set up, activate the anvio dev environment with:
```bash
conda activate anvio-dev
 ```
4. Now, given a list of Pfam accessions of interest, run the **anvi-script-pfam-accessions-to-hmms-directory** command to generate **Hidden Markov Models (HMMs)** corresponding to each Pfam sequence. This will allow us later on to scan through any contigs database and detect the presence of such sequences. An output folder for the HMMs also has to be provided, here we name it *HMM_AAB* (for Hidden Markov Model - Acetic Acid Bacteria).
```bash
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list PF03070 PF08042 PF13360 PF00171 PF00465 PF02317 PF00005 PF02887 PF00923 PF03971 PF00168 PF01161 PF13243 PF00330 PF17327 PF00196 PF00958 PF08240 PF00118  -O HMM_AAB
 ```
6. Now, given this HMM, can go through every contigs database (listed in external-genomes.txt) and retrieve corresponding sequences (sequences that get a hit given the HMM). The sequences are concatenated so that a phylogenetic tree can be drawn from them (aligned sequences are a prerequisite).
7. Now a phylogenetic tree can be drawn from the pfam-proteins.fa file using the anvi-gen-phylogenomic-tree command. However, this will generate a rough tree which is not entirely reliable.
8. Alternatively, can draw the tree using IQ-tree simply by using the following command: iqtree -s pfam-proteins.fa

Individual trees!!
# 1) Set up HMMs
pfams=("PF03070" "PF08042" "PF13360" "PF00171" "PF00465" "PF02317" "PF00005" "PF02887" "PF00923" "PF03971" "PF00168" "PF01161" "PF13243" "PF00330" "PF17327" "PF00196" "PF00958" "PF08240" "PF00118")
for pfam in "${pfams[@]}"; do anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list $pfam -O HMMs/${pfam}_HMM ; done

# 2) Annotating the genome with KOfam hits
for hmm in $(ls HMMs); do for file in *CONTIGS.db; do anvi-run-hmms -c "$file" -H HMMs/$hmm -T 4; done ; done

# 3) Estimating metabolism
cd HMMs
for hmm in $(ls); do anvi-get-sequences-for-hmm-hits --external-genomes ../external-genomes.txt -o  ../PFAM_seqs/${hmm%%_HMM}-proteins.fa --hmm-source $hmm --return-best-hit --get-aa-sequences --concatenate ; done
