##################################################
### Ana Cuesta Genome characteristic AAB paper ###
##################################################
##Clean your environment
rm(list=ls())

library("tidyverse")
library("readxl")
library(ggplot2)
library(ggrepel)
library(readxl)
library("ggpubr")
library("FSA")
library(cowplot)
library(ggpmisc)
library(grid)
library(gridExtra)
library(car)
library("multcompView")
library(dplyr)
library(ggplot2)
library(multcompView)
library(FSA)
library(ggplot2)

###Prepare the data
#set paths
basepath = "/Users/fmg115/Desktop/Phylo"
excelpath = paste(basepath,"/AAB Apr2024 BoxplotsOnly.xlsx",sep="")
all_features <- read_excel("/Users/fmg115/Desktop/Phylo/AAB Apr2024 BoxplotsOnly.xlsx", sheet = "AAB genomes")

#read data
all_features <- read_excel(excelpath, sheet = "AAB genomes", range = cell_cols("A:AF"))
#clean the data
all_features <- subset(all_features, Isolation_source != "NA")
all_features <- subset(all_features, Isolated_from != "Too small" & Isolated_from != "Too large" & Isolated_from != "Not viable genome")
all_features$Isolation_source <- factor(all_features$Isolation_source, levels=c('Insect','Fly','Plant','Ferment','Industrial','Acidophilic'))
abbreviate <- function(col) paste(substr(gsub( " .*$", "", col),start=0,stop=1), sub("^\\S+\\s+", '',  col), sep =  ". ")
all_features$Species <- lapply(all_features$Species, abbreviate)

#check the data
all_features
head(all_features)



#########################
### Association Plots ###
#########################
#geom_text_repel labels your data points and can be useful for a rough observation of data distribution in your plots

b <- ggplot(all_features, aes(x = Size, y = GC_content))

# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "loess") 

ggscatter(all_features, x = "Size", y = "GC_content",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal())

# Change color and shape by groups (Isolation source)
b + geom_point(aes(color = Isolation_source))+
  geom_smooth(aes(color = , fill = Isolation_source), method = "lm") +
  geom_rug(aes(color =Isolation_source)) +
  scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596"))+
  scale_fill_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596"))=
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 40)

# Add regression statistics 
b + geom_point(aes(color = Isolation_source))+
  geom_rug(aes(color =Isolation_source)) +
  geom_smooth(aes(color = Isolation_source), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596"))+
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 40) +
  ggpubr::stat_cor(aes(color = Isolation_source), label.x = 3)

# Split by groups
b + geom_point(aes(color = Isolation_source), alpha = 0.9)+
  geom_smooth(aes(color = Isolation_source, fill = Isolation_source), 
              method = "lm", fullrange = TRUE) +
  facet_wrap(~Isolation_source) +
  scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596"))+
  scale_fill_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596")) +
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 40) +
  theme_classic()

# Ellipses
b + geom_point(aes(color = Isolation_source), alpha = 0.9, size = 1)+
  stat_ellipse(aes(color = Isolation_source), type = "t", linewidth = 1)+
  scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF", "#289596"))+
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 40) +
  theme_bw()


##################################################
#### Association Plots between relevant clades ###
##################################################
#Remove the acidophillic group 
all_features <- subset(all_features, Isolation_source != "Acidophillic")
## Komagataeibacter vs Gluconacetobacter
#New data with the subset of clades that we want to include
excelpath = paste(basepath,"/KomagataeibactervsGluconacetobacter.xlsx",sep="")
KomGlu <- read_excel(excelpath, range = cell_cols("A:AG"))

b <- ggplot(KomGlu, aes(x = Size, y = GC_content))
# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "loess") 

ggscatter(KomGlu, x = "Size", y = "GC_content",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal())

b + geom_point(aes(color = Clade), alpha = 0.9, size = 1)+
  stat_ellipse(aes(color = Clade), type = "t", linewidth = 1)+
  scale_color_manual(values = c("#73C221", "#9933CD"))+
  #geom_text_repel(aes(label = Species, colour = Clade), size = 3, max.overlaps = 40) +
  theme_bw()
b
b <- b +theme_bw()
b <- b + scale_y_continuous(limits = c(20, 80))
b <- b + scale_x_continuous(limits = c(1000000, 6000000))
b <- b + theme(legend.position = "none")
print(b)

## Bombella vs Commensalibacter
#New data with the subset of clades that we want to include
excelpath = paste(basepath,"/BombellavsCommensalibacter.xlsx",sep="")
BombComm <- read_excel(excelpath, range = cell_cols("A:AG"))

b <- ggplot(BombComm, aes(x = Size, y = GC_content))
# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "loess") 

ggscatter(BombComm, x = "Size", y = "GC_content",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          ggtheme = theme_minimal())

b + geom_point(aes(color = Clade), alpha = 0.9, size = 1)+
  stat_ellipse(aes(color = Clade), type = "t", linewidth = 1)+
  scale_color_manual(values = c("#FFDE00", "#D8A438"))+
  #geom_text_repel(aes(label = Species, colour = Clade), size = 3, max.overlaps = 40) +
  theme_bw()
b
b <- b +theme_bw()
b <- b + scale_y_continuous(limits = c(20, 80))
b <- b + scale_x_continuous(limits = c(1000000, 6000000))
b <- b + theme(legend.position = "none")
print(b)




########################
### Genome size Plot ###
########################
p <- ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values =  c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) + scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16) +
  #geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 70,  position=pos) +
  xlab("Isolation source") + ylab("Genome size (in Mb)") + ggtitle("Genome sizes of Acetic Acid bacteria genomes") + 
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 18), axis.text = element_text(size = 12))
p
p <- p +theme_bw()
p <- p + theme(legend.position = "none")
print(p)

#########################
### GC content % Plot ###
#########################
p <- ggplot(all_features, aes(x=Isolation_source, y=as.numeric(GC_content), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) + scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16) +
  #geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 70,  position=pos) +
  xlab("Isolation source") + ylab("G-C content (in %)") + ggtitle("G-C contents of Acetic Acid bacteria genomes") + 
  theme(legend.position = "none", plot.title = element_text(size = 22))
p
p <- p +theme_bw()
p <- p + theme(legend.position = "none")
print(p)

###########################
### Average gene length ###
###########################
p <- ggplot(all_features, aes(x=Isolation_source, y=Average_genelenght, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, size = 2.5) +
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 1, max.overlaps = 3,  ) +
  scale_fill_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF"))+
  scale_color_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF"))+
  theme(legend.position = "none")
p
p <- p +theme_bw()
p <- p + theme(legend.position = "none")
print(p)

######################################
### Proportion of coding sequences ###
######################################
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF")) + scale_color_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF")) +
  geom_jitter(aes(colour = Isolation_source, size = 0.0001, alpha = .9), shape=16, position=pos) +
  #geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 70,  position=pos) +
  xlab("Isolation source") + ylab("CDS proportion (in %)") + ggtitle("Proportion of Coding Sequences (CDS) in Acetic Acid bacteria genomes") + 
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 18), axis.text = element_text(size = 12))+
  theme_bw()
p
p <- p +theme_bw()
p <- p + theme(legend.position = "none")
print(p)

#######################
### Number of genes ###
#######################
b <- ggplot(all_features, aes(x=Isolation_source, y=Total_gene_count, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, size = 2.5) +
  #geom_text_repel(aes(label = Species, colour = Isolation_source), size = 1, max.overlaps = 3,  ) +
  scale_fill_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF"))+
  scale_color_manual(values = c("#E6A401","#827E6F","#73C221", "#9933CD", "#FF40FF"))+
  theme(legend.position = "none")
b <- b +theme_bw()
b <- b + theme(legend.position = "none")
print(b)

#################################
### Temperature density plots ###
#################################
##Remove all the ones we don't know Temperature growth for
all_features <- subset(all_features, Ta != "NA")

ggplot(all_features, aes(x=as.numeric(Ta), color=Isolation_source)) + 
  
  # color property for changing color of plot
  # geom_density() function plots the density plot
  geom_density(alpha = 7, linewidth = 1.5) +
  scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF"))+
  theme_classic()



#############################################
### Statistics on genomic characteristics ###
#############################################

################################################## Genome size

##Test for normality
ggdensity(all_features$Size)
ggqqplot(all_features$Size)
shapiro.test(all_features$Size) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Size~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(Size~Isolation_source, data = all_features,
         method="bonferroni")


################################################## GC content

##Test for normality
ggdensity(all_features$GC_content)
ggqqplot(all_features$GC_content)
shapiro.test(all_features$GC_content) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(GC_content~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(GC_content~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Proportion of coding sequences

##Test for normality
ggdensity(all_features$CDS_prop)
ggqqplot(all_features$CDS_prop)
shapiro.test(all_features$CDS_prop) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(CDS_prop~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(CDS_prop~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Coding sequence length (all coding sequence length)

##Test for normality
ggdensity(all_features$CDS_length)
ggqqplot(all_features$CDS_length)
shapiro.test(all_features$CDS_length) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(CDS_length~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(CDS_length~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Average gene length

##Test for normality
ggdensity(all_features$Average_genelenght)
ggqqplot(all_features$Average_genelenght)
shapiro.test(all_features$Average_genelenght) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Average_genelenght~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(Average_genelenght~Isolation_source, data = all_features,
         method="bonferroni")

### Normalized data (by genome size)

##Test for normality
all_features <-all_features[!all_features$NORM_genelengthav=="NA", ]
ggdensity(all_features$NORM_genelengthav)
ggqqplot(all_features$NORM_genelengthav)
shapiro.test(all_features$NORM_genelengthav) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(NORM_genelengthav~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(NORM_genelengthav~Isolation_source, data = all_features,
               method="bonferroni")


################################################## Number of genes

##Test for normality
ggdensity(all_features$Total_gene_count)
ggqqplot(all_features$Total_gene_count)
shapiro.test(all_features$Total_gene_count) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Total_gene_count~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(Total_gene_count~Isolation_source, data = all_features,
         method="bonferroni")

### Normalized data (by genome size)

##Test for normality
all_features <-all_features[!all_features$NORM_genecount=="NA", ]
all_features$NORM_genecount <- as.numeric(all_features$NORM_genecount)
ggdensity(all_features$NORM_genecount)
ggqqplot(all_features$NORM_genecount)
shapiro.test(all_features$NORM_genecount) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(NORM_genecount~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(NORM_genecount~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Non-coding genes

##Test for normality
all_features <-all_features[!all_features$Non_coding=="NA", ]
all_features$Non_coding <- as.numeric(all_features$Non_coding)
ggdensity(all_features$Non_coding)
ggqqplot(all_features$Non_coding)
shapiro.test(all_features$Non_coding) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Non_coding~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(Non_coding~Isolation_source, data = all_features,
         method="bonferroni")

### Normalized data (by genome size)

##Test for normality
ggdensity(all_features$NORM_Non_coding)
ggqqplot(all_features$NORM_Non_coding)
shapiro.test(all_features$NORM_Non_coding) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(NORM_Non_coding~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(NORM_Non_coding~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Protein coding genes

##Test for normality
all_features <-all_features[!all_features$protein_coding=="NA", ]
all_features$protein_coding <- as.numeric(all_features$protein_coding)
ggdensity(all_features$protein_coding)
ggqqplot(all_features$protein_coding)
shapiro.test(all_features$protein_coding) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(protein_coding~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(protein_coding~Isolation_source, data = all_features,
         method="bonferroni")

### Normalized data (by genome size)

##Test for normality
ggdensity(all_features$NORM_Protein_conding)
ggqqplot(all_features$NORM_Protein_conding)
shapiro.test(all_features$NORM_Protein_conding) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(NORM_Protein_conding~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(NORM_Protein_conding~Isolation_source, data = all_features,
         method="bonferroni")


################################################## Pseudogenes

##Test for normality
all_features <-all_features[!all_features$Nb_pseudogenes=="NA", ]
all_features$Nb_pseudogenes <- as.numeric(all_features$Nb_pseudogenes)
ggdensity(all_features$Nb_pseudogenes)
ggqqplot(all_features$Nb_pseudogenes)
shapiro.test(all_features$Nb_pseudogenes) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Nb_pseudogenes~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(Nb_pseudogenes~Isolation_source, data = all_features,
         method="bonferroni")

### Normalized data (by genome size)

##Test for normality
ggdensity(all_features$NORM_Pseudogenes)
ggqqplot(all_features$NORM_Pseudogenes)
shapiro.test(all_features$NORM_Pseudogenes) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(NORM_Pseudogenes~Isolation_source, data = all_features) ##significant differences based on isolation source
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(NORM_Pseudogenes~Isolation_source, data = all_features,
         method="bonferroni")

###########################################################################################
### Statistics on genomic characteristics, pairwise comparisons between relevant groups ###
###########################################################################################

############################################################################ Komagtaeibacter vs Gluconacetobacter ##################################################

excelpath = paste(basepath,"/KomagateibactervsGluconacetobacter.xlsx",sep="")

######################### Are Komagataeibacter and Gluconacetobacter different in Size? 
##Test for normality
ggdensity(KomGlu$Size)
ggqqplot(KomGlu$Size)
shapiro.test(KomGlu$Size) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox_test (Size~Genus, data = KomGlu, 
                     exact = FALSE) ##significant differences based on clades
wilcoxonZ(x= KomGlu$Size[KomGlu$Genus=='Gluconacetobacter'], KomGlu$Size[KomGlu$Genus=='Komagaiteibacter'],
          paired = TRUE, exact = FALSE, correct = FALSE)


######################### Are Komagataeibacter and Gluconacetobacter different in GC content? 
##Test for normality
ggdensity(KomGlu$GC_content)
ggqqplot(KomGlu$GC_content)
shapiro.test(KomGlu$GC_content) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox.test (GC_content~Clade, data = KomGlu) ##significant differences based on clades
wilcoxonZ(x= KomGlu$GC_content[KomGlu$Clade=='Gluconacetobacter'], KomGlu$GC_content[KomGlu$Clade=='Komagaiteibacter'])


######################### Are Komagataeibacter and Gluconacetobacter different in Ta ranges? 
KomGlu$Ta <- as.numeric(KomGlu$Ta)

##Test for normality
ggdensity(KomGlu$Ta)
ggqqplot(KomGlu$Ta)
shapiro.test(KomGlu$Ta) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox.test (Ta~Clade, data = KomGlu) ##significant differences based on clades
wilcoxonZ(x= KomGlu$Ta[KomGlu$Clade=='Gluconacetobacter'], KomGlu$Ta[KomGlu$Clade=='Komagaiteibacter'])



################################################## Bombella vs Commensalibacter ##################################################

excelpath = paste(basepath,"/BombellavsCommensalibacter.xlsx",sep="")

######################### Are Bombella and Commensalibacter different in Size? 
##Test for normality
ggdensity(BombComm$Size)
ggqqplot(BombComm$Size)
shapiro.test(BombComm$Size) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox.test (Size~Clade, data = BombComm) ##significant differences based on clades
wilcoxonZ(x= BombComm$Size[BombComm$Clade=='Commensalibacter'], BombComm$Size[BombComm$Clade=='Bombella'])


######################### Are Bombella and Commensalibacter different in GC content? 
##Test for normality
ggdensity(BombComm$GC_content)
ggqqplot(BombComm$GC_content)
shapiro.test(BombComm$GC_content) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox.test (GC_content~Clade, data = BombComm) ##significant differences based on clades
wilcoxonZ(x= BombComm$GC_content[BombComm$Clade=='Commensalibacter'], BombComm$GC_content[BombComm$Clade=='Bombella'])


######################### Are Bombella and Commensalibacter different in Ta ranges? 
BombComm$Ta <- as.numeric(BombComm$Ta)
##Test for normality
ggdensity(BombComm$Ta)
ggqqplot(BombComm$Ta)
shapiro.test(BombComm$Ta) ## p<0.05 = Differ from normality, so non-parametric test

####Man U Whitney for non-parametric data (only two categories)
wilcox.test (Ta~Clade, data = BombComm) ##significant differences based on clades
wilcoxonZ(x= BombComm$Ta[BombComm$Clade=='Commensalibacter'], BombComm$Ta[BombComm$Clade=='Bombella'])



#############################################################################
### Statistics on genomic characteristics, comparisons between all caldes ###
#############################################################################

######################### How much do clades differ in Size? 
##Remove those species without a clade assigned
AllClades <- subset(all_features, Clades != "NA")
p <- ggplot(AllClades, aes(x=Clades, y=as.numeric(Size), fill=Clades)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Clades, alpha = .9), shape=16) +
  #geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 70,  position=pos) +
  ggtitle("Genome size") + 
  scale_fill_manual(values = c("#B1E5FF", "#37D482", "#FF0066", "#1770D4", "#9F613F", "#FF9900", "#D4FB79",
                               "#D935D9", "#21720F", "#4BFF23", "#D1C4FF", "#8EFFA8"))+
  scale_color_manual(values = c("#B1E5FF", "#37D482", "#FF0066", "#1770D4", "#9F613F", "#FF9900", "#D4FB79",
                                "#D935D9", "#21720F", "#4BFF23", "#D1C4FF", "#8EFFA8"))+
  theme(legend.position = "none", plot.title = element_text(size = 22))

##Test for normality
ggdensity(all_features$Size)
ggqqplot(all_features$Size)
shapiro.test(all_features$Size) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(Size~Clades, data = all_features) ##significant differences based on Species
####Dunn's Post-hoc test for mutiple comparisons
All_Comparisons <-dunnTest(Size~Clades, data = all_features,
                           method="bonferroni")

All_Comparisons


######################### How much do clades differ in GC content? 
q <- ggplot(all_features, aes(x=Clades, y=as.numeric(GC_content), fill=Clades)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Clades, alpha = .9), shape=16) +
  #geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 70,  position=pos) +
  ggtitle("Genome size") + 
  scale_fill_manual(values = c("#B1E5FF", "#37D482", "#FF0066", "#1770D4", "#9F613F", "#FF9900", "#D4FB79",
                               "#D935D9", "#21720F", "#4BFF23", "#D1C4FF", "#8EFFA8"))+
  scale_color_manual(values = c("#B1E5FF", "#37D482", "#FF0066", "#1770D4", "#9F613F", "#FF9900", "#D4FB79",
                                "#D935D9", "#21720F", "#4BFF23", "#D1C4FF", "#8EFFA8"))+
  theme(legend.position = "none", plot.title = element_text(size = 22))

##Test for normality
ggdensity(all_features$GC_content)
ggqqplot(all_features$GC_content)
shapiro.test(all_features$GC_content) ## p<0.05 = Differ from normality, so non-parametric test

####Kruskal-Wallis for non-parametric data
kruskal.test(GC_content~Clades, data = all_features) ##significant differences based on Species
####Dunn's Post-hoc test for mutiple comparisons
dunnTest(GC_content~Clades, data = all_features,
         method="bonferroni")

