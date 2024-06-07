library("tidyverse")
library("readxl")
library(ggplot2)
library(readxl)
library(FSA)

# COMPLEMENTARY PLOTS
# Set paths
basepath = "/Users/thibaultbret/"
excelpath = paste(basepath,"Documents/Work/AAB Apr2024 BoxplotsOnly.xlsx",sep="")

# Read data
all_features <- read_excel(excelpath, sheet = "AAB genomes", range = cell_cols("A:Z"))
# Remove unneeded rows
all_features <- subset(all_features, Isolation_source != "NA" & Isolation_source != "Ferment/Plant" & Isolation_source != "Human disease" & Isolation_source != "Sourdough")
all_features <- subset(all_features, Isolated_from != "Not verified" & Isolated_from != "NA")
# Set categories for Isolation_source row
all_features$Isolation_source <- factor(all_features$Isolation_source, levels=c('Insect','Fly','Plant','Ferment','Ferment/Industrial'))
# Abbreviate species names for clearer labels except those that end with "sp."
abbreviate <- function(col) paste(substr(gsub( " .*$", "", col),start=0,stop=1), sub("^\\S+\\s+", '',  col), sep =  ". ")
all_features <- all_features %>% mutate(Species = ifelse(!grepl(' sp\\.', Species), abbreviate(Species), Species))


# GC content
kruskal.test(GC_content ~ Isolation_source, data = all_features)
dunnTest(GC_content ~ Isolation_source, data = all_features)

# Genome size
kruskal.test(Size ~ Isolation_source, data = all_features)
dunnTest(Size ~ Isolation_source, data = all_features)

# CDS prop
kruskal.test(CDS_prop ~ Isolation_source, data = all_features)
dunnTest(CDS_prop ~ Isolation_source, data = all_features)

# Gene count
kruskal.test(`Total gene count` ~ Isolation_source, data = all_features)
dunnTest(`Total gene count` ~ Isolation_source, data = all_features)

# non-coding gene count
kruskal.test(as.numeric(`Nb Non-coding genes`) ~ Isolation_source, data = all_features)
dunnTest(as.numeric(`Nb Non-coding genes`) ~ Isolation_source, data = all_features)

# pseudogenes
kruskal.test(`Nb pseudogenes` ~ Isolation_source, data = all_features)
dunnTest(as.numeric(`Nb pseudogenes`) ~ Isolation_source, data = all_features)


# Pearson's correlation
cor.test(all_features$GC_content, all_features$Size, method = "pearson")

