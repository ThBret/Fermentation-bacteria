library("tidyverse")
library("readxl")
library(ggplot2)
library(ggrepel)
library(readxl)

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


### Scatter plot GC content / CDS prop vs genome size with ellipses ###
# Read data (Acidophilic) and combine with AAB data
Acidophilic = subset(read_excel(excelpath, sheet = "Acidophilic bacteria", range = cell_cols("A:J")), select = c("Species","Assembly Accession","Isolation_source","Size","GC_content","CDS_prop","Isolated_from"))
all_features_with_acidophilic = bind_rows(subset(all_features, select = c("Species","Assembly Accession","Isolation_source","Size","GC_content","CDS_prop","Isolated_from")), Acidophilic)
all_features_with_acidophilic <- all_features_with_acidophilic %>% mutate(Species = ifelse(!grepl(' sp\\.', Species), abbreviate(Species), Species))

# Set palette                            acidophilic / ferment / industrial / fly / insect / plant
Isolation_source_palette_acidophilic <- c("#299596","#9933CC","#FF00FF","#CCCCCC",'#E6A401',"#73C221")
# Extract labels
Isolated_from_labels <- all_features_with_acidophilic %>% pull(Isolated_from)

### GC content vs genome size with labels (isolated from) ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=GC_content, color=Isolation_source)) +
  geom_point(size = 2, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  geom_text_repel(aes(label = Isolated_from, colour = Isolation_source), size = 5, max.overlaps = 16) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("GC content (in %)") + labs(color = "Isolation source") +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) +
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-GC-isolationlabels.png", p, width = 25, height = 15)

### CDS prop vs genome size with labels (isolated from)  ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=CDS_prop, color=Isolation_source)) +
  geom_point(size = 3, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  geom_text_repel(aes(label = Isolated_from, colour = Isolation_source), size = 5, max.overlaps = 16) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("CDS proportion (in %)") + labs(color = "Isolation source") + guides(fill = none) +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) +
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-CDS-isolationlabels.png", p, width = 25, height = 15)

### GC content vs genome size with labels (species names) ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=GC_content, color=Isolation_source)) +
  geom_point(size = 2, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 16) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("GC content (in %)") + labs(color = "Isolation source") +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) +
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-GC-specieslabels.png", p, width = 25, height = 15)

### CDS prop vs genome size with labels (speces names)  ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=CDS_prop, color=Isolation_source)) +
  geom_point(size = 3, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 16) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("CDS proportion (in %)") + labs(color = "Isolation source") + guides(fill = none) +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) +
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-CDS-specieslabels.png", p, width = 25, height = 15)

### GC content vs genome size without labels ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=GC_content, color=Isolation_source)) +
  geom_point(size = 2, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("Guanine-Cytosine content (in %)") + labs(color = "Isolation source") +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) + 
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-GC-nolabels.png", p, width = 25, height = 15)

### CDS prop vs genome size without labels  ###
p <- ggplot(all_features_with_acidophilic, aes(x=Size/1000000, y=CDS_prop, color=Isolation_source)) +
  geom_point(size = 3, alpha = 0.5) + stat_ellipse(lwd = 2, alpha = 0.5) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("CDS proportion (in %)") + labs(color = "Isolation source") + guides(fill = none) +
  guides(fill = guide_legend(override.aes = list(color = Isolation_source_palette_acidophilic))) +
  scale_fill_manual(values = Isolation_source_palette_acidophilic) + scale_color_manual(values = Isolation_source_palette_acidophilic) + theme_light() +
  theme(plot.title = element_text(family = "Times New Roman", size = 40, face = "bold", hjust = .5), axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30), legend.key.size = unit(1.5, 'cm'),
        legend.key.width = unit(3, 'cm'), legend.title = element_text(size=20, face = "bold", family = "Times New Roman"), legend.text = element_text(size=25, family = "Times New Roman", face = "bold"), legend.position = "inside", legend.position.inside = c(.85,.2)) + 
  guides(fill=guide_legend(title="Isolation source"))

ggsave("/Users/thibaultbret/scatterplot-CDS-nolabels.png", p, width = 25, height = 15)







# Set palette                     insect / fly / plant / ferment / industrial
Isolation_source_palette <- c('#E6A401',"#CCCCCC","#73C221","#9933CC","#FF00FF")

### Boxplot GC content with labels (isolated from) ###
Isolated_from_labels <- all_features %>% pull(Isolated_from)
pos <- position_jitter(width = 0.3, seed = 2)
p <- ggplot(all_features, aes(x=Isolation_source, y=GC_content, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"gcplot-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot GC content with labels (species names) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=GC_content, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"gcplot-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot GC content without labels ###
p <- ggplot(all_features, aes(x=Isolation_source, y=GC_content, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"gcplot-nolabels.png",sep=""), p, width = 18, height = 14)



### Boxplot CSD length with labels (isolated from) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_length/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDSlen-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot CDS length with labels (species names) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_length/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDSlen-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot CDS length without labels ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_length/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDSlen-nolabels.png",sep=""), p, width = 18, height = 14)



# Remove NA values
all_features_no_NA <- all_features %>% drop_na(`Nb Non-coding genes`)
Isolated_from_labels_no_NA <- all_features_no_NA %>% pull(Isolated_from)
### Boxplot Nb non-coding genes with labels (isolated from) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb Non-coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels_no_NA, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"noncodinggenes-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb non-coding genes with labels (species names) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb Non-coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"noncodinggenes-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb non-coding genes without labels ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb Non-coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"noncodinggenes-nolabels.png",sep=""), p, width = 18, height = 14)



# Remove NA values
all_features_no_NA <- all_features %>% drop_na(`Nb protein coding genes`)
Isolated_from_labels_no_NA <- all_features_no_NA %>% pull(Isolated_from)
### Boxplot Nb protein-coding genes with labels (isolated from) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb protein coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels_no_NA, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"proteincodinggenes-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb protein-coding genes with labels (species names) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb protein coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"proteincodinggenes-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb protein-coding genes without labels ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb protein coding genes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"proteincodinggenes-nolabels.png",sep=""), p, width = 18, height = 14)



# Remove NA values
all_features_no_NA <- all_features %>% drop_na(`Nb pseudogenes`)
Isolated_from_labels_no_NA <- all_features_no_NA %>% pull(Isolated_from)
### Boxplot Nb pseudogenes with labels (isolated from) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb pseudogenes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels_no_NA, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"pseudogenes-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb pseudogenes with labels (species names) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb pseudogenes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"pseudogenes-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Nb pseudogenes without labels ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(`Nb pseudogenes`), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"pseudogenes-nolabels.png",sep=""), p, width = 18, height = 14)



### Boxplot Genome size with labels (isolated from) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genomesize-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Genome size with labels (species names) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genomesize-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Genome size without labels ###
p <- ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genomesize-nolabels.png",sep=""), p, width = 18, height = 14)



### Boxplot CDS prop with labels (isolated from) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDS-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot CDS prop with labels (species names) ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDS-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot CDS prop without labels ###
p <- ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"CDS-nolabels.png",sep=""), p, width = 18, height = 14)



### Boxplot Total gene count with labels (isolated from) ###
# Remove NA values
all_features_no_NA <- all_features %>% drop_na(`Total gene count`)
Isolated_from_labels_no_NA <- all_features_no_NA %>% pull(Isolated_from)

p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=`Total gene count`, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Isolated_from_labels_no_NA, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genecount-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Total gene count with labels (species names) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=`Total gene count`, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 8, max.overlaps = 8,  position=pos) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genecount-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot Total gene count without labels ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=`Total gene count`, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_classic() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"genecount-nolabels.png",sep=""), p, width = 18, height = 14)








### Boxplot growth temperature with labels (isolated from) ###
# Remove NA values
all_features[all_features == "NA"] <- NA
all_features_no_NA <- all_features %>% drop_na(`Ta`)
Isolated_from_labels <- all_features_no_NA %>% pull(Isolated_from)

p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(Ta), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) +
  geom_text_repel(aes(label = Isolated_from_labels, colour = Isolation_source), size = 3, max.overlaps = 8,  position=pos) +
  xlab("Isolation source") + ylab("Temperature (in °C)") + ggtitle("Growth temperature of Acetic Acid Bacteria genomes") + theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 20), axis.text = element_text(size = 15))

ggsave(paste(basepath,"tempplot-isolationlabels.png",sep=""), p, width = 18, height = 14)

### Boxplot growth temperature with labels (species names) ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(Ta), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) +
  geom_text_repel(aes(label = Species, colour = Isolation_source), size = 3, max.overlaps = 8,  position=pos) +
  xlab("Isolation source") + ylab("Temperature (in °C)") + ggtitle("Growth temperature of Acetic Acid Bacteria genomes") + theme_light() +
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 18), axis.text = element_text(size = 15))

ggsave(paste(basepath,"tempplot-specieslabels.png",sep=""), p, width = 18, height = 14)

### Boxplot growth temperature without labels ###
p <- ggplot(all_features_no_NA, aes(x=Isolation_source, y=as.numeric(Ta), fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values = Isolation_source_palette) + scale_color_manual(values = Isolation_source_palette) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=pos) + theme_light() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), legend.position = "none", axis.title = element_text(size = 35, family = "Times New Roman", face = "bold"), axis.text = element_text(size = 30))

ggsave(paste(basepath,"tempplot-nolabels.png",sep=""), p, width = 18, height = 14)
