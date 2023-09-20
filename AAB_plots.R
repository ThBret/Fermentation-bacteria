library("tidyverse")
library("readxl")
library(ggplot2)

#phylo trees
library(phytools)
pfam_tree <- ape::read.tree('/Users/thibaultbret/pfam-phylogenomic-tree.txt')
sccg_tree <- ape::read.tree('/Users/thibaultbret/ALIGN.treefile')
cophylo_tree <- cophylo(pfam_tree, sccg_tree)
## plot co-phylogenies
plot(cophylo_tree,link.type="curved",link.lty="solid",link.col=make.transparent("red",0.25))

#read data
library(readxl)
all_features <- read_excel("complete_df_ACM.xlsx")
head(all_features)

#Remove NAs & industrial
all_features <- subset(all_features, Isolation_source != "NA" & Isolation_source != "Industrial")
all_features$Isolation_source <- factor(all_features$Isolation_source, levels=c('Insect','Plant','Ferment'))

#Summary statistics
summary = all_features %>%
  group_by(Isolation_source) %>%
  summarise("Mean GC content (in %)" = mean(GC_content), "Mean CDS proportion (in %)" = mean(CDS_prop), "Mean genome size (in Mb)" = mean(Size)/1000000)

view(summary)
#scatter plot GC content (no trend lines)
ggplot(all_features, aes(x=Size/1000000, y=GC_content, color=Isolation_source)) +
  geom_point() + stat_ellipse(geom = "polygon", aes(fill = Isolation_source), alpha = 0.25) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("GC content (in %)") + labs(color = "Isolation source") + guides(fill = FALSE)

#scatter plot CDS proportion (no trend lines)
ggplot(all_features, aes(x=Size/1000000, y=CDS_prop, color=Isolation_source)) +
  geom_point() + stat_ellipse() + stat_ellipse(geom = "polygon", aes(fill = Isolation_source), alpha = 0.25) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000-.5, max(all_features$Size)/1000000)) +
  ylab("Proportion of CDS (in %)") + labs(color = "Isolation source") + guides(fill = FALSE)

#scatter plot GC content (with trend lines)
ggplot(all_features, aes(x=Size/1000000, y=GC_content, color=Isolation_source)) +
  geom_point(alpha = .6) +
  geom_smooth(method=lm, formula = y ~ x, aes(fill=Isolation_source)) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000, max(all_features$Size)/1000000)) +
  ylab("GC content (in %)") + labs(color = "Isolation source") + guides(fill = FALSE)

#scatter plot CDS proportion (with trend lines)
ggplot(all_features, aes(x=Size/1000000, y=CDS_prop, color=Isolation_source)) +
  geom_point(alpha = .6) +
  geom_smooth(method=lm, formula = y ~ x, aes(fill=Isolation_source)) +
  scale_x_continuous(name="Genome size (in Mb)", limits=c(min(all_features$Size)/1000000, max(all_features$Size)/1000000)) +
  ylab("Proportion of CDS (in %)") + labs(color = "Isolation source") + guides(fill = FALSE)

#boxplot GC content
ggplot(all_features, aes(x=Isolation_source, y=GC_content, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=position_jitter(0.2)) +
  xlab("Isolation source") + ylab("GC content (in %)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842")) + scale_color_manual(values = c('#E6A401',"#73C221","#E27842"))

ggplot(all_features, aes(x=Isolation_source, y=GC_content, fill=Isolation_source)) + geom_boxplot() +
  xlab("Isolation source") + ylab("GC content (in %)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842"))

#boxplot CDS proportion
ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=position_jitter(0.2)) +
  xlab("Isolation source") + ylab("Proportion of CDS (in %)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842")) + scale_color_manual(values = c('#E6A401',"#73C221","#E27842"))

ggplot(all_features, aes(x=Isolation_source, y=CDS_prop, fill=Isolation_source)) + geom_boxplot() +
  xlab("Isolation source") + ylab("Proportion of CDS (in %)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842"))

#boxplot Genome size
ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot(alpha=.5) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16, position=position_jitter(0.2)) +
  xlab("Isolation source") + ylab("Genome size (in Mb)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842")) + scale_color_manual(values = c('#E6A401',"#73C221","#E27842"))

ggplot(all_features, aes(x=Isolation_source, y=Size/1000000, fill=Isolation_source)) + geom_boxplot() +
  xlab("Isolation source") + ylab("Genome size (in Mb)") + theme(legend.position = "none") +
  scale_fill_manual(values = c('#E6A401',"#73C221","#E27842"))

