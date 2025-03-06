############################################
### Ana Cuesta CAZyme analyses AAB paper ###
############################################
##Clean your environment
rm(list=ls())

library(vegan)
library(devtools)
library(tibble)
library(dplyr)
library(effectsize)
library(pairwiseAdonis)
library(phyloseq)
library(phyloseqGraphTest)
library(readxl)
library(beeswarm)
library("ggbeeswarm")
library(ggpubr)
library(ggrepel)
library(ggVennDiagram)
library("venn")
library("ggvenn")
library("VennDiagram")
library("writexl")

## Loading in the data and formatting it ##
#=========================================#
basepath = "/Users/fmg115/Desktop/CAZYmes"

## Metadata ## with acidophilic and human diseases
metadata = read_xlsx("/Users/fmg115/Desktop/CAZYmes/metadata_new.xlsx")
metadata$isolation = as.factor(metadata$isolation)
metadata <- metadata[!metadata$sample == "Granulibacter_bethesdensis",]

## CAZYme family profiles
family_profiles_abun_wide = readRDS("/Users/fmg115/Desktop/CAZYmes/family_profiles_abun_wide.rds")
CAZymeprofiles <- as.data.frame(family_profiles_abun_wide)
write_xlsx(family_profiles_abun_wide, "/Users/fmg115/Desktop/CAZYmes/CAZYMESFINAL.xlsx")

## Check the data
# Print the first few rows of the data frame
head(family_profiles_abun_wide)
# Inspect the structure of the data frame
str(family_profiles_abun_wide)

# Clean and prepare the data
family_profiles_numeric <- family_profiles_abun_wide[, sapply(family_profiles_abun_wide, is.numeric)]

# Remove  outgroups 
family_profiles_abun_wide <- family_profiles_abun_wide[!family_profiles_abun_wide$bin == "Granulibacter_bethesdensis",]
family_profiles_abun_wide <- family_profiles_abun_wide[!family_profiles_abun_wide$bin == "Acidocella_aromatica",]
family_profiles_abun_wide <- family_profiles_abun_wide[!family_profiles_abun_wide$bin == "Roseomonas_mucosa",]



#################################
### Calculate alpha diversity ###
#################################

### Observed richness ###

obs_df_cazyme = family_profiles_abun_wide %>%
  mutate_if(is.numeric, ~1 * (. != 0)) %>%
  estimateR(.) %>%
  as.data.frame() %>%
  dplyr::slice(.data=., 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "bin") %>%
  left_join(metadata, by = c("bin" = "sample")) %>%
  column_to_rownames(var = "bin") %>%
  rename(obs_cazyme = S.obs)

## Anova ##
obs_cazyme.aov = aov(obs_cazyme ~ isolation, data = obs_df_cazyme)
summary(obs_cazyme.aov)
TukeyHSD(obs_cazyme.aov, which = "isolation")
cohens_f(obs_cazyme.aov)

## Plot Observed richness
o <- ggplot(obs_df_cazyme, aes(x=isolation, y=obs_cazyme, fill=isolation)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values =  c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) + scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16) +
  xlab("Isolation source") + ylab("Obverved richness") + ggtitle("Observed richness of CAZyme families") + 
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 18), axis.text = element_text(size = 12))
o
o <- o +theme_bw()
o <- o + theme(legend.position = "none")
print(o)



### Shannon richness ###

shannon_df_cazyme = as.data.frame(vegan::diversity(family_profiles_abun_wide, index = "shannon")) %>%
  rename(shannon_cazyme = `vegan::diversity(family_profiles_abun_wide, index = "shannon")`) %>%
  rownames_to_column() %>%
  left_join(metadata, by = c("bin" = "sample")) %>%
  column_to_rownames()

## Anova ##
shannon_cazyme.aov = aov(Shannon_Diversity ~ isolation, data = family_profiles_abun_wide)
summary(shannon_cazyme.aov)
TukeyHSD(shannon_cazyme.aov, which = "isolation")
cohens_f(shannon_cazyme.aov)

## Plot Shannon diversity
s <- ggplot(shannon_df_cazyme, aes(x=isolation, y=Shannon_Diversity, fill=isolation)) + geom_boxplot(alpha=.5) +
  scale_fill_manual(values =  c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) + scale_color_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD", "#FF40FF")) +
  geom_jitter(aes(colour = Isolation_source, alpha = .9), shape=16) +
  xlab("Isolation source") + ylab("Shannon diversity") + ggtitle("Shannon diversity of CAZyme families") + 
  theme(legend.position = "none", plot.title = element_text(size = 22, face = "bold", hjust = .5), axis.title = element_text(size = 18), axis.text = element_text(size = 12))
s
s <- s +theme_bw()
s <- s + theme(legend.position = "none")
print(s)


#### Alpha diversity analysis end ####

################################
### Calculate beta diversity ###
################################
set.seed(123)

nmds = metaMDS(family_profiles_abun_wide, distance = "bray")
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds)$sites) #extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
#add columns to data frame 
data.scores$bin = family_profiles_abun_wide$bin
data.scores$isolation <- metadata$isolation
head(data.scores)

xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes(colour = isolation))+ 
  stat_ellipse(aes(colour = isolation))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "isolation", y = "NMDS2")  + 
  scale_colour_manual(values = c("#E6A401","#C4C2BB","#73C221", "#9933CD")) 

xx

### Beta diversity statistics
adonis.res = adonis2(family_profiles_abun_wide ~ isolation, data = metadata, permutations = 9999, method = "bray", by = "terms")
adonis.res

adonis_pairwise.res = pairwise.adonis2(family_profiles_abun_wide ~ isolation , data = metadata, method = "bray", nperm = 9999)
adonis_pairwise.res
