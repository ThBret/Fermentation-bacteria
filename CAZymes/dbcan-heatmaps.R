setwd("/Users/thibaultbret/dbcan-output")
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(phytools)
library(grid)

#### Read in data ####
#====================#
## Metadata ##
metadata = read.delim("/Users/thibaultbret/isolationgroups.txt", header = FALSE, col.names = c("bin","isolation_source"))

## Cazyme family ##
family_profiles_abun_wide = readRDS("/Users/thibaultbret/family_profiles_abun_wide.rds")

## Cazyme substrate ##
substrate_profiles_abun_wide = readRDS("/Users/thibaultbret/substrate_profiles_abun_wide.rds")

## Tree ##
tree <- read.tree("/Users/thibaultbret/Documents/Work/sccg-tree-noSRR.treefile")
tree$tip.label <- gsub("_", " ", tree$tip.label)


#### Cazyme family heatmap ####
#=============================#

### Preparing datafarmes for heatmap ###
family_profiles_abun_wide <- family_profiles_abun_wide[!grepl("^SRR", family_profiles_abun_wide$bin), ]
family_profiles_abun_wide_df = data.frame(family_profiles_abun_wide, row.names = "bin")
# Reformat species names
rownames(family_profiles_abun_wide_df) <- gsub("_", " ", row.names(family_profiles_abun_wide_df))
# Reorder family_profiles_abun_wide_df rows based on phylogeny
family_profiles_abun_wide_df <- family_profiles_abun_wide_df[rev(tree$tip.label), ]

# Annotation rows
annotations.rows = metadata %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "bin") %>%
  arrange(desc(isolation_source))
rownames(annotations.rows) <- gsub("_", " ", rownames(annotations.rows))
annotations.rows$isolation_source = factor(annotations.rows$isolation_source, levels = c("Ferment", "Plant", "Insect","Fly","Industrial_ferment","Acidophilic","Human_disease"))

colours = list(
  isolation_source = c(Fly = "#CCCCCC", Industrial_ferment = "#FF00FF", Insect = "#E6A401", Ferment = "#9933CC", Plant = "#73C221", Acidophilic = "#299596", Human_disease = "#000000"))

## log10 transfrm data ##
family_profiles_abun_wide_log10 = log10(family_profiles_abun_wide_df)
family_profiles_abun_wide_log10[family_profiles_abun_wide_log10 == "-Inf"] = -(1-as.numeric(range(family_profiles_abun_wide_log10[family_profiles_abun_wide_log10 != "-Inf"])[1]))# removing "inf" spawned from log10()

# Reorder annotations.rows based on row names of family_profiles_abun_wide_log10
annotations.rows <- annotations.rows[rev(tree$tip.label), , drop = FALSE]


# Plot heat map
pdf(file = "/Users/thibaultbret/cazyme-family-profiles-4k.pdf", width = 15, height = 20)
pheatmap(as.matrix(family_profiles_abun_wide_log10), 
         fontsize = 10,
         fontsize_col = 6,
         annotation_legend = FALSE, 
         show_colnames = TRUE,
         #annotation_col = annotations.cols,
         annotation_row = annotations.rows,
         annotation_colors = colours,
         border_color = FALSE,
         cutree_cols = 4,
         cluster_rows = FALSE,
         #annotation_names_row = FALSE,
         #row_split = annotations.rows$isolation_source,
         #legend_labels = c("Absence", "Very low presence", "Low presence", "Moderate presence", "High presence"),
         legend_breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)
         #cluster_row_slices = FALSE)
)
dev.off()




#### Cazyme substrate heatmap ####
#===============================#

### Preparing datafarmes for heatmap ###
substrate_profiles_abun_wide <- substrate_profiles_abun_wide[!grepl("^SRR", substrate_profiles_abun_wide$bin), ]
substrate_profiles_abun_wide_df = data.frame(substrate_profiles_abun_wide, row.names = "bin")
# Reformat species names
rownames(substrate_profiles_abun_wide_df) <- gsub("_", " ", row.names(substrate_profiles_abun_wide_df))
# Reorder substrate_profiles_abun_wide_df rows based on phylogeny
substrate_profiles_abun_wide_df <- substrate_profiles_abun_wide_df[rev(tree$tip.label), ]

# Filter out substrates with few hits (less than 20)
column_sums <- colSums(substrate_profiles_abun_wide_df, na.rm = TRUE)
columns_less_than_20 <- names(column_sums[column_sums < 20])
filtered_substrate_profiles_abun_wide_df <- substrate_profiles_abun_wide_df[, !names(substrate_profiles_abun_wide_df) %in% columns_less_than_20]

# Filter out rows with no data
row_sums <- rowSums(filtered_substrate_profiles_abun_wide_df, na.rm = TRUE)
non_zero_rows <- row_sums != 0
filtered_substrate_profiles_abun_wide_df <- filtered_substrate_profiles_abun_wide_df[non_zero_rows, ]

## log10 transform data ##
substrate_profiles_abun_wide_log10 = log10(filtered_substrate_profiles_abun_wide_df)
substrate_profiles_abun_wide_log10[substrate_profiles_abun_wide_log10 == -Inf] <- -(1 - min(substrate_profiles_abun_wide_log10[substrate_profiles_abun_wide_log10 != -Inf], na.rm = TRUE)) # removing "inf" spawned from log10()

# Remove annotations for rows that were filtered out 
filtered_annotations_rows <- annotations.rows[rownames(substrate_profiles_abun_wide_log10), , drop = FALSE]

# Plot heat map
pdf(file = "/Users/thibaultbret/cazyme-substrate-profiles.pdf", width = 15, height = 20)
pheatmap(as.matrix(substrate_profiles_abun_wide_log10), 
         fontsize = 10,
         annotation_legend = FALSE, 
         show_colnames = TRUE,
         #annotation_col = annotations.cols,
         annotation_row = filtered_annotations_rows,
         annotation_colors = colours,
         border_color = FALSE,
         cluster_rows = FALSE,
         #cutree_cols = 9,
         #annotation_names_row = FALSE,
         #legend_labels = c("Absence", "Very low presence", "Low presence", "Moderate presence", "High presence"),
         legend_breaks = c(-1, -0.5, 0, 0.5, 1, 1.5)
         #cluster_row_slices = FALSE)
)

dev.off()
