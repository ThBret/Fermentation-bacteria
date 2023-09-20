library(ape)
library(phangorn)

# Setup
directory_path <- "/Users/thibaultbret/Documents/Work/Individual trees/"

pfams_dict <- c("PF00005"="ABC transporter",
  "PF00196"="Bacterial regulatory proteins, luxR family",
  "PF00923"="Transaldolase/Fructose-6-phosphate aldolase",
  "PF02887"="Pyruvate kinase, alpha/beta domain",
  "PF08240"="Alcohol dehydrogenase GroES-like domain",
  "PF00118"="TCP-1/cpn60/GroEL chaperonin family",
  "PF00330"="Aconitase family (aconitate hydratase)",
  "PF00958"="GMP synthase C terminal domain",
  "PF03070"="TENA/THI-4/PQQC family",
  "PF13243"="Squalene-hopene cyclase C-terminal domain",
  "PF00171"="Aldehyde dehydrogenase family",
  "PF00465"="Iron-containing alcohol dehydrogenase",
  "PF01161"="Phosphatidylethanolamine-binding protein",
  "PF08042"="PqqA family",
  "PF13360"="PQQ-like domain")

pfams = c("PF00005","PF00196","PF00923","PF02887","PF08240","PF00118","PF00330","PF00958",
          "PF03070","PF13243","PF00171","PF00465","PF01161","PF08042","PF13360")

Ferments = c("Acetobacter_cerevisiae","Acetobacter_estunensis","Acetobacter_aceti","Gluconobacter_oxydans","Gluconobacter_japonicus","Acidomonas_methanolica","Acetobacter_fabarum")
Plants = c("Saccharibacter_floricola","Gluconobacter_kanchanaburiensis","Gluconobacter_cadivus","Gluconobacter_albidus","Gluconobacter_aidae","Gluconobacter_frateurii","Gluconobacter_thailandicus","Asaia_platycodi","Kozakia_baliensis","Neoasaia_chiangmaiensis","Gluconacetobacter_diazotrophicus","Gluconacetobacter_liquefaciens")
Insects = c("Bombella_apis","Saccharibacter_sp","Parasaccharibacter_apium","Bombella_intestini","Bombella_favorum","Bombella_mellum","Gluconobacter_sphaericus","Commensalibacter_intestini","Commensalibacter_sp","Commensalibacter_papalotli","Acetobacter_oryzifermentans","Acetobacter_pomorum","Oecophyllibacter_saccharovorans")

# Load data
tree_files <- paste(pfams, "-proteins.fa.treefile", sep = "")

pfam_trees <- lapply(tree_files, function(file) read.tree(file.path(directory_path, file)))

# Open a PDF device for plotting
pdf("/Users/thibaultbret/trees.pdf", width = 10, height = 10)  # Adjust width and height as needed

# Plot trees
for (i in 1:length(pfam_trees)) {
  # Create a vector of colours corresponding to the label groups
  label_colours <- rep("black", length(pfam_trees[[i]]$tip.label))
  label_colours[pfam_trees[[i]]$tip.label %in% Ferments] <- "#E27842"
  label_colours[pfam_trees[[i]]$tip.label %in% Insects] <- "#E6A401"
  label_colours[pfam_trees[[i]]$tip.label %in% Plants] <- "#73C221"

  # Plot the tree
  plotBS(pfam_trees[[i]], edge.color = "black", tip.color = label_colours, main = paste(pfams[i]," (",pfams_dict[pfams[i]],")", sep = ""))
}

# Close the PDF device
dev.off()



# PLOT INDIVIDUAL TREE
pdf("/Users/thibaultbret/trees.pdf", width = 10, height = 10)  # Adjust width and height as needed

plot_indiv_tree <- function(input_tree, id){
  tree_file <- read.tree(input_tree)
  label_colours <- rep("black", length(tree_file$tip.label))
  label_colours[tree_file$tip.label %in% Ferments] <- "#E27842"
  label_colours[tree_file$tip.label %in% Insects] <- "#E6A401"
  label_colours[tree_file$tip.label %in% Plants] <- "#73C221"
  # Plot the tree
  if(missing(id)){
  plotBS(tree_file, edge.color = "black", tip.color = label_colours, main = "Single-copy core gene tree", sep = "")
  } else {
  plotBS(tree_file, edge.color = "black", tip.color = label_colours, main = paste(id," (",pfams_dict[id],")", sep = ""))
  }
}

#plot_indiv_tree("BIS-PF00005-proteins.fa.treefile","PF00005")
plot_indiv_tree("/Users/thibaultbret/ALIGN-bp.treefile")

dev.off()


